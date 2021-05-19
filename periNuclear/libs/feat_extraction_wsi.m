function feat_extraction_wsi(folder_matpatches,folder_pyepistroma,folder_matcellmask,folder_savepath,quality)
placeholder = 'placeholder';
imgList=dir([folder_matpatches '*.mat']);
indx = randperm(numel(imgList));
numFiles=length(imgList);
imgList = {imgList(:).name};
qualitysaturationLim = quality.saturationLim;
qualityredChannelLim = quality.redChannelLim;
qualityblurLimit = quality.blurLimit;
qualityimgarea = quality.imgarea;
parfor nn=1:numFiles
    n= indx(nn);
    [~,imgName]=fileparts(imgList{n});
    % ToDo: Check if the nuclei and ES are already being processed, to be
    % used with another feature set
    imgFile=[folder_matpatches imgList{n}];
    outputFolder=[folder_savepath 'dataset_output/' imgName];
    maskFolder=[outputFolder '/png_binmask/png_cellmask/'];
    maskEFolder=[outputFolder '/png_binmask/png_cellepimask/'];
    maskSFolder=[outputFolder '/png_binmask/png_cellstromask/'];
    maskBFolder=[outputFolder '/png_binmask/png_cellbundmask/'];
    ESmaskFolder=[outputFolder '/png_binmask/png_epithmask/'];
    SSmaskFolder=[outputFolder '/png_binmask/png_stromask/'];
    EBmaskFolder=[outputFolder '/png_binmask/png_boundaryepistromask/'];
    featlocFolder=[outputFolder '/mat_perifeatures/'];
    disp([outputFolder])
    mkdir(outputFolder);
    mkdir(maskFolder);
    mkdir(maskEFolder);
    mkdir(maskSFolder);
    mkdir(maskBFolder);
    mkdir(ESmaskFolder);
    mkdir(SSmaskFolder);
    mkdir(EBmaskFolder);
    mkdir(featlocFolder);
    % open .mat file
    imgFile_load = load(imgFile);
    numtilestruct = length(imgFile_load.tileStruct);
    indx2 = randperm(numel(imgFile_load.tileStruct));
    fprintf('Processing file %s and tiles %d\n',imgName,numtilestruct);
    try
        curTile_ESmask_struct=h5read([folder_pyepistroma imgName '_ESmask'],'/mask');
    catch
        disp('disrupted or missing Epistroma Mask');
        filePh_error = fopen([folder_savepath 'errorlist_patches.txt'],'a');
        fprintf(filePh_error,'Disrupted epistroma mask %s with tiles %d\n',imgName,numtilestruct);
        fclose(filePh_error);
        % ToDo: delete the created folder
        rmdir(outputFolder,'s')
        continue;
    end
    curTile_Nmask_struct=h5read([folder_matcellmask imgName '_mask'],'/mask');
    for i=1:numtilestruct
        featFile=sprintf('%s/%s_%d.mat',featlocFolder,imgName,indx2(i));
        maskFile=sprintf('%s/%s_%d.png',maskFolder,imgName,indx2(i));
        maskEFile=sprintf('%s/%s_%d.png',maskEFolder,imgName,indx2(i));
        maskSFile=sprintf('%s/%s_%d.png',maskSFolder,imgName,indx2(i));
        maskBFile=sprintf('%s/%s_%d.png',maskBFolder,imgName,indx2(i));
        ESmaskFile=sprintf('%s/%s_%d.png',ESmaskFolder,imgName,indx2(i));
        SSmaskFile=sprintf('%s/%s_%d.png',SSmaskFolder,imgName,indx2(i));
        EBmaskFile=sprintf('%s/%s_%d.png',EBmaskFolder,imgName,indx2(i));
        curTile=imgFile_load.tileStruct(indx2(i)).data;
        if (getSaturationMetric(curTile)>qualitysaturationLim && ...
                getRedMetric(curTile)>qualityredChannelLim && ...
                blurMetric(curTile)>qualityblurLimit && ...
                getAreaTissue(curTile)>qualityimgarea)
            if exist(featFile,'file')~=2
                parsave(featFile, placeholder);
                % Check if Cell mask and Epistroma mask have same quantity
                if size(curTile_ESmask_struct) == size(curTile_Nmask_struct)
                    fprintf('Processing tile %s_%d\n',imgName,indx2(i));
                    curTile_ESmask = curTile_ESmask_struct(:,:,indx2(i))'; % For a single mask
                    curTile_Nmask = curTile_Nmask_struct(:,:,indx2(i))'; % For a single mask
                    [nucleiCentroids,nucleiCentroids_stro,nucleiCentroids_epi, nucleiCentroids_bund,...
                        feat_perinuclei, feat_epi_perinuclei, feat_stro_perinuclei, feat_bund_perinuclei] = get_periNuclear(curTile,curTile_ESmask,curTile_Nmask,maskFile,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile);
                    parsave_cellfeat_peri(featFile,nucleiCentroids,nucleiCentroids_stro,nucleiCentroids_epi,nucleiCentroids_bund,...
                        feat_perinuclei, feat_epi_perinuclei, feat_stro_perinuclei, feat_bund_perinuclei);
                else
                    filePh_error = fopen([folder_savepath 'errorlist_patches.txt'],'a');
                    fprintf(filePh_error,'Mismatch cell and epistroma mask %s_%d\n',imgName,indx2(i));
                    fclose(filePh_error);
                    continue
                end
            else
                %disp('mat tile already exists')
                continue
            end
        end
    end
end
%ToDo: include method to delete the 1Kb files
parfor n=1:numFiles
    [~,imgName]=fileparts(imgList{n});
    % ToDo: Check if the nuclei and ES are already being processed, to be
    % used with another feature set
    %imgFile=[folder_matpatches imgList{n}];
    outputFolder=[folder_savepath 'dataset_output/' imgName];
    featlocFolder=[outputFolder '/mat_perifeatures/*.mat'];
    files = dir(featlocFolder);
    for ii = 1:length(files)
        if files(ii).bytes==183 % file with 'placeholder' word
            delete(fullfile(files(ii).folder, files(ii).name))
        end
    end
end
fprintf('Done!\n');
% Save the geenral feature description
 description_periNuclei = {'mean value of periNuclei area', 'median of periNuclei area', 'variance of periNuclei area', 'standard deviation of periNuclei area', ...
    'mean value of periNuclei area to nuclei area ratio', 'mean value of periNuclei area to nuclei area ratio', 'mean value of periNuclei area to nuclei area ratio', 'mean value of periNuclei area to nuclei area ratio',...
    'Haralick:mean intensity angular_second_moment_energy','Haralick:standard deviation intensity angular_second_moment_energy',...
    'Haralick:mean intensity contrast_energy','Haralick:standard deviation intensity contrast_energy',...
    'Haralick:mean intensity correlation','Haralick:standard deviation intensity correlation',...
    'Haralick:mean intensity contrast_var','Haralick:standard deviation intensity contrast_var',...
    'Haralick:mean intensity inverse_difference_moment_homogeneity','Haralick:standard deviation intensity inverse_difference_moment_homogeneity',...
    'Haralick:mean intensity sum_average','Haralick:standard deviation intensity sum_average',...
    'Haralick:mean intensity sum_variance_aprox','Haralick:standard deviation intensity sum_variance_aprox',...
    'Haralick:mean intensity sum_entropy','Haralick:standard deviation intensity sum_entropy',...
    'Haralick:mean intensity entropy_cut_out_zeros','Haralick:standard deviation intensity entropy_cut_out_zeros',...
    'Haralick:mean intensity difference_variane,','Haralick:standard deviation intensity difference_variane,',...
    'Haralick:mean intensity difference_entropy','Haralick:standard deviation intensity difference_entropy',...
    'Haralick:mean intensity information_measure_of_correlation_I','Haralick:standard deviation intensity information_measure_of_correlation_I',...
    'Haralick:mean intensity information_measure_of_correlation_II','Haralick:standard deviation intensity information_measure_of_correlation_II',...
    'Haralick:mean intensity maximal_correlation_coefficient','Haralick:standard deviation intensity maximal_correlation_coefficient',...
    'entropy of periNuclei pixel intensity'};
save([folder_savepath 'periNuclei_feature_description.mat'],'description_periNuclei');
end