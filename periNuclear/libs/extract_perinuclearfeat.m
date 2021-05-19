function [feats_periNuclei] =  extract_perinuclearfeat(img, img2, cellRadius, Area, binfo)
%cytoplasm area (periNuclei)
%the ratio of cytoplasm area (periNuclei) to nuclei area
%haralick features in periNuclei region
%entropy of the perinuclei region
%Notice: the haralick input is the masked rgb images, so it's necessary to
%generate the mask of peri-nuclei region. Needs to treat the area and
%texture feature separately, since the texture will calculate based on the
%co-occurence matrix of entire matrix. The area should calculated based per
%object
periAreas = pi*cellRadius^2-Area;
areaRatioCyto2Nuclei = periAreas./Area;
%calculating haralick feature
info.dist = 1;
info.win = 1;
info.grays = 256;
%n = 0;
if ndims(img) < 3
    gimg = img;
elseif ndims(img) == 3
    gimg = rgb2gray(img); % assume img is rgb
else
    fprintf('Unidentified image format')
end
%mask(:,:,1) = (gimg ~= max(max(gimg)));
grays = info.grays;
%win = info.win;
%dist = info.dist;
%himg = uint16(rescale_range(gimg,0,grays-1));
%% Adapted code for the new haralick
%centroids = cat(1, binfo.Centroid);
harFeat=[];
bbx_rad = 5;
parfor i=1:length(binfo)
    nucleus=binfo(i);
    bbox = nucleus.BoundingBox;
    bbox = [round(bbox(1)) round(bbox(2)) (bbox(3) - 1) (bbox(4) - 1)];
    try
        %roi = img2(bbox(2)-bbx_rad : bbox(2) + bbox(4)+bbx_rad, bbox(1)-bbx_rad : bbox(1) + bbox(3)+bbx_rad, :);
        roi = parimgboxroi(img2, bbox, bbx_rad)
    catch merror
        if ~isempty(merror)
            %roi = img2(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3), :);
            roi = parimgboxroi(img2, bbox, bbx_rad)
        end
    end
    %merror = [];
    %imshow(roi)
    gray=rgb2gray(roi);
    % Haralick features
    glcm = graycomatrix(gray);    
    harFeat=[harFeat;haralickTextureFeatures(glcm,1:14)'];
end
%%
[harfeta_size,~] = size(harFeat);
if harfeta_size == 1
    % Instead of using std of 1 single value, which is 0, the whole values
    % are being passed
    Texturefeats = [harFeat harFeat];
else
    Texturefeats = [mean(harFeat) std(harFeat)];
end
J =  entropy(img);
if isempty(binfo) == 1
    feats_periNuclei  = zeros(1,37);
else
    feats_periNuclei = [mean(periAreas), median(periAreas), var(periAreas), std(periAreas),...
    mean(areaRatioCyto2Nuclei), median(areaRatioCyto2Nuclei), var(areaRatioCyto2Nuclei), std(areaRatioCyto2Nuclei),...
    Texturefeats, J];
end
end