function [nucleiCentroids,nucleiCentroids_stro,nucleiCentroids_epi, nucleiCentroids_bund, feat_perinuclei, feat_epi_perinuclei, feat_stro_perinuclei, feat_bund_perinuclei] = get_periNuclear(I,ES,M,maskFile,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile)
%GETPERINUCLEARFEATURES
% Clean image by watershed
% Check if the mask is already processed
if exist(maskFile,'file')~=2
    [M] = get_cleanMask(I,M,maskFile);
    M = imbinarize(M);
else
    % load mask
    M = imread(maskFile);
    M = imbinarize(M);
end
%%Adapt the ES mask to the nuclei
if [exist(maskEFile,'file') && exist(maskSFile,'file') && exist(maskBFile,'file')] == 0
    [M_stro,M_epi,M_bund] = get_epistroma(I,M,ES,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile);
else
    % load mask
    M_epi = imread(maskEFile);
    M_stro = imread(maskSFile);
    M_bund = imread(maskBFile);
end
% Get cell centroids - All cells
grayImg=rgb2gray(I);
regionProperties = regionprops(M,grayImg,'Centroid');
nucleiCentroids = cat(1, regionProperties.Centroid);
% Extract based on all tissue cells
 [feat_perinuclei] = get_basicfeat(M,I);
% Get cell centroids - All cells
regionProperties_epi = regionprops(M_epi,grayImg,'Centroid');
nucleiCentroids_epi = cat(1, regionProperties_epi.Centroid);
regionProperties_stro = regionprops(M_stro,grayImg,'Centroid');
nucleiCentroids_stro = cat(1, regionProperties_stro.Centroid);
regionProperties_bund = regionprops(M_bund,grayImg,'Centroid');
nucleiCentroids_bund = cat(1, regionProperties_bund.Centroid);
%% Extract features based on Epithelium
 [feat_epi_perinuclei] = get_basicfeat(M_epi,I);
%% Extract features based on Stroma
 [feat_stro_perinuclei] = get_basicfeat(M_stro,I);
%% Extract features based on epiStroma boundary
 [feat_bund_perinuclei] = get_basicfeat(M_bund,I);
 % Get the description_stro_perinuclei
end