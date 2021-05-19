function [M_stro,M_epi, M_bund] = get_epistroma(I,M,ES,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile)
[q_r,q_m,~] = size(I);
ES = imresize(ES,[q_r q_m],'Antialiasing',false);
ES = uint8(255*mat2gray(ES)); 
ES = uint8(ES);
level = graythresh(ES);
ES = imbinarize(ES,level);
ES1 = imfill(ES,18,'holes');
ES1 = bwareaopen(ES1, 50);
%figure
%imshow(ES1)
%
%se = strel('disk',3);
%ES1 = ~imopen(~ES1,se);
%ES1 = imfill(ES1,18,'holes');
%ES1= ~imfill(~ES1,8,'holes');
% Epithelium
M_epi = M.*ES1;
M_epi = bwareaopen(M_epi, 50);
% Stroma
ESs = ~ES1;
M_stro = M.*ESs;
M_stro= bwareaopen(M_stro, 50);
% Boundary
ES_Bound = boundarymask(ES1);
ES_Bound1 = imdilate(ES_Bound, strel('disk',40));
M_bund = M.*ES_Bound1;
M_bund= bwareaopen(M_bund, 50);
% Save
% Tissue maks
imwrite(ES1,ESmaskFile);
imwrite(ESs,SSmaskFile);
imwrite(ES_Bound1,EBmaskFile);
% Cell masks
imwrite(M_epi,maskEFile);
imwrite(M_stro,maskSFile);
imwrite(M_bund,maskBFile);
end