function bounds = mask2bounds_clean(img,isobjectblack)
%% takes in an image mask and extracts boundaries and centroids from it
% addpath(genpath('I:\code\lcib_1.0'))
% 
% uiopen('C:\Users\George\Desktop\54987_Progressor1_seg.jpg',1)
% 
% uiopen('C:\Users\gcl25\Downloads\x54987_NonProgressor1_seg.tif',1)

% img = x54987_NonProgressor1_seg;
% img = x54987_Progressor1_seg;

% threshold - Only if the Image comes from the DL Nuclei Model
%pres = max(double(img(:)));
%img(img > pres/2) = pres;
%img(img < pres/2) = 0;

%rgb2gray(img);
%binary = rgb2gray(img);
%if nargin < 2
%    mask = logical(~img); % object is black
%else
%    mask = logical(img); % object is white
%end

%size(mask)
% If the objects are black and the background is white!
% mask = logical(~img); % object is black
% If the objects are white and the background is black!
mask = logical(img); % object is white
conncomp = bwconncomp(mask);
L = labelmatrix(conncomp);

% if matlabpool('size') == 0
%     matlabpool(4)
% end

for i = 1:conncomp.NumObjects
    
    % get connected components from mask
    boundvec = conncomp.PixelIdxList{i};

    % George's old method
%     % prep masks for each connected component
%     temp = false(size(mask)); % initialize
%     vec = temp(:); % vectorize
%     vec(boundvec) = 1; % set mask
%     bimg = reshape(vec,size(mask)); %new mask

% Patrick's new method
bimg = L == i;

    % calculate bounds
    %[bnds] = trace_bound(bimg);
%     [bnds] = trace_bound_c(bimg);
%     bounds(i).r_in = bnds(:,1)';
%     bounds(i).c_in = bnds(:,2)';
    
    bimg = imdilate(bimg,strel('square',3));

    %[bnds] = trace_bound(bimg);
    [bnds] = trace_bound_c(bimg);
    bounds(i).r = bnds(:,1)';
    bounds(i).c = bnds(:,2)';
    
    %plot([bounds(i).r; bounds(i).r(1)], [bounds(i).c; bounds(i).c(1)])
    
    % calculate centroids
    S = regionprops(bimg,'Centroid');
    bounds(i).centroid_r = S.Centroid(1);
    bounds(i).centroid_c = S.Centroid(2);
    %plot(bounds(i).centroid_r,bounds(i).centroid_c)
end

% for i = 1:numel(bounds)
%     plot([bounds(i).r bounds(i).r(1)], [bounds(i).c bounds(i).c(1)]); hold on
%     plot(bounds(i).centroid_c,bounds(i).centroid_r,'r.');
% end


