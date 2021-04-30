function feat_extraction_main(varargin)
% Paths
folder_matpatches = varargin{1};
folder_pyepistroma = varargin{2};
folder_matcellmask = varargin{3};
folder_savepath = varargin{4};
% Quality parameters for Nuclear and peri-nuclear features
tileSize = 2000; % Size of the patch
quality.redChannelLim=100; % ignore the current tile if the mean intensity of R channel is less than this threshold
quality.tissueLim=.1; % process the current tile only if it contains at least this percentage of tissue
quality.blurLimit=.2; % process the current tile only if its blurriness score is below this value
quality.saturationLim=tileSize*20;
quality.imgarea = 0.2*(tileSize*tileSize);
feat_extraction_wsi(folder_matpatches,folder_pyepistroma,folder_matcellmask,folder_savepath,quality)
end