function [feat_perinuclei] = get_basicfeat(M,HEtile)
nucleiMaskTile = logical(M);
maskperim = bwperim(nucleiMaskTile);
binfo = regionprops(maskperim,'ConvexHull', 'Centroid', 'ConvexArea', 'BoundingBox');
centers = cat(1,binfo.Centroid);
cellRadius = 15;
if size(centers,1) == 0
    tmpmask = zeros(size(nucleiMaskTile));
else
    tmpmask = createCirclesMask(size(nucleiMaskTile), centers, cellRadius*ones(size(centers,1),1));
end
maskPeriNuclei = tmpmask-nucleiMaskTile;
maskPeriNuclei = uint8(repmat(maskPeriNuclei,1,1,3));
[feat_perinuclei] = extract_perinuclearfeat(HEtile.*maskPeriNuclei, HEtile.*uint8(~nucleiMaskTile), cellRadius, [binfo.ConvexArea], binfo);
img = HEtile.*maskPeriNuclei;
img2 = HEtile.*uint8(~nucleiMaskTile);
Area =  [binfo.ConvexArea];

end