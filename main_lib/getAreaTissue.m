function BW_area=getAreaTissue(I)
level = graythresh(I);
BW = imbinarize(I,level);
BW = BW(:,:,1).*BW(:,:,2).*BW(:,:,3);
BW_area = sum(sum(~BW));
end
