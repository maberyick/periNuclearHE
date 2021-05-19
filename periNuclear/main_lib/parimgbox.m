function [gray,R,G,B] = parimgbox(image, bbox)
    roi = image(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3),:);
    gray=rgb2gray(roi);
    R=roi(:,:,1);
    G=roi(:,:,2);
    B=roi(:,:,3);
end