function [roi] = parimgboxroi(image, bbox, bbx_rad)
    roi = image(bbox(2)-bbx_rad : bbox(2) + bbox(4)+bbx_rad, bbox(1)-bbx_rad : bbox(1) + bbox(3)+bbx_rad, :);
end