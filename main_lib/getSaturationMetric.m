function sat=getSaturationMetric(I)
hsv = rgb2hsv(I);
sat=sum(sum(hsv(:,:,2)>.1));
end
