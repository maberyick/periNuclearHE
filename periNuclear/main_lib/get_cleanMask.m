function [M] = get_cleanMask(I,M,maskFile)
    M = imbinarize(M);
    M = logical(M);
    M = bwareafilt(M,[30 300]);
    M = bwareaopen(M,10);
    M1 = getWatershedMask(I,true,6,10);
    M1 = logical(M1);
    M1 = bwareafilt(M1,[30 300]);
    M1 = bwareaopen(M1,10);
    M = M+M1;
    imwrite(M,maskFile);
end