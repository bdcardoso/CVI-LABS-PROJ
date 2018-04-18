Irgb = imread('Moedas1.jpg');
Igray = rgb2gray(Irgb);
Ibw = im2bw(Igray, 0.9);
se = strel('square',2);
Imorph = imerode(Ibw,se);
Iarea = bwareaopen(Imorph,130);
Idiff = Imorph - Iarea;
Ifinal = bwareaopen(Idiff,70);
stat = regionprops(Ifinal,'BoundingBox');
str = sprintf('number of detected nails : %d', length(stat));
imshow(Irgb); title(str); hold on;
for cnt = 1 : length(stat)
    
    bb = stat(cnt).BoundingBox;
    rectangle('position',bb,'edgecolor',rand(1,3));    
    
end