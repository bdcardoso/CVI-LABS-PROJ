close all, clear all;

img = imread('rabbitBW.jpg');
figure, hold on, imshow(img);
%pause
se = strel('disk',3);

for k=1:30
    k
    img = imerode(img, se);
    imshow(img); drawnow
    pause(.2)
end

