clc;
clear all;
[f, p] = uigetfile('*.jpg'); %retrive files
I = imread([p f]);
[pathstr,name,ext] = fileparts([p f]);

img=imresize(I,[512, 512]);

      img = rgb2gray(img);

if strcmp(name,'8') || strcmp(name,'9') || strcmp(name,'10') || strcmp(name,'11') || strcmp(name,'12') || strcmp(name,'13.jpg') || strcmp(name,'14') || strcmp(name,'15') || strcmp(name,'16')
    G = fspecial('gaussian',[6 6],2);
    img= imfilter(img,G,'same');
else
    G = fspecial('gaussian',[3 3],2);
    img= imfilter(img,G,'same');
end
