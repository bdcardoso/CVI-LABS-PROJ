close all; clear all;
myPath = 'z:\Desktop\CVI\Lab5\Pedestrian\';

N=100;

img2 =imread(strcat(myPath,'ped7c1350.tif'));
figure, imshow(img2);
[nlin ncol dummy]=size(img2);
npixels = nlin*ncol;
hr=imhist(img2(:,:,1),N)/npixels;
hg=imhist(img2(:,:,2),N)/npixels;
hb=imhist(img2(:,:,3),N)/npixels;
H1 = [hr' hg' hb'];
figure, bar(H1);