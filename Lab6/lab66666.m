clear all, close all;

I =imread('tire.tif');
J = histeq(I, 64);
imshow(I)
figure, imshow(J)