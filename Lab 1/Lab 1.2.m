img = imread ('eight.tif');
figure, imshow(img);

imageN = imnoise(img, 'salt & pepper', 0.09);
figure, imshow(imageN);
K =filter2 (fspecial ('average', 4), imageN)/255;
figure, imshow(K);
imageF2 = medfilt2(imageN);
figure, imshow(imageF2);
