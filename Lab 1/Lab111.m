
%img = imread('veiculoGray.jpg);
img = imread ('lena512.bmp');
figure, imshow(img);

image = imnoise(img, 'gaussian',0, 0.05);
figure, imshow(image);

K =filter2 ( image, [3 3]);
figure, imshow(K);

%image = imnoise(img, 'salt & pepper', 0.05);
%figure, imshow(image);

%K1 =imfilert2 (image, fspecial ('average', 9), imageN)/255;
%figure, imshow(K);
%imageF2 = medfilt2(imageN);
%figure, imshow(imageF2);
