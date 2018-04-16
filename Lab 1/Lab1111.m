img = imread ('lena512.bmp');
figure, imshow(img);

image = imnoise(img, 'salt & pepper', 0.05);
figure, imshow(image);

K1 = imfilter(image, fspecial('average', 3));

K2 = medfilt2(image, [3 3]);
figure, imshow(K1)
figure, imshow(K2)