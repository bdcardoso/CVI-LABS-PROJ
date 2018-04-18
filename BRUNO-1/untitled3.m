
div = 0;

N = 20;

S = imread('Moedas3.jpg');
%S = im2bw(S);
%S2 = ~S;
%imshow(S2);

d= imbinarize(S)

denoiseImg1 = sum (d,3)/N;
figure, imagesc(denoiseImg1); colormap gray




%props = regionprops(d, 'Area');

%figure, imshow(props)
%S3 = bwlabel(S2);
%imagesc(S3);



%K1 = imfilter(S, fspecial('average', 3));

%K2 = medfilt2(S, [3 3]);
%figure, imshow(K1)
%figure, imshow(K2)

