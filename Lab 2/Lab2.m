imgg = imread('veiculoGray.jpg');

div = 0;

N = 40;

[L C] = size(imgg);
%MyImg = zeros(L, C, N);
MyImg(1:L,1:C,1) = imnoise(imgg, 'salt & pepper', 0.7);

for i=1:N
    %noise = round(randn(size(imgg))*div);
    %imshow(mat2gray(noise));
    %image(:,:,i) = max(min(imgg+uint8(noise), 255),0);
    %image(:,:,i) = imnoise(imgg, m'gaussian', 0, 0.12);
    
    
    image(:,:,i)=imnoise(imgg,'salt & pepper', 0.22);
    figure(1); 
    imshow(image(:,:,i));
end

denoiseImg1 = sum (image,3)/N;
figure, imagesc(denoiseImg1); colormap gray

denoiseImg2 = median (image,3);
figure, imagesc(denoiseImg2); colormap gray

    