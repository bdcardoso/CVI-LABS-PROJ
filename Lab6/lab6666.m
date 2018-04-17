close all; clear all;
path = 'z:\Desktop\CVI\MATERIAL\Coffe\';

frameIdComp = 4;
str = ['%s%.' num2str(frameIdComp) 'd.%s'];

nFrame = 1048;
step = 1;

str1 = sprintf(str, path, 1, 'jpg');
I = imread(str1);
[L, C, Z] = size(I);

vid4D = zeros([240 320 3 nFrame/step]);
figure; hold on

for k=1 : step : -nFrame
    k
    str1 = sprintf(str,path,k,'jpg');
    img = imread (str1);
    vid4D(:,:,:,k) = img;
    imshow(img); drawnow
    %pause(.2)
    %disp('');
    
end
bkg = median(vid4D,4);
figure; imagesc(uint8(bkg));
    

Bkg = zeros(size(I));
alfa = 0.02;

for k = 1 : step : nFrame
    k
    str1 = sprintf(str,path,k,'jpg');
    img = imread (str1);
    
    Y = img;
    Bkg = alfa * double(Y) + (1-alfa) * double(Bkg);
    imshow(uint8(Bkg)); drawnow
end

    