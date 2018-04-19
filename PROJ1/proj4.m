clear all; 
close all;

img = imread('Moedas4.jpg');
figure, imshow(img);

%I = rangefilt(img,[0 1 0;1 1 1;0 1 0]);
I = imread('Moedas4.jpg');
figure, imshow(I);

figure, imhist(I(:,:,1));
T1 = graythresh(I);

redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);

maxGrayLevelR = max(redChannel(:));
minGrayLevelR = min(redChannel(:));
% Convert percentage threshold into an actual number.
thrper = 0.75;
thresholdLevel = minGrayLevelR + thrper*(maxGrayLevelR  - minGrayLevelR);
binaryImageR = redChannel  > thresholdLevel;

maxGrayLevelG = max(greenChannel(:));
minGrayLevelG = min(greenChannel(:));
% Convert percentage threshold into an actual number.
thrper = 0.75;
thresholdLevel = minGrayLevelG + thrper*(maxGrayLevelG  - minGrayLevelG);
binaryImageG = greenChannel  > thresholdLevel;

maxGrayLevelB = max(blueChannel(:));
minGrayLevelB = min(blueChannel(:));

thresholdLevel = minGrayLevelB + thrper*(maxGrayLevelB  - minGrayLevelB);
binaryImageB = blueChannel  > thresholdLevel;

%red_BW = imfill(binaryImageR,'holes');
%green_BW = imfill(binaryImageG,'holes');
%blue_BW = imfill(binaryImageB,'holes');
%figure, imshow(red_BW); title('red')
%figure, imshow(green_BW); title('green')
%figure, imshow(blue_BW); title('blue')

%b_bw = im2bw(red_BW + green_BW + ~blue_BW,0.0001);
%figure, imshow(b_bw); title('some')

h = fspecial ('average');
imageF = imfilter(I,h);

It = im2bw(imageF,T1);
It2 = imfill(It,'holes');
figure, imshow(It2);

[lb num] = bwlabel(It2); %estava BW3
%figure,
%subplot(1,3,1); imshow(mat2gray(lb)); title('Etiquetas');
%subplot(1,3,1); imshow(label2rgb(lb));title('Labels')
stats = regionprops(lb);
areas = [stats.Area];
[dummy indM] = max(areas);
imgBr = (lb == indM);
%subplot(1,3,2);imshow(imgBr);title('Maior �rea');
%subplot(1,3,3);imshow(I);title('C�rebro');
%alternativa
%subplot(1,3,4);imshow(double(imgg).*(imgBr));title('C�rebro');
n_ojb=0
areas = [];
figure; imshow(lb); hold on;
for k=1: num
    areas = [areas length(find(lb==k))]
    if all(areas(k)>500)
        plot(stats(k).Centroid(1),stats(k).Centroid(2), '.r','markersize',25);
        drawnow;
        
    end

end
%figure; imshow(lb);title(''); hold on;
[val ind] = max(areas)

    
imgBr1 = (lb == ind);
%plot(1,1);imshow(imgBr1);title('Maior �rea');
%plot(1,2);imshow(I);title('C�rebro');


%plot(stats(ind).Centroid(1),stats(ind).Centroid(2), '.r','markersize',50);
%    drawnow;

%pause
