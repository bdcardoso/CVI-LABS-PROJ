clear all; 
close all;

I = imread('Moedas1.jpg');
figure, imshow(I);

figure, imhist(I(:,:,1));
T1 = graythresh(I);

% Extract the individual red, green, and blue color channels.
% (Note: rgbImage might be after you've taken the log of it.)
redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);

max1 = max(redChannel(:));
max2 = max(greenChannel(:));
max3 = max( blueChannel(:));

min1 = min(redChannel(:));
min2 = min(greenChannel(:));
min3 = min ( blueChannel(:));

maxGrayLevelR = min([max1 max2 max3]);
minGrayLevelR = max([min1 min2 min3]);
% Convert percentage threshold into an actual number.
thrper = 0.75;
thresholdLevel = minGrayLevelR + thrper*(maxGrayLevelR  - minGrayLevelR);
binaryImageR = redChannel  > thresholdLevel;

It = im2bw(I,T1);
It2 = imfill(binaryImageR,'holes');
figure, imshow(It2);

%It = im2bw(I,T1);
%figure, imshow(It);

[lb num] = bwlabel(It2); %estava BW3
%figure,
%subplot(1,3,1); imshow(mat2gray(lb)); title('Etiquetas');
%subplot(1,3,1); imshow(label2rgb(lb));title('Labels')
stats = regionprops(lb);
areas = [stats.Area];
[dummy indM] = max(areas);
imgBr = (lb == indM);
%subplot(1,3,2);imshow(imgBr);title('Maior área');
%subplot(1,3,3);imshow(I);title('Cérebro');
%alternativa
%subplot(1,3,4);imshow(double(imgg).*(imgBr));title('Cérebro');
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
%plot(1,1);imshow(imgBr1);title('Maior área');
%plot(1,2);imshow(I);title('Cérebro');


%plot(stats(ind).Centroid(1),stats(ind).Centroid(2), '.r','markersize',50);
%    drawnow;

%pause
