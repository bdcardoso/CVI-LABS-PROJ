close all; 
clear all; 
clc;

I1 = imread('Moedas4.jpg');
figure, imshow(I1);

I=rgb2gray(I1);


[h,w]=size(I);
%figure;imshow(I);

c = edge(I, 'canny',0.4);  

I1 = cat(1,I1,c);
figure, imshow(I1);
%I1 = imsharpen(I1,'Radius', 2, 'Amount',1);
%figure, imshow(I1);

lab_he = rgb2lab(I1);
ab = lab_he(:,:,2:3);
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

nColors = 3;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
% - -

pixel_labels = reshape(cluster_idx,nrows,ncols);

figure, imshow(pixel_labels,[]) , title('image labeled by cluster index')


minD = min(pixel_labels(:));
maxD = max(pixel_labels(:));

mapped_image = (double(pixel_labels) - minD) ./ (maxD - minD);

ncmap = size(colormap, 1);

mapped_image = mapped_image .* ncmap;

if ncmap == 2
   mapped_image = mapped_image >= 0.5;  %logical
elseif ncmap <= 256
   mapped_image = uint8(mapped_image);
else
   mapped_image = uint16(mapped_image);
end

imwrite( mapped_image, 'YourData.tif' )



%imwrite(pixel_labels, 'YourData.jpg' );

I2 = imread('YourData.tif');
%I2 = rgb2gray(I2);
T1 = graythresh(I2);
c = edge(I2, 'canny',0.01);  

%d2 = imfill(c, 'holes');  

SE = strel('square',3);
SE = [0 1 0
      1 1 1
      0 1 0];
BW = imdilate(c,SE);

%BW2 =imerode(BW,SE);

d2 = imfill(BW, 'holes');  

[lb num] = bwlabel(c); %estava BW3

stats = regionprops(lb);
areas = [stats.Area];
[dummy indM] = max(areas);
imgBr = (lb == indM);
n_ojb=0
areas = [];
figure; imshow(lb); hold on;



%figure; imshow(lb);title(''); hold on;
[val ind] = max(areas)

%TESTING TESTING TESTING - LIGAR OBJECTOS 

ma = cellstr(num2str(classlabel(:)));    %num2str can process vector

boundaries = bwboundaries(d2);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 3);
    plot(stats(k).Centroid(1),stats(k).Centroid(2), '.r','markersize',25);
    
    plot(x(i),y(i),'Marker',ma{k},'MarkerSize',fontsize,'Color',colors(7));
    
end
hold off;

fprintf('BOUNDARIES NUMBER IS %d \n',numberOfBoundaries);


% Define object boundaries
numberOfBoundaries = size(boundaries, 1)
% message = sprintf('Found %d boundaries', numberOfBoundaries);
% uiwait(helpdlg(message));
% Find minimum distance between each pair of boundaries
for b1 = 1 : numberOfBoundaries
	for b2 = 1 : numberOfBoundaries
		if b1 >= b2
			% Can't find distance between the region and itself
			continue;
		end
		boundary1 = boundaries{b1};
		boundary2 = boundaries{b2};
		boundary1x = boundary1(:, 2);
		boundary1y = boundary1(:, 1);
		x1=1;
		y1=1;
		x2=1;
		y2=1;
		overallMinDistance = inf; % Initialize.
		% For every point in boundary 2, find the distance to every point in boundary 1.
		for k = 1 : size(boundary2, 1)
			% Pick the next point on boundary 2.
			boundary2x = boundary2(k, 2);
			boundary2y = boundary2(k, 1);
			% For this point, compute distances from it to all points in boundary 1.
			allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
			% Find closest point, min distance.
			[minDistance(k), indexOfMin] = min(allDistances);
			if minDistance(k) < overallMinDistance
				x1 = boundary1x(indexOfMin);
				y1 = boundary1y(indexOfMin);
				x2 = boundary2x;
				y2 = boundary2y;
				overallMinDistance = minDistance(k);
			end
		end
		% Find the overall min distance
		minDistance = min(minDistance);
		% Report to command window.
        
		fprintf('The minimum distance from region %d to region %d is %.3f pixels\n', b1, b2, minDistance);

		% Draw a line between point 1 and 2
        color = char('r','g','b','y','m','c','w','k');


        
        
		line([x1, x2], [y1, y2], 'Color', 'y', 'LineWidth', 1);
	end
end

%annotation('Found %d regions\n', b1,[2 3 300 200])
fprintf('Found %d regions\n', b1);

 