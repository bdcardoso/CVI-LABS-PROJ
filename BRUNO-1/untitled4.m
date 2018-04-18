close all; 
clear all; 
clc;

  

I =imread('Moedas3.jpg');

I=rgb2gray(I);


[h,w]=size(I);
figure;imshow(I);

c = edge(I, 'canny',0.4);  
figure; imshow(c);         

%se = strel('disk',0);      
%I2 = imdilate(c,se);       
%imshow(I2);  


d2 = imfill(c, 'holes');  
figure, imshow(d2);        

Label=bwlabel(d2,4);


[lb num] = bwlabel(d2); %estava BW3

stats = regionprops(lb);
areas = [stats.Area];
[dummy indM] = max(areas);
imgBr = (lb == indM);
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

%TESTING TESTING TESTING - LIGAR OBJECTOS 

boundaries = bwboundaries(d2);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 3);
end
hold off;

% Define object boundaries
numberOfBoundaries = size(boundaries, 1)
% message = sprintf('Found %d boundaries', numberOfBoundaries);
% uiwait(helpdlg(message));
% Find minimum distance between each pair of boundaries
for b1 = 1 : numberOfBoundaries
	for b2 = 2 : numberOfBoundaries
		if b1 == b2
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
		line([x1, x2], [y1, y2], 'Color', 'y', 'LineWidth', 3);
	end
end


%DETECCAO DE MOEDAS DAQUI PARA BAIXO

a1=(Label==1);
a2=(Label==2);
a3=(Label==3);
a4=(Label==4);
a5=(Label==5);
a6=(Label==6);
a7=(Label==7);
a8=(Label==8);
a9=(Label==9);
a10=(Label==10);


D1 = bwdist(~a1);           
figure, imshow(D1,[]),        
[xc1 yc1 r1]=distance(D1);
f1=whichcoin(r1)
disp(r1)


D2 = bwdist(~a2);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D2,[]),      %  
[xc2 yc2 r2]=distance(D2);
f2=whichcoin(r2)
disp(r2)

D3 = bwdist(~a3);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D3,[]),      %  
[xc3 yc3 r3]=distance(D3);
f3=whichcoin(r3)
disp(r3)

D4 = bwdist(~a4);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D4,[]),      %  
[xc4 yc4 r4]=distance(D4);
f4=whichcoin(r4)
disp(r4)

D5 = bwdist(~a5);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D5,[]),      %  
[xc5 yc5 r5]=distance(D5);
f5=whichcoin(r5)
disp(r5)

D6 = bwdist(~a6);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D6,[]),      %  
[xc6 yc6 r6]=distance(D6);
f6=whichcoin(r6)
disp(r6)

D7 = bwdist(~a7);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D7,[]),      %  
[xc7 yc7 r7]=distance(D7);
f7=whichcoin(r7)
disp(r7)

D8 = bwdist(~a8);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D8,[]),      %  
[xc8 yc8 r8]=distance(D8);
f8=whichcoin(r8)
disp(r8)

D9 = bwdist(~a9);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D9,[]),      %  
[xc9 yc9 r9]=distance(D9);
f9=whichcoin(r9)
disp(r9)

D10 = bwdist(~a10);           % computing minimal euclidean distance to non-white pixel 
figure, imshow(D10,[]),      %  
[xc10 yc10 r10]=distance(D10);
f10=whichcoin(r10)
disp(r10)



