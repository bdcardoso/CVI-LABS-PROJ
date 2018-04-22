close all; 
clear all; 
clc;



%
%
%SEC√ß√ÉO DOS TESTES1
%
%

%get the input image
% I=imread('Moedas4.jpg');
% 
% imshow(I),title('Image:1');
% 
% %change the color space
% cform = makecform('srgb2lab');
% K=applycform(I,cform);
% 
% figure;imshow(K),title('Image:2');
% 
% %equalise brightness to get skin area
% %K=J(:,:,2);% 2nd page of 3-d vector j
% 
% figure;imshow(K),title('Image:3');
% 
% L=graythresh(K(:,:,2));% find appropriate gray thresh value
% BW1=im2bw(K(:,:,2),L);% convert to binary image based on threshold
% 
% figure;imshow(BW1),title('Image:4');
% 
% bw2=imfill(BW1,'holes');% fill patches with holes
% figure;imshow(bw2)
% bw3 = bwareaopen(bw2,1890); %opens area greater than 1890
% cc=bwconncomp(bw3)% connected comp for finding the density of people in image
% density=cc.NumObjects / (size(bw3,1) * size(bw3,2))
% figure;imshow(bw3)
% labeledImage = bwlabel(bw3, 8);%same as connected components
% figure;imshow(labeledImage)
% blobMeasurements = regionprops(labeledImage,'all');%measure all properties of the image
% numberOfPeople = size(blobMeasurements, 1)% count the number of people
% % draw bounding boxes
% imagesc(I);
% 
% hold on;
% title('Original with bounding boxes');
% for k = 1 : numberOfPeople % Loop through all blobs.
% % Find the mean of each blob.
% % directly into regionprops.
% thisBlobsBox = blobMeasurements(k).BoundingBox; % Get list of pixels in current blob.
% x1 = thisBlobsBox(1);%1st side
% y1 = thisBlobsBox(2);%2nd side
% x2 = x1 + thisBlobsBox(3);%3rd side
% y2 = y1 + thisBlobsBox(4);%4th side
% x = [x1 x2 x2 x1 x1];
% y = [y1 y1 y2 y2 y1];
% %subplot(3,4,2);
% plot(x, y, 'LineWidth', 2);
% end


%
%
%SEC√ß√ÉO DOS TESTES2
%
%
% 
% A=imread('Moedas4.jpg');
% B=im2bw(A);
% C=imfill(B,'holes');
% 
% label=bwlabel(C);
% max(max(label))
% im1= (label==1);
% 
% for j=1:max(max(label))
%     [row, col] = find(label==j);
%     len=max(row)-min(row)+2;
%     breadth=max(col)-min(col)+2;
%     target=uint8(zeros([len breadth] ));
%     sy=min(col)-1;
%     sx=min(row)-1;
% 
%     for i=1:size(row,1)
%         x=row(i,1)-sx;
%         y=col(i,1)-sy;
%         target(x,y)=A(row(i,1),col(i,1));
%     end
%     
%     mytitle=strcat('Object Number:',num2str(j));
%     figure,imshow(target);title(mytitle);
% end


%
%
%SEC√ß√ÉO DOS TESTES3
%
%


%https://www.mathworks.com/help/images/examples/marker-controlled-watershed-segmentation.html


% I=imread('Moedas4.jpg');
% 
% T99=rgb2gray(I);
% T2=edge(T99,'Canny');
% 
% hy = fspecial('prewitt');
% hx = hy';
% Iy = imfilter(double(T2), hy, 'replicate');
% Ix = imfilter(double(T2), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% figure
% imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
% 
% 
% 
% se = strel('disk', 20);
% Io = imopen(T2, se);
% figure
% imshow(Io), title('Opening (Io)')
% 
% Ie = imerode(T2, se);
% Iobr = imreconstruct(Ie, T2);
% figure
% imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
% % 
% % Ioc = imclose(Io, se);
% % figure
% % imshow(Ioc), title('Opening-closing (Ioc)')
% % 
% % 
% % se = strel('disk', 20);
% % Io = imopen(T2, se);
% % figure
% % imshow(Io), title('Opening (Io)')
% % 
% % Iobrd = imdilate(Iobr, se);
% % Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
% % Iobrcbr = imcomplement(Iobrcbr);
% % figure
% % imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
% 
% 
% 
% 
% 
% %kernel = -1*ones(3);
% %kernel(2,2) = 17;
% %T = imfilter(T2, kernel);
% 
% %T2=im2bw(T,graythresh(T));
% %T2=~T2;
% 
% %B = bwboundaries(T2);
% B=bwboundaries(Iobr);
% imshow(T2)
% text(10,10,strcat('\color{green}Objects Found:',num2str(length(B))))
% hold on
% 
% for k = 1:length(B)
% boundary = B{k};
% plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 0.2)
% end
% 
% 
% 
% %
% %
% % % NOVO TESTE
% %
% %

% I = imread('Moedas3.jpg');
% imshow(I)
% axis off
% title('Original Image')
% 
% thresh = multithresh(I,2);
% seg_I = imquantize(I,thresh);
% 
% 
% 
% 
% 
% 
% RGB = label2rgb(seg_I); 	 
% figure;
% imshow(RGB)
% axis off
% title('RGB Segmented Image')

% %
% %
% % FIM SECCAO DE TESTES
% %
% %


%PARTE FUNCIONAL
I = imread('Moedas4.jpg');

background = imopen(I,strel('disk',150));

I = imsubtract(I,background);

figure, imshow(I);

I=rgb2gray(I);

SE = strel('square',3);
SE = [0 1 0
      1 1 1
      0 1 0];
  

I =imerode(I,SE);

[h,w]=size(I);

figure;imshow(I);
T1 = graythresh(I) ;
%fprintf('%d',T1);


%c = edge(I, 'canny', T1);  
%T1 = imgaussfilt(T1, 5);
c = edge(I, 'canny',T1); 

%EDITA-ME !!!!!

figure; imshow(c); 


BW = imdilate(c,SE);

BW2 =imerode(BW,SE);

%figure;imshow(BW2);

%se = strel('disk',0);      
%I2 = imdilate(c,se);       
%imshow(I2); 
% bw_a = padarray(I,[1 1],1,'pre');
% bw_a_filled = imfill(bw_a,'holes');
% bw_a_filled = bw_a_filled(2:end,2:end);
% imshow(bw_a_filled)
% bw_b = padarray(padarray(I,[1 0],1,'pre'),[0 1],1,'post');
% bw_b_filled = imfill(bw_b,'holes');
% bw_b_filled = bw_b_filled(2:end,1:end-1);
% imshow(bw_b_filled);
% bw_c = padarray(I,[1 1],1,'post');
% bw_c_filled = imfill(bw_c,'holes');
% bw_c_filled = bw_c_filled(1:end-1,1:end-1);
% imshow(bw_c_filled)
% bw_d = padarray(padarray(I,[1 0],1,'post'),[0 1],1,'pre');
% bw_d_filled = imfill(bw_d,'holes');
% bw_d_filled = bw_d_filled(1:end-1,2:end);
% imshow(bw_d_filled)
% bw_filled = bw_a_filled | bw_b_filled | bw_c_filled | bw_d_filled;
% BW2 = padarray(BW2,[1 1],1,'pre');
% BW2 = padarray(padarray(BW2,[1 0],1,'pre'),[0 1],1,'post');
% BW2 = padarray(BW2,[1 1],1,'post');
% BW2 = padarray(padarray(BW2,[1 0],1,'post'),[0 1],1,'pre');
% figure;imshow(BW2);title("TESTE");

d1 = imfill(c, 'holes');
d2 = imfill(BW2, 'holes');  
figure, imshow(d1);        

Label=bwlabel(d2,4);


[lb num] = bwlabel(d2); %estava BW3
[lb1 num1] = bwlabel(d1); %estava BW3

if num < num1
    num = num1;
    lb = lb1;
end

stats = regionprops(lb,'Perimeter','Area','Centroid','BoundingBox');
Areas = [stats.Area];
Perimetros = [stats.Perimeter]
[dummy indM] = max(Areas);
imgBr = (lb == indM);
n_ojb=0
areas = [];
figure; imshow(lb); hold on;
i = 0
for k=1: num
    areas = [Areas length(find(lb==k))]
    if all(areas(k)>100)
        i = i+1;
        text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(i),'HorizontalAlignment','center','Color','r');
    end

end

%figure; imshow(lb);title(''); hold on;
[val ind] = max(Areas)

%TESTING TESTING TESTING - LIGAR OBJECTOS 

boundaries = bwboundaries(d2);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 1);
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
        color = char('r','g','b','y','m','c','w');


        
        
		line([x1, x2], [y1, y2], 'Color', 'y', 'LineWidth', 1);
    end
end

   

%annotation('Found %d regions\n', b1,[2 3 300 200])
fprintf('Found %d regions\n', b1);

% Imprimir Areas

%N„o imprime o "6" nesta imagem pois ele n„o est· a fazer o calculo da area
%bem para esse elemento

figure; imshow(lb);title('Area'); hold on;

for k=1: num
    areas = [areas length(find(lb==k))]
    if all(areas(k)>200)
        text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(areas(k)),'HorizontalAlignment','center','Color','r');
    end

end

hold off;

% Imprimir ¡rea por ordem

figure; imshow(lb);title('¡rea Ordem'); hold on;

areasord = sort(Areas,'descend');

for k=1: num
    areas = [areas length(find(lb==k))]
    if all(areas(k)>200)
        i = find(areasord == areas(k));
        text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(i),'HorizontalAlignment','center','Color','r');
    end

end

hold off;

% Imprimir Perimetros

figure; imshow(lb);title('Perimeter'); hold on;

perisord = sort(Perimetros,'descend');

for k=1: num
    areas = [Perimetros length(find(lb==k))]
    if all(areas(k)>100)
        text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(areas(k)),'HorizontalAlignment','center','Color','r');
    end

end

hold off;

% Imprimir Perimetros por ordem

figure; imshow(lb);title('Perimeter por Ordem'); hold on;

perisord = sort(Perimetros,'descend');

for k=1: num
    areas = [Perimetros length(find(lb==k))]
    if all(areas(k)>100)
        i = find(perisord == areas(k));
        text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(i),'HorizontalAlignment','center','Color','r');
    end

end

hold off;

% Imprimir Areas


%DETECCAO DE MOEDAS DAQUI PARA BAIXO

% a1=(Label==1);
% a2=(Label==2);
% a3=(Label==3);
% a4=(Label==4);
% a5=(Label==5);
% a6=(Label==6);
% a7=(Label==7);
% a8=(Label==8);
% a9=(Label==9);
% a10=(Label==10);
% 
% 
% D1 = bwdist(~a1);           
% %figure, imshow(D1,[]),        
% [xc1 yc1 r1]=distance(D1);
% f1=whichcoin(r1)
% disp(r1)
% 
% 
% D2 = bwdist(~a2);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D2,[]),      %  
% [xc2 yc2 r2]=distance(D2);
% f2=whichcoin(r2)
% disp(r2)
% 
% D3 = bwdist(~a3);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D3,[]),      %  
% [xc3 yc3 r3]=distance(D3);
% f3=whichcoin(r3)
% disp(r3)
% 
% D4 = bwdist(~a4);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D4,[]),      %  
% [xc4 yc4 r4]=distance(D4);
% f4=whichcoin(r4)
% disp(r4)
% 
% D5 = bwdist(~a5);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D5,[]),      %  
% [xc5 yc5 r5]=distance(D5);
% f5=whichcoin(r5)
% disp(r5)
% 
% D6 = bwdist(~a6);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D6,[]),      %  
% [xc6 yc6 r6]=distance(D6);
% f6=whichcoin(r6)
% disp(r6)
% 
% D7 = bwdist(~a7);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D7,[]),      %  
% [xc7 yc7 r7]=distance(D7);
% f7=whichcoin(r7)
% disp(r7)
% 
% D8 = bwdist(~a8);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D8,[]),      %  
% [xc8 yc8 r8]=distance(D8);
% f8=whichcoin(r8)
% disp(r8)
% 
% D9 = bwdist(~a9);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D9,[]),      %  
% [xc9 yc9 r9]=distance(D9);
% f9=whichcoin(r9)
% disp(r9)
% 
% D10 = bwdist(~a10);           % computing minimal euclidean distance to non-white pixel 
% %figure, imshow(D10,[]),      %  
% [xc10 yc10 r10]=distance(D10);
% f10=whichcoin(r10)
% disp(r10)
% 
% 
% 
