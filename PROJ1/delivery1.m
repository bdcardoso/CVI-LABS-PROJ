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
I = imread('Moedas1.jpg');

figure, imshow(I),title("Original Image");

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

figure;imshow(I);title("Image Subtract with erode");
T1 = graythresh(I) ;
%fprintf('%d',T1);


%c = edge(I, 'canny', T1);  
%T1 = imgaussfilt(T1, 5);
c = edge(I, 'canny',T1); 

%EDITA-ME !!!!!

figure; imshow(c);title("Image with edge"); 


BW = imdilate(c,SE);

BW = imdilate(BW,SE);

BW = imdilate(BW,SE);

BW2 =imerode(BW,SE);


bw_a = padarray (BW, [1 1], 1, 'pre' );
bw_a_filled = imfill (bw_a, 'holes' );
bw_a_filled = bw_a_filled (2: end, 2: end);
 
bw_b = padarray (padarray (BW, [1 0], 1, 'pre' ), [0 1], 1, 'post' );
bw_b_filled = imfill (bw_b, 'holes' );
bw_b_filled = bw_b_filled (2: end, 1: end-1);
 
bw_c = padarray (BW, [1 1], 1, 'post' );
bw_c_filled = imfill (bw_c, 'holes' );
bw_c_filled = bw_c_filled (1: end-1,1: end-1);

bw_d = padarray (padarray (BW, [1 0], 1, 'post' ), [0 1], 1, 'pre' );
bw_d_filled = imfill (bw_d, 'holes' );
bw_d_filled = bw_d_filled (1: end-1,2: end);
BW3 = bw_a_filled |  bw_b_filled |  bw_c_filled |  bw_d_filled;


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
d3 = imfill(BW3, 'holes');

figure, imshow(d1),title("Image D1 fill"); 
figure, imshow(d2),title("Image D2 fill");
figure, imshow(d3),title("Image D2 fill"); 

Label=bwlabel(d2,4);

[lb num] = bwlabel(d3);
[lb1 num1] = bwlabel(d2); %estava BW3
[lb2 num2] = bwlabel(d1); %estava BW3

if num < num1
    num = num1;
    lb=lb1;
end

if num < num2
    num = num2;
    lb =lb2;
end

%fiz 3 dilates, faÁo 3 erodes

figure;imshow(lb);

stats = regionprops(lb,'Perimeter','Area','Centroid','BoundingBox');
Areas = [stats.Area];
Perimetros = [stats.Perimeter]
[dummy indM] = max(Areas);
imgBr = (lb == indM);
n_ojb=0
areas = [];
figure; imshow(lb);title("Centroid and boundaries"); hold on;
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

%Contar dinheiro

%Moedas aprox = Matrix n*2   n = moedas , 2 = [area perimetro]

Moedas=[10838 368
        14175 420
        17615 468
        15180 435
        19375 490
        23050 535
        21095 512];
    
Valor = [0.01
         0.02
         0.05
         0.10
         0.20
         0.50
         1];

Dinheiro = 0
for k=1: num
    areas = [Areas length(find(lb==k))]
    peri = [Perimetros length(find(lb==k))]
    if all(areas(k)>100)
        for i=1:size(Moedas,1)
            difArea = abs(Moedas(i) - areas(k))
            difPeri = abs(Moedas(i+7) - peri(k))
            if(difArea <=1000 && difPeri <= 10)
                Dinheiro = Dinheiro + Valor(i);
            end
        end
    end

end

fprintf('Dinheiro Total: %.3g Ä\n', Dinheiro);


% Detectar Sharpness e ordenar por Sharpness




