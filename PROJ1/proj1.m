close all; 
clear all; 
clc;

%PARTE FUNCIONAL
file_name='Moedas3.jpg';

I = imread(file_name);

figure, imshow(I),title('Original Image');

background = imopen(I,strel('disk',150));

I = imsubtract(I,background);

%figure, imshow(I);

I=rgb2gray(I);

SE = strel('square',3);
SE = [0 1 0
      1 1 1
      0 1 0];
  

I =imerode(I,SE);


[h,w]=size(I);

%figure;imshow(I);title('Image Subtract with erode');
T1 = graythresh(I) ;
%fprintf('%d',T1);

%c = edge(I, 'canny', T1);  
%T1 = imgaussfilt(T1, 5);
c = edge(I, 'canny',T1); 

%EDITA-ME !!!!!

%figure; imshow(c);title('Image with edge'); 


BW = imdilate(c,SE);

BW = imdilate(BW,SE);

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


d1 = imfill(c, 'holes');
d2 = imfill(BW2, 'holes'); 
d3 = imfill(BW3, 'holes');

%figure, imshow(d1),title('Image D1 fill'); 
%figure, imshow(d2),title('Image D2 fill');
%figure, imshow(d3),title('Image D3 fill'); 

Label=bwlabel(d2,4);

[lb num] = bwlabel(d3);
[lb1 num1] = bwlabel(d2); %estava BW3
[lb2 num2] = bwlabel(d1); %estava BW3
flag = 1;
al = num+2;
if num > num1 || num > num2
    if al > num1
        num = num1;
        lb=lb1;
        flag = 0;
    end

    if al > num2
        num = num2;
        lb =lb2;
        flag = 0;
    end
end

%fiz 3 dilates, faço 3 erodes
if flag == 1
    lb =imerode(lb,SE);
    lb =imerode(lb,SE);
    lb =imerode(lb,SE);
end


%figure;imshow(lb);

stats = regionprops(lb,'Perimeter','Area','Centroid','BoundingBox', 'ConvexHull');
Areas = [stats.Area];
Perimetros = [stats.Perimeter];
[dummy indM] = max(Areas);
imgBr = (lb == indM);
n_ojb=0;
areas = [];
figure; imshow(lb);title('Centroid and boundaries'); hold on;
i = 0;
for k=1: num
    areas = [Areas length(find(lb==k))];
    if all(areas(k)>100)
        i = i+1;
        plot(stats(k).Centroid(1),stats(k).Centroid(2), '.r','markersize',25);
    end

end

%figure; imshow(lb);title(''); hold on;
[val ind] = max(Areas);

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
numberOfBoundaries = size(boundaries, 1);
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
%fprintf('Found %d regions\n', b1);

% Imprimir Areas

%Não imprime o '6' nesta imagem pois ele não está a fazer o calculo da area
%bem para esse elemento

% figure; imshow(lb);title('Area'); hold on;
% 
% for k=1: num
%     areas = [areas length(find(lb==k))];
%     if all(areas(k)>200)
%         
%         
%         text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(areas(k)),'HorizontalAlignment','center','Color','r');
%     end
% 
% end
% 
% hold off;

% Imprimir Área por ordem

figure; imshow(lb);title('Área Ordem'); hold on;

areasord = sort(Areas,'descend');

for k=1: num
    areas = [areas length(find(lb==k))];
    if all(areas(k)>200)
        i = find(areasord == areas(k));
        texto =  strcat(strcat(num2str(i),' - '), num2str(areas(k)));
        text(stats(k).Centroid(1),stats(k).Centroid(2),texto,'HorizontalAlignment','center','Color','r');
    end

end

hold off;

% Imprimir Perimetros
% 
% figure; imshow(lb);title('Perimeter'); hold on;
% 
% perisord = sort(Perimetros,'descend');
% 
% for k=1: num
%     areas = [Perimetros length(find(lb==k))]
%     if all(areas(k)>100)
%         text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(areas(k)),'HorizontalAlignment','center','Color','r');
%     end
% 
% end
% 
% hold off;

% Imprimir Perimetros por ordem

figure; imshow(lb);title('Perimeter por Ordem'); hold on;

perisord = sort(Perimetros,'descend');

for k=1: num
    areas = [Perimetros length(find(lb==k))];
    if all(areas(k)>100)
        i = find(perisord == areas(k));
        texto =  strcat(strcat(num2str(i),' - '), num2str(areas(k)));
        text(stats(k).Centroid(1),stats(k).Centroid(2),texto,'HorizontalAlignment','center','Color','r');
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

Dinheiro = 0;
for k=1: num
    areas = [Areas length(find(lb==k))];
    peri = [Perimetros length(find(lb==k))];
    if all(areas(k)>100)
        for i=1:size(Moedas,1)
            difArea = abs(Moedas(i) - areas(k));
            difPeri = abs(Moedas(i+7) - peri(k));
            if(difArea <=1000 && difPeri <= 10)
                %fprintf('Moeda: %.3g €\n', Valor(i));
                Dinheiro = Dinheiro + Valor(i);
            end
        end
    end

end

fprintf('Dinheiro Total: %.3g €\n', Dinheiro);


% % Detectar Sharpness e ordenar por Sharpness
% G1 = imread('Moedas4.jpg');
% G = double(rgb2gray(G1));
% [Gx, Gy]=gradient(G);
% 
% S=sqrt(Gx.*Gx+Gy.*Gy);
% sharpness=sum(sum(S))./(numel(Gx));


figure; imshow(lb);title('Sharpeness por Ordem'); hold on;

image = lb;
sharp_list = [];
imagem = imread(file_name);


for k=1: num
    i = (lb == k);
    image = image | i;
    ImCrp = imagem.*repmat(uint8(i),[1,1,3]);
    G = double(rgb2gray(ImCrp));
    [Gx, Gy]=gradient(G);
    S=sqrt(Gx.*Gx+Gy.*Gy);
    sharpness=sum(sum(S))./(numel(Gx));
    sharp_list(k)=sharpness;
    %image = bwperim(image3);
    
end

ImCrp = imagem.*repmat(uint8(image),[1,1,3]);
imshow(ImCrp);
%imshow(rgb2gray(imagem).*uint8(image));

order_list=sort(sharp_list,'descend');

for k=1: num
    i = (find(sharp_list == order_list(k)));
    texto =  strcat(strcat(num2str(k),' - '), num2str(sharp_list(i)));
    text(stats(k).Centroid(1),stats(k).Centroid(2),texto,'HorizontalAlignment','center','Color','r');

end

%[lbw nume] = bwlabel( double(rgb2gray(imagem)));

hold off;



figure, imshow(ImCrp);hold on;
hold on
[yCenter, xCenter] = ginput(1);
hold off

%
b = bwboundaries(lb);
%
[L, num_Obj] = bwlabel(lb, 4);
%
for k = 1:num_Obj;
    Obj = L ==k;

      bb = b{k};

      X_obj = bb(:, 1);
      Y_obj = bb(:, 2);

      Selec{k} = inpolygon(xCenter,yCenter,X_obj,Y_obj); 
  end

Selec;
%
Selec = [Selec{:}];
[value,index] = max(Selec);
%
Obj = L ==index;
%
figure, imshow(Obj);

areas_list = [];
pert_list=[];
[lbw nume] = bwlabel(Obj);
imshow(lbw);
imshow(lbw);
stats1 = regionprops(lbw,'Perimeter','Area','Centroid','BoundingBox', 'ConvexHull');
Areas1 = [stats1.Area];
Perimetros1=[stats1.Perimeter]; 

figure; imshow(lb);title('Área Close to Object'); hold on;

for k=1: num
    areas_list(k) = abs( Areas1 - Areas(k));
    
end

areasord = sort(areas_list,'ascend');

for k=1: num
    if all(areas(k)>200)
        %i = find(perisord == areas(k));
        i = find(areasord == areas_list(k));
        texto =  strcat(strcat(num2str(i),' - '), num2str(areasord(i)));
        text(stats(k).Centroid(1),stats(k).Centroid(2),texto,'HorizontalAlignment','center','Color','r');
    end

end

hold off;

figure; imshow(lb);title('Perimetro Close to Object'); hold on;

for k=1: num
    pert_list(k) = abs( Perimetros1 - Perimetros(k));
    
end

pertsord = sort(pert_list,'ascend');

for k=1: num
    if all(areas(k)>200)
        %i = find(perisord == areas(k));
        i = find(pertsord == pert_list(k));
        texto =  strcat(strcat(num2str(i),' - '), num2str(pertsord(i)));
        text(stats(k).Centroid(1),stats(k).Centroid(2),texto,'HorizontalAlignment','center','Color','r');
    end

end

hold off;
