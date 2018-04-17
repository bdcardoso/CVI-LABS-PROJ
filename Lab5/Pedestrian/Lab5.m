clear all;
close all;

thr = 50;
minArea = 35;
baseNum = 1350;
seqLength = 100;

imgbk = imread('ped7c0000.tif');

se = strel ('disk',3);

figure; hold on;
for i=0:seqLength
    imgfr = imread(sprintf('ped7c%.4d.tif', baseNum+i));
    imshow(imgfr); %hold on;
    
    imgdif = ...
        (abs ( double (imgbk(:,:,1))-double(imgfr(:,:,1)))>thr) | ...
        (abs ( double (imgbk(:,:,2))-double(imgfr(:,:,2)))>thr) | ...
        (abs ( double (imgbk(:,:,3))-double(imgfr(:,:,3)))>thr);
    
    bw = imclose(imgdif,se);
    %imshow(bw);
    

    
   [lb num] = bwlabel(bw);
   regionProps = regionprops(lb, 'area', 'FilledImage', 'Centroid');
   inds = find([regionProps.Area]> minArea);
   
   %hold on;
   
   regnum = length(inds);
   
   
   %for i=1:length(inds)
   for j=1:regnum
       %props = regionprops(double(regionProps(inds(i)).FilledImage),...
       %    'Orientation','MajorAxisLength','MinorAxisLength');
       %ellipse(props.MajorAxisLength/2,props.MinorAxisLength/2,...
       %    -props.Orientation*pi/180,...
       %    regionProps(inds(i)).Centroid(1),...
       %    regionProps(inds(i)).Centroid(2),'r');
       
       [lin col] = find(lb == inds(j));
       upLPoint = min ([lin col]);
       dWindow = max([lin col]) - upLPoint + 1;
       rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)], 'EdgeColor', [1 1 0], ...
           'linewidth', 2);
       
%        if exist('propsT')
%            propsT = [propsT props];
%        else
%            propsT = props;
%        end
   end
   drawnow;
end
       
       
       
       
       
       
   
    
    
