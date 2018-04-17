I = imread('Moedas1.jpg');
figure,imshow(I) 
axis off
title('RGB Image');
threshRGB = multithresh(I,7);
threshForPlanes = zeros(3,7);			

for i = 1:3
    threshForPlanes(i,:) = multithresh(I(:,:,i),7);
end
value = [0 threshRGB(2:end) 255]; 
quantRGB = imquantize(I, threshRGB, value);
quantPlane = zeros( size(I) );

for i = 1:3
    value = [0 threshForPlanes(i,2:end) 255]; 
    quantPlane(:,:,i) = imquantize(I(:,:,i),threshForPlanes(i,:),value);
end

quantPlane = uint8(quantPlane);
imshowpair(quantRGB,quantPlane,'montage') 
axis off
title('Full RGB Image Quantization        Plane-by-Plane Quantization')
dim = size( quantRGB );
quantRGBmx3   = reshape(quantRGB,   prod(dim(1:2)), 3);
quantPlanemx3 = reshape(quantPlane, prod(dim(1:2)), 3);

colorsRGB   = unique(quantRGBmx3,   'rows' );
colorsPlane = unique(quantPlanemx3, 'rows' );

disp(['Unique colors in RGB image            : ' int2str(length(colorsRGB))]);

disp(['Unique colors in Plane-by-Plane image : ' int2str(length(colorsPlane))]);
