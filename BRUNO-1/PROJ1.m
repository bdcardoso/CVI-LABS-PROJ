clc 
clear all
close all

a = imread('Moedas4.jpg');

figure, imshow(a), title('First import')

b = rgb2gray(a);

figure, imshow(b), title('now in grayscale')
c = 255-b

figure, imshow(c);title('Now original in B')

d= imbinarize(c)
%d = im2bw(c);
e = imfill(d, 'holes');

figure, imshow(e), title('now filled')

f = bwlabel(e);

figure, imshow(f), title('now labelled')

vislabels(f),title('each object labeled')

g = regionprops(f, 'Area', 'BoundingBox');


g(1)

area_values = [g.Area]

idx = find((7000<= area_values) & (area_values <= 9000))

h = ismember(f, idx);

figure, imshow(h), title('Area between 7000 and 9000')

i = regionprops(f, 'Area', 'centroid');

score= (min(sqrt([i.Area]),[i.Perimeter]/4)./(max(sqrt([i.Area]),[i.Perimeter]/4))).^2;

figure, imshow(c), title('squareness of each object')
for cnt = 1:length(i)
    text(i(cnt).Centroid(1), i(cnt).Centroid(2),num2str(score(cnt)),'FontSize',15,'color','red');
end


%AUX FUNCTION

function vislabels(L)
%VISLABELS Visualize labels of connected components
%   VISLABELS is used to visualize the output of BWLABEL.
%
%   VISLABELS(L), where L is a label matrix returned by BWLABEL,
%   displays each object's label number on top of the object itself.
%
%   Note: VISLABELS requires the Image Processing Toolbox.
%
%   Example
%   -------
%       bw = imread('text.png');
%       L = bwlabel(bw);
%       vislabels(L)
%       axis([1 70 1 70])

%   Steven L. Eddins
%   Copyright 2008 The MathWorks, Inc.

% Form a grayscale image such that both the background and the
% object pixels are light shades of gray.  This is done so that the
% black text will be visible against both background and foreground
% pixels.

background_shade = 200;
foreground_shade = 240;
I = zeros(size(L), 'uint8');
I(L == 0) = background_shade;
I(L ~= 0) = foreground_shade;

% Display the image, fitting it to the size of the figure.
imageHandle = imshow(I, 'InitialMagnification', 'fit');

% Get the axes handle containing the image.  Use this handle in the
% remaining code instead of relying on gca.
axesHandle = ancestor(imageHandle, 'axes');

% Get the extrema points for each labeled object.
s = regionprops(L, 'Extrema');

% Superimpose the text label at the left-most top extremum location
% for each object.  Turn clipping on so that the text doesn't
% display past the edge of the image when zooming.
hold(axesHandle, 'on');
for k = 1:numel(s)
   e = s(k).Extrema;
   text(e(1,1), e(1,2), sprintf('%d', k), ...
      'Parent', axesHandle, ...
      'Clipping', 'on', ...
      'Color', 'b', ...
      'FontWeight', 'bold');
end
hold(axesHandle, 'off');

end





