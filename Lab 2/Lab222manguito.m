I = imread('eight.tif');
figure, imshow(I); hold on,
c = [222 272 300 270 221 194];
r = [21 21 75 121 121 75];

c = [ c c(1)]
r = [ r r(1)]

plot(c,r,'*b-');

BW = roipoly(I,c,r);
figure;imshow(BW);

ImCrp = I.*uint8(BW);

figure;imshow(ImCrp);