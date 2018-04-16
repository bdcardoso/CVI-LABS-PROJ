clear all, close all
img = ones(300,300);
figure;
imshow(img); axis on; hold on;


top = [ 50  55  60 65 70 75 80 
        100 100 100 100 100 100 100];

L_side = [ 50 50
            105 110];
        
R_side = [ 80 80
            105 110];

troncoL = [ 60 60 60 60 60 60
    110 115 120 125 130 135 ]

troncoR = [ 70 70 70 70 70 70
    110 115 120 125 130 135 ]

shape= 50+[top R_side troncoR troncoL L_side]