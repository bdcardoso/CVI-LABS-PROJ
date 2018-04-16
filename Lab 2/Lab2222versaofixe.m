close all, clear all

A = zeros(9);
A(1:6, 1:6) = 1
B = zeros(9);
B(5:9,5:9) = 1

A = zeros(9);
A(1:6,1:6) = 1;
B = zeros(9);
B(6:9,6:9) = 1

figure;imagesc(A); colormap gray
figure;imagesc(B); colormap gray


figure;imagesc(~A); colormap gray
figure;imagesc(~B); colormap gray

C2 = A .* B;
figure;imagesc(C2); colormap gray

C3 = A | B
figure;imagesc(C3); colormap gray

C4 = A