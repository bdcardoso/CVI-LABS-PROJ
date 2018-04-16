clear all, close all

imgg = imread('veiculoGray.jpg');
imshow(imgg); hold on;
N=0; but=1;

while (but==1 | but == 32)
    [ci, li, but]=ginput(1);
    but
    
    if but==1 %add point
        N = N+1;
        cp(N) = ci;
        lp(N)= li;
        plot(ci, li, 'r.', 'MarkerSize', 18); drawnow;
        if N > 1
            plot(cp(:), lp(:), 'r.-', 'MarkerSize', 8); drawnow;
        end
    end
end
cp=cp'; lp = lp';
BW = roipoly(imgg, cp, lp);

ImCrp = imgg.*uint8(BW);
imagesc(ImCrp); colormap gray

%pause


