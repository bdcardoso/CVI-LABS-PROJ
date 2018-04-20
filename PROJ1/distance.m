function [centx,centy,r]=distance(D);

[w h]=size(D');

mx=max(max(D));
r=mx;
for i=1:h
    for j=1:w
        if D(i,j)==mx;
            centx=j;
            centy=i;
        end
    end
end