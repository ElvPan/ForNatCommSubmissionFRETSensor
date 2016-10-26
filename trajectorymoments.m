function [trajmom,xc,yc]=trajectorymoments(x,y,p,q);

%EP 2012
%this function computes coordintes (x,y) moments of 
%a particle trajectory for order p and q 



    m00=sum(sum((x.^0).*(y.^0)));

    xc=sum(sum((x.^1).*(y.^0)))/m00;
    yc=sum(sum((x.^0).*(y.^1)))/m00;

    trajmom=sum(sum(((x).^p).*((y).^q)));


    