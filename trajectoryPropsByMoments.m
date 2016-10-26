function [xc,yc,theta,width,length,mp20,mp02,mp11]=trajectoryPropsByMoments(x,y);
%%

%x=squeeze(resnew(index,1));
%y=squeeze(resnew(index,2));
%area of the pattern 
[mu00,xc,yc]=trajectorymoments(x,y,0,0);
%centroid position in x
%[mu10,x0,y0]=trajectorymoments(x,y,1,0);
%xc=mu10./mu00;
%centroid position in y
%[mu01,x0,y0]=trajectorymoments(x,y,0,1);
%yc=mu01./mu00;

%define m11 m02 m20 to define theta, 
%width and length of the pattern
[mu11,x0,y0]=trajectorymoments(x,y,1,1);
[mu02,x0,y0]=trajectorymoments(x,y,0,2);
[mu20,x0,y0]=trajectorymoments(x,y,2,0);
% 

mp20=mu20./mu00;
mp02=mu02./mu00;
mp11=mu11./mu00;


% one way

  a=mp20-(xc.^2);
  b=2*(mp11-(xc.*yc));
  c=mp02-(yc.^2);
  theta=0.5*atan2(b,(a-c));
  width=(6*(a+c-((b.^2)+(a-c).^2).^0.5)).^0.5;
  length=(6*(a+c+((b.^2)+(a-c).^2).^0.5)).^0.5;
% 
% 
% a=mp20-(xc.^2);
% b=(mp11-(xc.*yc));
% c=mp02-(yc.^2);
% theta=0.5*atan2(2*b,(a-c));
% width=(0.5*((a+c)-(4*b^2+(a-c)^2)^0.5))^0.5;
% length=(0.5*((a+c)+(4*b^2+(a-c)^2)^0.5))^0.5;
% 



