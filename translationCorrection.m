%this function translates the image by "dx" and "dy" in x and y directions
%respectively

%input: prev_im=the image  to be translated
%        Him, Wim = width and hight of the image
%          dx,dy= displacement in x and y directions
%output: newIm=translated image

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015

function [newIm]=translationCorrection(prev_im,Him, Wim ,dx,dy)

[m,n,l]=size(prev_im);

if l==3
 newIm=zeros(Him,Wim,3);
 tempIm=zeros(Him+1,Wim+1,3);
 tempIm(1:Him,1:Wim,:)=prev_im;
 x=1:Wim;
 y=1:Him;
 oldx=x+(round(dx));
 oldy=y+(round(dy));
 oldx(oldx<1)=Wim+1;
 oldy(oldy<1)=Him+1;
 oldx(oldx>Wim)=Wim+1;
 oldy(oldy>Him)=Him+1;
 
 newIm(y,x,1)=tempIm(oldy,oldx,1);
 newIm(y,x,2)=tempIm(oldy,oldx,2);
 newIm(y,x,3)=tempIm(oldy,oldx,3);
else
 newIm=zeros(Him,Wim);
 tempIm=zeros(Him+1,Wim+1);
 tempIm(1:Him,1:Wim)=prev_im;
 x=1:Wim;
 y=1:Him;
 oldx=x+(round(dx));
 oldy=y+(round(dy));
 oldx(oldx<1)=Wim+1;
 oldy(oldy<1)=Him+1;
 oldx(oldx>Wim)=Wim+1;
 oldy(oldy>Him)=Him+1;
 
 newIm(y,x)=tempIm(oldy,oldx);

    
end

       
end