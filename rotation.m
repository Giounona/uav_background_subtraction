%this function applies affine transformation to the imput image

%input: prev_im=the image  to be transformed
%        Him, Wim = width and hight of the image
%        dx,dy= displacement in x and y directions
%        rot_ang=gotation angle
%        r_h=scaling factor
%output: newIm=transformed image

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015


function [newIm]=rotation(prev_im, rot_ang, Him, Wim ,dx,dy,r_h)

 newIm=zeros(Him,Wim,3);
 newIm1=zeros(Him,Wim);
 newIm2=zeros(Him,Wim);
 newIm3=zeros(Him,Wim);

      for x=1:Wim
          for y=1:Him
              
             x2=x-(Wim+1)/2;
             y2=y-(Him+1)/2;
             oldx=round(r_h*x2*cos(rot_ang)-r_h*y2*sin(rot_ang)+(Wim+1)/2)+(round(dx));
             oldy=round(r_h*y2*cos(rot_ang)+r_h*x2*sin(rot_ang)+(Him+1)/2)+(round(dy));

             if ~(oldx<1 || oldy<1 || oldx>Wim || oldy>Him)
             newIm1(y,x)=prev_im(oldy,oldx,1);
             newIm2(y,x)=prev_im(oldy,oldx,2);
             newIm3(y,x)=prev_im(oldy,oldx,3);
             
             end
          end
      end
      
newIm(:,:,1)=newIm1;
newIm(:,:,2)=newIm2;
newIm(:,:,3)=newIm3;
end
           
