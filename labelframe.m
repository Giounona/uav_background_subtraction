%this function labels the frame with foreground regions

%input: im=the frame to be labeled
%       mask=the foreground mask
%output: labeled=labeled frame

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015



function [labeled]=labelframe(im, mask)

   temp1=im(:,:,1);
   temp2=im(:,:,2);
   temp3=im(:,:,3);
   temp1(mask==1)=1;
   temp2(mask==1)=0;
   temp3(mask==1)=0;
   im(:,:,1)=temp1;
   im(:,:,2)=temp2;
   im(:,:,3)=temp3;

labeled=im;