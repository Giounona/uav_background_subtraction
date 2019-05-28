%this function calculates the PSNR metric for two images

%input: curr_im, transformed_im=the image  to be matched
%        Him, Wim = width and hight of the image
%output: PSNR

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015



function [ PSNR ] = psnr_calc( curr_im, transformed_im,Wim,Him)

%      transformed_im(transformed_im==0)=curr_im((transformed_im==0)); 
%      c_diff=(rgb2gray(uint8(curr_im.*255))-rgb2gray(uint8(transformed_im*255))).^2;                
%      MSerror=sum(sum(c_diff))/(Wim*Him);
%      PSNR=10*log(255^2/MSerror);
%      
     transformed_im(transformed_im==0)=curr_im((transformed_im==0)); 
     c_diff=(rgb2gray((curr_im))-rgb2gray((transformed_im))).^2;                
     MSerror=sum(sum(c_diff))/(Wim*Him);
     PSNR=10*log(1/MSerror);


end

