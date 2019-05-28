%This function finds the translation that offers the maximum mutual
%information

%input: im1,im2=the images ot be matched
%       r_c=ratio Wm/Wim;
%output: tx,ty=the translation that offers the maximum mutual information
%        maxcor= the maximum mutual information achieved
%        finalim= the final translated image

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015


function [ tx,ty,maxcor, finalim] = findMaxInf(im1,im2,r_c)


shift=round(0.7/r_c);
if shift>30
    shift=4;
end

%uncomment the following lines to apply mutual information on gradient images
% im1=rgb2gray(im1);
% im2=rgb2gray(im2);

% [fx,fy]=gradient(im1);
% im1=sqrt(fx.^2+fy.^2);
% [fx,fy]=gradient(im2);
% im2=sqrt(fx.^2+fy.^2);
% 


[m,n,l]=size(im2);
maxcor=0;
for i=-shift:1:shift
    for j=-shift:1:shift
        %im_translated=circshift(im1,[i,j]);
        [im_translated]=translationCorrection(im1, m, n ,j,i);
        [M] = MI_GG(uint8(im_translated*255),uint8(im2*255));  %mutual information method
        %[M]=region_mutual_information(im_translated,im2);%region mutual information method
        %M = MI_EMPCA(rgb2gray(im_translated), rgb2gray(im2));%mutual  information with EMPCA
        if maxcor<M         
           maxcor=M;
           finalim=im_translated;
           tx=j;
           ty=i; 
        end  
%         figure(55), imshow(c_diff)
%         pause(0.5)
    end
end

