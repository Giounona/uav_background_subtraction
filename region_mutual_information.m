%this function finds region mutual information between two images

%input: matA,matB=the image  to be matched
%       
%output: RMI=Region mutual information

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015


    function [RMI]=region_mutual_information(matA, matB)
    r =1;
    if size(matA)~= size(matB)
       disp('Matrices passed to Region MI must be of the same size.');
       quit;
    end
    
   
    [rows, cols ]=size( matA);
    N = (rows - 2*r) * (cols - 2*r);
    d = 2 * (2*r + 1)^2;

 %   # Given two images, A and B for each corresponding pair of pixels [Aij, Bij]
  %  #create a vector vij ... and a matrix P ...
     P = generate_P( matA, matB,r);

   % # Subtract the mean from the points so that they are centered
    %#at the origin
    px = sum(P,2) ./ (N);
    PX=px * ones(N,1)';
    %PX = (px.*N)';
    P0 =(P- PX);

   % # Calculate the covariance of the points C = 1/N P0 P0^T
    C = (1/(N)* (P0* P0'));

    %# Estimate the joint entropy as H_g(C)
    HgC =log((2*pi*exp(1))^(d/2) * det(C));
    HgCA = log((2*pi*exp(1))^(d/2) * det(C(1:d/2, 1:d/2)));
    HgCB = log((2*pi*exp(1))^(d/2) * det(C(d/2+1:end, d/2+1:end)));

    RMI = HgCA + HgCB - HgC;
    end

    function [ P ] = generate_P( matA, matB,r)

    [rows, cols] = size(matA);
     N = (rows - 2*r) * (cols - 2*r);
     d = 2 * (2*r + 1)^2;

     P = zeros(d, N);
   
    BA = im2col(matA,[3 3]);
    BB = im2col(matB,[3 3]);
    P=[BA;BB];
      
    end


