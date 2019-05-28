function RMI = MI_EMPCA(matA,matB)
% function M = MI_GG(X,Y)
% Compute the mutual information of two images: X and Y, having
% integer values.
% 
% INPUT:
% X --> first image 
% Y --> second image (same size of X)
%
% OUTPUT:
% M --> mutual information of X and Y
%
% Written by GIANGREGORIO Generoso. 
% DATE: 04/05/2012
% E-MAIL: ggiangre@unisannio.it
%__________________________________________________________________________

  r =1;
  if size(matA)~= size(matB)
       disp('Matrices passed to Region MI must be of the same size.');
       quit;
  end
    
   BA = im2col(matA,[3 3]);
   BB = im2col(matB,[3 3]);
   [evecA,evalA] = empca(BA',1);
   [evecB,evalB] = empca(BB',1);
   
   P=[evecA';evecB'];
   [d,N]=size(P);


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
    
    
% X = double(evecA(:,1));
% Y = double(evecB(:,1));
% 
% X_norm = X - min(X(:)) +1; 
% Y_norm = Y - min(Y(:)) +1;
% 
% matAB(:,1) = X_norm(:);
% matAB(:,2) = Y_norm(:);
% h = accumarray(matAB+1, 1); % joint histogram
% 
% hn = h./sum(h(:)); % normalized joint histogram
% y_marg=sum(hn,1); 
% x_marg=sum(hn,2);
% 
% Hy = - sum(y_marg.*log2(y_marg + (y_marg == 0))); % Entropy of Y
% Hx = - sum(x_marg.*log2(x_marg + (x_marg == 0))); % Entropy of X
% 
% arg_xy2 = hn.*(log2(hn+(hn==0)));
% h_xy = sum(-arg_xy2(:)); % joint entropy
% M = Hx + Hy - h_xy; % mutual information