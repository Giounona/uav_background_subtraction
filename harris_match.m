%this function finds the matching harris points between two images

%input: img1, img2=the imageg  to be matched
% 
%output: matchLoc1 matchLoc2=the location of matching points in the two
%images

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015



function [matchLoc1 matchLoc2] = harris_match(img1, img2)

%uses vlFeat to calculate harris


img1s = im2single(rgb2gray(img1));
[loc1,des1] = vl_covdet(img1s, 'method', 'MultiscaleHarris','PeakThreshold',0.0000001);
     
img2s = im2single(rgb2gray(img2));
[loc2,des2] = vl_covdet(img2s, 'method', 'MultiscaleHarris','PeakThreshold',0.0000001);

des1=des1';
des2=des2';
loc1=loc1';
loc2=loc2';


% For efficiency in Matlab, it is cheaper to compute dot products between
%  unit vectors rather than Euclidean distances.  Note that the ratio of 
%  angles (acos of dot products of unit vectors) is a close approximation
%  to the ratio of Euclidean distances for small angles.
%
% distRatio: Only keep matches in which the ratio of vector angles from the
%   nearest to second nearest neighbor is less than distRatio.
distRatio = 0.6;   

% For each descriptor in the first image, select its match to second image.
des2t = des2';                          % Precompute matrix transpose
matchTable = zeros(1,size(des1,1));
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;        % Computes vector of dot products
   [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

   % Check if nearest neighbor has angle less than distRatio times 2nd.
   if (vals(1) < distRatio * vals(2))
      matchTable(i) = indx(1);
   else
      matchTable(i) = 0;
   end
end
% save matchdata matchTable
%}

% Create a new image showing the two images side by side.
img3 = appendimages(img1,img2);
img3=uint8(img3.*255);
%Show a figure with lines joining the accepted matches.
% figure('Position', [100 100 size(img3,2) size(img3,1)]);
% colormap('gray');
% imagesc(img3);
% hold on;
% cols1 = size(img1,2);
% for i = 1: size(des1,1)
%   if (matchTable(i) > 0)
%     line([loc1(i,1) loc2(matchTable(i),1)+cols1], ...
%          [loc1(i,2) loc2(matchTable(i),2)], 'Color', 'c');
%   end
% end
% hold off;
 num = sum(matchTable > 0);
fprintf('Found %d matches.\n', num);

idx1 = find(matchTable);
idx2 = matchTable(idx1);
x1 = loc1(idx1,1);
x2 = loc2(idx2,1);
y1 = loc1(idx1,2);
y2 = loc2(idx2,2);

matchLoc1 = [x1,y1];
matchLoc2 = [x2,y2];

end