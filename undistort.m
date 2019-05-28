%this function removes lence distortion form the image
%input:Idistorted=distorted image
%output:I=undistorted image
%the calobration parameterd were obtained with camera calibration routine
%form OpenCV library

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015

function I = undistort(Idistorted)

% fx = params.fx;
% fy = params.fy;
% cx = params.cx;
% cy = params.cy;
% k1 = params.k1;
% k2 = params.k2;
% k3 = params.k3;
% p1 = params.p1;
% p2 = params.p2;

fx =  726.698443892279 ;
fy = 727.300752756491  ;
cx = 255.5;
cy =  255.5;
k1 = -0.203155439721038;
k2 = 0;
k3 = 0;
p1 = 0;
p2 = 0;


K = [fx 0 cx; 0 fy cy; 0 0 1];

I = zeros(size(Idistorted));
[i j] = find(~isnan(I));

% Xp = the xyz vals of points on the z plane
Xp = inv(K)*[j i ones(length(i),1)]';

% Now we calculate how those points distort i.e forward map them through the distortion
r2 = Xp(1,:).^2+Xp(2,:).^2;
x = Xp(1,:);
y = Xp(2,:);

x = x.*(1+k1*r2 + k2*r2.^2+k3*r2.^3) + 2*p1.*x.*y + p2*(r2 + 2*x.^2);
y = y.*(1+k1*r2 + k2*r2.^2+k3*r2.^3) + 2*p2.*x.*y + p1*(r2 + 2*y.^2);

% u and v are now the distorted cooridnates
u = reshape(fx*x + cx,size(I));
v = reshape(fy*y + cy,size(I));

% Now we perform a backward mapping in order to undistort the warped image coordinates
I = interp2(Idistorted, u, v);