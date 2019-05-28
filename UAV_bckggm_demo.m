%this is the main fuction ato run the background subtraction for UAV

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015


clear all;
close all;
addpath ../.
cfg.input.dir ='./2014-10-07/photos/';
cfg.output.dir='./2014-10-07/results/';
fid=fopen('./2014-10-07/telemetry/telemetry_data.txt','r');
maskim=rgb2gray(imread('./2014-10-07/mask.jpg'));
addpath(genpath('./third party/'))

%*****telemetry related data varied according to the selected dataset
W_sensor=38.8033;
ASLg=546.091003418; %2014-10-04 
%ASLg=553.721008301;%2014-10-07
%ASLg=546.99597168;%2015
%%%%**********************************

cfg.input.basename = '%06d';
cfg.input.fmt = 'jpg';
cfg.output.fmt = 'png';

  
cfg.exp.alpha = 0.07;%0.1%the higher alpha the faster the learnign rate, the lower alpha the more noise and fuller images
cfg.exp.K =4;
cfg.exp.var = 0.07;%0.06 % initial variance  
cfg.exp.minvar=cfg.exp.var*0.1;% the initial values are so small becaise we are using the normalized colour space where the values are of the range 0.a

start=400;%the starting frame
bckgdStats=0.08;
hop=1;%how many frames to skip. DO NOT INCREASE IT.
cfg.exp.idx = [start:hop:5000]; 

cfg.exp.thr = 2.5; 
cfg.exp.T = 0.6; 

model = [];
Mymodel=[];
idx = [];
addpath('GCmex1.3');
scaling = 0.5;


count=1;
format long g %this is required because for lat and long we need high precision
for rawind=1:start%reads the telemetry file
InputText = textscan(fid,'%s',9,'delimiter',';');
end    
     

level = graythresh(maskim);
BW = im2bw(maskim,level);
mask=BW;


% main background modeling cycle
for i=cfg.exp.idx,
    
  fprintf('processing frame: %3d \b\n',i)
  im = imread([cfg.input.dir,sprintf(cfg.input.basename,i),'.',cfg.input.fmt]);% read from file
  im=double(im)/255;
  Iundist(:,:,1) = undistort(im(:,:,1)); %remove lence distortion
  Iundist(:,:,2) = undistort(im(:,:,2));
  Iundist(:,:,3) = undistort(im(:,:,3));
  im=Iundist;

  im=imresize(im,scaling);
  mask=imresize(BW,scaling);
  [Him,Wim,depth]=size(im);
  im=imrotate(im,180); % this is rotated because the IMU sensor is at 180 degrees with the camera. The hardware should be corrected 
 
  if count==1;%obtain the telemetry data to initialize the background model
      InputText = textscan(fid,'%s',9,'delimiter',';');  
      curr_data=getTelemetryData(InputText,W_sensor, ASLg); 
      curr_im=im;
      rot_ang=0; dx=0; dy=0; r_h=1;
      [model,idx,res,Mymodel] = bckggm(im,model,idx,cfg.exp,rot_ang,dx,-dy,Mymodel,r_h,mask);
  else
      for j=1:hop
        InputText = textscan(fid,'%s',9,'delimiter',';');  
      end
      
      prev_data=curr_data;
      prev_im=curr_im;
      curr_im=im;
      curr_data=getTelemetryData(InputText,W_sensor, ASLg); 
      r_c=curr_data.Wm/Wim;
      r_h=(curr_data.Height/prev_data.Height);
      
      de=(curr_data.e-prev_data.e)/(r_c);%easting and northing values
      dn=(curr_data.n-prev_data.n)/(r_c); 
      
      Yawdiff=curr_data.Yaw-prev_data.Yaw;
      
%       angle=(pi/180)*(360-prev_data.Yaw);%anticlocwise rotation
%       dx=(de*cos(angle)+dn*sin(angle));
%       dy=(-de*sin(angle)+dn*cos(angle));   
      
      angle=(pi/180)*(prev_data.Yaw);%clocwise rotation
      dy=(dn*cos(angle)+de*sin(angle));%displacement over y axis
      dx=(-dn*sin(angle)+de*cos(angle));   %displacement over x axis
      rot_ang=(pi/180)*(sign(Yawdiff)*min(abs(Yawdiff),360-abs(Yawdiff)));     
      
     %affine matrix that describes the motion of the UAV
     af_mat=[ r_h*cos(rot_ang)   r_h*sin(rot_ang) -round(dy);
                -r_h*sin(rot_ang)   r_h*cos(rot_ang)  round(dx) ;
                 0                    0                        1];             

     transformed_imA=affine_warp(prev_im,af_mat,'nearest');
     
     [fx,fy]=gradient(rgb2gray(curr_im));
     curr_grad=sqrt(fx.^2+fy.^2);
     [fx,fy]=gradient(rgb2gray(transformed_imA));
     transf_grad=sqrt(fx.^2+fy.^2);
   
     [outputA GregA] =dftregistration(fft2(curr_grad),fft2(transf_grad));% global registration to correct the errors from the GPS location
     tx=outputA(4);
     ty=outputA(3);
     shift=round(2.5/r_c);
     if abs(tx)<shift && abs(ty)<shift
        tx=dx-outputA(4);
        ty=-dy-outputA(3);
     else
        tx=dx;
        ty=-dy;
     end
     
     % Updateing the backgroun model. The initial background model "model"
     % is replaced by the backgrounf model "Mymodel" which allows easier
     % application of affine transformation (in case of videos that no mask is required)
     [model,idx,res,Mymodel] = bckggm(im,model,idx,cfg.exp,rot_ang,tx,ty,Mymodel,r_h,mask);    
     

     [im]=labelframe(im, res.segm); %colour in red the foreground objects
     A=zeros(Him,3*Wim,3);
     A(:,1:Wim,:)=im;
     A(:,Wim+1:2*Wim,1)=res.segm;
     A(:,Wim+1:2*Wim,2)=res.segm;
     A(:,Wim+1:2*Wim,3)=res.segm;
     A(:,2*Wim+1:3*Wim,:)=res.bckg;
     imwrite(A,[cfg.output.dir, num2str(i),'.jpg'],'jpg','Quality',90);

 end
 count=count+1;
  
end

