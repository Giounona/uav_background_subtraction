%
%This function compaires the aerial video frame registration methods
%Install vlfeat before running the experiments http://www.vlfeat.org
%
%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015


clear all;
close all;
addpath ../.
cfg.input.dir ='./2014-10-07/photos/';
fid=fopen('./2014-10-07/telemetry/telemetry_data.txt','r');
cfg.input.basename = '%06d';
cfg.input.fmt = 'jpg';
cfg.output.fmt = 'png';

addpath(genpath('./third party/'))

  
%*****telemetry related data, varied according to the selected dataset
W_sensor=38.8033;
ASLg=546.091003418; %2014-10-04 
%ASLg=553.721008301;%2014-10-07
%ASLg=546.99597168;%2015
%%%%**********************************


start=400;
hop=1;
cfg.exp.idx = [start:hop:700]; % [2:10:354];   % 212 for Erlangen data, 354 for advbgst, 1-4956 ETH_EFloor
scaling = 0.5;

count=1;

format long g
for rawind=1:start
InputText = textscan(fid,'%s',9,'delimiter',';');
end    
     
counts=[];
psnrsMy=[];
psnrsSIFT=[];
psnrsHARRIS=[];
mssimsMy=[];
mssimsSIFT=[];
mssimsHARRIS=[];

for i=cfg.exp.idx,
    
  fprintf('processing frame: %3d \b\n',i)
  im = imread([cfg.input.dir,sprintf(cfg.input.basename,i),'.',cfg.input.fmt]);% read from file
  im=double(im)/255;
  Iundist(:,:,1) = undistort(im(:,:,1));
  Iundist(:,:,2) = undistort(im(:,:,2));
  Iundist(:,:,3) = undistort(im(:,:,3));
  im=Iundist;

  im=imresize(im,scaling);
  [Him,Wim,depth]=size(im);
  im=imrotate(im,180);
 
  
  if count==1;
      InputText = textscan(fid,'%s',9,'delimiter',';');  
      curr_data=getTelemetryData(InputText,W_sensor, ASLg); 
      curr_im=im;
      accum_im=curr_im;
      rot_ang=0; dx=0; dy=0; r_h=1;
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
      
      de=(curr_data.e-prev_data.e)/(r_c);
      dn=(curr_data.n-prev_data.n)/(r_c); 
      
      Yawdiff=curr_data.Yaw-prev_data.Yaw;
      
%       angle=(pi/180)*(360-prev_data.Yaw);%anticlocwise rotation
%       dx=(de*cos(angle)+dn*sin(angle));
%       dy=(-de*sin(angle)+dn*cos(angle));   
      
      angle=(pi/180)*(prev_data.Yaw);%clocwise rotation
      dy=(dn*cos(angle)+de*sin(angle));
      dx=(-dn*sin(angle)+de*cos(angle));   
       
      rot_ang=(pi/180)*(sign(Yawdiff)*min(abs(Yawdiff),360-abs(Yawdiff)));      


      %[accum_im]=rotation(accum_im, rot_ang, Him, Wim ,dx,-dy,r_h); 
      af_mat=[ r_h*cos(rot_ang)   r_h*sin(rot_ang) -round(dy);
                -r_h*sin(rot_ang)   r_h*cos(rot_ang)  round(dx) ;
                 0                    0                        1];             
      transformed_imA=affine_warp(prev_im,af_mat,'nearest');

      [fx,fy]=gradient(rgb2gray(curr_im));
      curr_grad=sqrt(fx.^2+fy.^2);
      [fx,fy]=gradient(rgb2gray(transformed_imA));
      transf_grad=sqrt(fx.^2+fy.^2);
   
      [outputA GregA] =dftregistration(fft2(curr_grad),fft2(transf_grad));
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

      af_mat3=[ r_h*cos(rot_ang)   r_h*sin(rot_ang) ty;
                -r_h*sin(rot_ang)   r_h*cos(rot_ang)  tx ;
                 0                    0                        1];
      imMy=affine_warp(prev_im,af_mat3,'bicubic');     
      % [ tx,ty,minerr, finalim] = findMaxInf(transformed_imA,curr_im,r_c);
    
 %********************SIFT features**************************************
      [matchLoc1 matchLoc2] = siftMatch(uint8(curr_im.*255), uint8(prev_im.*255));
      %median filtering of matched sift points
      [ matchLoc1, matchLoc2 ] = filter_much_points( matchLoc1, matchLoc2);
      [H corrPtIdx] = findHomography(matchLoc2',matchLoc1');
      tform = maketform('projective', H');
      imSIFT = imtransform(prev_im,tform,'bicubic',...
                                  'XData',[1 Wim], 'YData',[1 Him],'UData' , [0 Wim] , 'VData' , [0 Him]);
      [fx,fy]=gradient(rgb2gray(imSIFT));
      transf_grad=sqrt(fx.^2+fy.^2);
      [outputS GregA] =dftregistration(fft2(curr_grad),fft2(transf_grad));
      [imSIFT]=translationCorrection(imSIFT,Him, Wim ,-outputS(4),-outputS(3));
 %************************************Harris features*********************
      [ muchloc1, muchloc2 ] = harris_match( curr_im, prev_im);
      [H3 corrPtIdx] = findHomography(muchloc2',muchloc1');
      tform = maketform('projective', H3');
      imHARRIS = imtransform(prev_im,tform,'bicubic',...
                                  'XData',[1 Wim], 'YData',[1 Him],'UData' , [0 Wim] , 'VData' , [0 Him]);
      [fx,fy]=gradient(rgb2gray(imHARRIS));
      transf_grad=sqrt(fx.^2+fy.^2);
      [outputH GregA] =dftregistration(fft2(curr_grad),fft2(transf_grad));
      [imHARRIS]=translationCorrection(imHARRIS,Him, Wim ,-outputH(4),-outputH(3));
                              
       
%      figure(1), imshow(curr_im);
%      figure(2), imshow(imHARRIS);
%      figure(4), imshow(imMy);
%      figure(5), imshow(imSIFT);
   
   %*************calculation of MSSIM and PSNR metrics*******************  
      [mssimMy, ssim_map] = ssim_index(uint8(rgb2gray(curr_im)*255), uint8(rgb2gray(imMy)*255));
      [mssimHARRIS, ssim_map] = ssim_index(uint8(rgb2gray(curr_im)*255), uint8(rgb2gray(imHARRIS)*255));
      [mssimSIFT, ssim_map] = ssim_index(uint8(rgb2gray(curr_im)*255), uint8(rgb2gray(imSIFT)*255));
     
      [ psnrSIFT ] = psnr_calc( curr_im,imSIFT,Wim,Him);
      [ psnrMY ] = psnr_calc( curr_im, imMy,Wim,Him);
      [ psnrimHARRIS ] = psnr_calc( curr_im,  imHARRIS,Wim,Him);
      counts=[counts; count];


      psnrsMy=[psnrsMy,psnrMY];
      psnrsSIFT=[psnrsSIFT,psnrSIFT];
      psnrsHARRIS=[psnrsHARRIS,psnrimHARRIS];
      mssimsMy=[mssimsMy,mssimMy];
      mssimsSIFT=[mssimsSIFT,mssimSIFT];
      mssimsHARRIS=[mssimsHARRIS,mssimHARRIS];
    
   %plot the results
      if i==start+200 
         
        figure(1),plot(counts, mssimsSIFT,'ro','LineWidth',2),hold on,...
            plot(counts, mssimsMy,'c+','LineWidth',2),hold on, plot(counts,mssimsHARRIS,'g*','LineWidth',2)
        
        figure(2),plot(counts,psnrsSIFT,'ro','LineWidth',2),hold on,...
            plot(counts,psnrsMy,'c+','LineWidth',2),plot(counts,psnrsHARRIS,'g*','LineWidth',2)
         
      end
      
    %    [accum_im]=rotation(accum_im, rot_ang, Him, Wim ,tx,ty,r_h);  
    %    accum_im(accum_im==0)=im((accum_im==0)); 
    %    h=figure(4);subplot(1,3,1), imshow((rgb2gray(curr_im)-rgb2gray(imMy))), title(num2str(i))
    %    subplot(1,3,2),imshow(accum_im);
    %    subplot(1,3,3),imshow((imMy)), title(num2str(i))
    %    saveas(h,['./2014-10-07/results/', num2str(i),'f.jpg']);

  

 end
 count=count+1;
  
end

