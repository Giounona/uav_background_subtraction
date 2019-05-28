%this function extract the telementry data from the given text line

%input: InputText= one line form the telemetry text file
%       W_sensor = width of the sensor 
%       ASLg= the minimum  hight obtained from telemetry data
%output: data= structure that contains all the relevant telemetry data

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015

function [data]=getTelemetryData(InputText,W_sensor, ASLg)
  format long g
  data.framenum=str2double(cell2mat(InputText{1}(1))); %frame number
  data.Lat=str2double(cell2mat(InputText{1}(3)));%lattitude
  data.Lon=str2double(cell2mat(InputText{1}(4)));%longitude
  data.ASLobs=str2double(cell2mat(InputText{1}(5))); %hight
  data.Yaw=str2double(cell2mat(InputText{1}(9)));%yaw ungle

  data.FOV=(pi/180)*(W_sensor/2);%field of view
  data.ASLg=ASLg; %hight from ground
  
  data.Height=data.ASLobs-data.ASLg;
  data.Wm=2*data.Height*tan(data.FOV); %the width covered by the field of view in meters
  [data.e, data.n, data.zoneNumber ]=LLtoUTM(data.Lat, data.Lon); %convertion of latitude and longitude to easting and nirthing
  