%this function applies median filtering to matched sift points

%input:  matchLoc1, matchLoc2= initial location of SIFT points in image1
%
%output: matchLoc1, matchLoc2= final location of SIFT points in image1


%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015


function [ matchLoc1, matchLoc2 ] = filter_much_points( matchLoc1, matchLoc2)

[m,n]=size(matchLoc1);
flag=zeros(m,1);
thresh=zeros(m,1)+50;
vector=[matchLoc2(:,1)-matchLoc1(:,1), matchLoc2(:,2)-matchLoc1(:,2)];


for i=1:m
    
    x1=matchLoc1(i,1);
    y1=matchLoc1(i,2);
    d=sqrt((matchLoc1(:,1)-x1).^2+(matchLoc1(:,2)-y1).^2);
    flag(d<thresh)=1;
    if(sum(flag)>3)
    [idx,val]=find(flag==1);
    neighbour=vector(idx,:);    
    mad1=mad(neighbour);
    median1=median(neighbour);
        if abs(median1-[matchLoc2(i,1)-x1,matchLoc2(i,2)-y1])>1.4826*mad1
            dx=median1(1);
            dy=median1(2);
            matchLoc2(i,1)=x1+dx;
            matchLoc2(i,2)=y1+dy;
        end        
    end
    
    
end


end

