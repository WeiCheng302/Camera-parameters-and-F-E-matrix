function [ output_args ] = Fundamental_matrix_students(im_name1, im_name2)
%UNTITLED Summary of this function goes here
%   template for assignment 2b
%   Unoptimized code, however easy to read...
%   Illustrates F-matrix

close all
format short g
skip=0;

%Options how to start: 
%>>Fundamental_matrix_students 
%this command automatically opens
%images, but you need to measure by yourself
%>>Fundamental_matrix_students(1)
%this command goes to case 2 and it automatically opens
%images, and gives you pre-measured (inaccurate) image measurements
%use this until your code works
%>>Fundamental_matrix_students('yourimage1,jpg','yourimage2.jpg')
%you can use your own images and make measurements by yourself
%notice that if you don't have circular targets, you need to modify code

switch nargin
    case 1 %for pre-measured demo
        skip=1;
        im1=imread('left.jpg');
        im2=imread('right.jpg');
    case 2 %for any images
        'input detected'
        im1=imread(im_name1);
        im2=imread(im_name2);
    otherwise %for pre-selected images
        im1=imread('left.jpg');
        im2=imread('right.jpg');
end

%open images in one figure using subplots
figure(1)
im_size=size(im1)
subplot(1,2,1)
imshow(im1);
subplot(1,2,2)
imshow(im2);

%initialize
t_points1=zeros(8,3);
t_points2=zeros(8,3);
%collect 8 tie points

if skip==0 %skip image measurements, if pre-measured option is selected
    for i=1:8 %measure from images
        'make correspondin point measurement - one measurement: left first then right'
        'adjust zoom and press space when ready to measure'
        pause %enables zooming
        subplot(1,2,1)

        %left image observation
        %[x,y]=ginput(1) %mouse input %uncomment this and comment the next
        %line to get rid of automatic measurements of a circle
        [x y]=center_of_target(im1,30,40);
        t_points1(i,:)=[x y 1];
        hold on
        filledCircle([x,y],3,3000,'b')
        hold off
        subplot(1,2,2)
        %right image observation
        %[x,y]=ginput(1) %uncomment this and comment the next
        %line to get rid of automatic measurements of a circle
        [x y]=center_of_target(im2,30,40);
        t_points2(i,:)=[x y 1];
        hold on
        filledCircle([x,y],3,3000,'r')
        hold off
    end
else %demo for images left.jpg right.jpg (inaccurate measurements i.e. 
     %give a wrong location of epipolar point, but works to test code
     %quickly
    t_points1=[945.74       536.03            1;
       1998.8       645.59            1;
       3211.6       1094.8            1;
       3182.4       1861.2            1;
       1453.4       2174.3            1;
       816.07       2278.8            1;
       774.13       1103.2            1;
         1451       1142.5            1];
     t_points2=[1055.9       582.61            1;
       1970.6       404.04            1;
       3242.8       520.57            1;
       3244.3       1402.7            1;
       1503.7       2079.8            1;
       968.01       2276.8            1;
       924.96         1158            1;
       1437.9       1038.7            1]
end
  'original observations'
  t_points1  
  t_points2
  
  %normalize image observations
  'normalized image observations and transformantion metrices'
  [n_t_points1 T1]=normalize_image_observations(t_points1)
  [n_t_points2 T2]=normalize_image_observations(t_points2)
 
 
  
  %Make design matrix A from x2'Fx1=0
  %============================================================
  nof_obs=size(n_t_points1);
  A=zeros(nof_obs(1),9);
  for i=1:nof_obs(1)
      x1=n_t_points1(i,1);
      x2=n_t_points2(i,1);
      y1=n_t_points1(i,2);
      y2=n_t_points2(i,2);
      %TODO: build the design matrix A according to lecture 3, slide 25 
      A(i,:)=[x2*x1,x2*y1,x2,y2*x1,y2*y1,y2,x1,y1,1];
  end
  A
  %A
  A2=zeros(8,8)
  for i=1:8
     A2(i,:)= [A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),A(i,8)];
  end
  A2
  %TODO: Complete the SVD solution of homogeneous system of equations
  %lecture 3, slides from 39 to 43
  L=[1,1,1,1,1,1,1,1]';
  X=inv(A2'*A2)*A2'*L;
  X(9,1)=1;
  X
  E=[X(1,1),X(2,1),X(3,1);X(4,1),X(5,1),X(6,1);X(7,1),X(8,1),X(9,1)]
  [U,V,D]=svd(A)
  %TODO:solve parametars of fundamental matrix f11-f33 with SVD
  %TODO:solution vector is the last column of V
  F=[D(1,9),D(2,9),D(3,9);D(4,9),D(5,9),D(6,9);D(7,9),D(8,9),D(9,9)]
  %TODO:add results (from the solution vector) to F-matrix ("reshape" also works...)
   
  'normalized F-matrix'
  %TODO:ensure the condition: det(F)=0 by setting s3=0
  [U1,S1,V1]=svd(F)
  S1(3,3)=0  
  %TODO:normalized Fundamental matrix
  F2=U1*S1*V1'
  
  %denormalize F-matrix
  'denormalized F-matrix'
  F=T2'*F2*T1
  
%==========================================================
  for i=1:10
    'mark point from right image and see epipolar line on left image'
    'zoom if needed - then press "space"'
    pause
    subplot(1,2,2)
    [x,y]=ginput(1)
    %t_pointsR(i,:)=[x y 1];
    hold on
    filledCircle([x,y],3,3000,'b')
    hold off
    obs_right=[x y 1]';
    epipolarline=F'*obs_right
    e_line_final=[-epipolarline(1)/epipolarline(2) -epipolarline(3)/epipolarline(2)]  
    subplot(1,2,1)
    hold on
    x=1:im_size(2);
    y=e_line_final(1)*x+e_line_final(2);
    plot(x,y, 'g-')
    hold off 
  end
  

end

function [ x2,y2 ] = center_of_target( image,tolerance,distance)
%for color images
%   Detailed explanation goes here
%image=imread(image_n);

[x,y]=ginput(1);
x=round(x) %ensure integer
y=round(y) %ensure integer
%create check_matrix
seed=image(y,x,:);

c_matrix=zeros(distance*2+1+2); %for checking neighborhood
cx=distance+1+1; %center of c_matrix
cy=distance+1+1; %center of c_matrix
%mark centerpoint
c_matrix(cx,cy)=1;

x2=0;
y2=0;
x_n=0;
y_n=0;
for s=1:distance %limit search for threshold distance
  for i=-s:s
    for j=-s:s        
        if c_matrix(j+cy,i+cx)==0 %if not already accepted
          if sum(sum(c_matrix(j+cy-1:j+cy+1,i+cx-1:i+cx+1)))>0 %if there is selected neighbor
            if image(y+j,x+i,1)>seed(1)-tolerance&&image(y+j,x+i,2)>seed(2)-tolerance&&image(y+j,x+i,3)>seed(3)-tolerance %compare with gray values of seed point
                c_matrix(j+cy,i+cx)=1; %update check matrix
                y2=y2+y+j; %accumulate coordinate
                y_n=y_n+1; %number of observations
                x2=x2+x+i; %accumulate coordinate
                x_n=x_n+1; %number of observations
            end
          end
        end
    end
  end
end
x2=x2/x_n; %average
y2=y2/y_n; %average

end

function h = filledCircle(center,r,N,color)
%---------------------------------------------------------------------------------------------
% FILLEDCIRCLE Filled circle drawing
% 
% filledCircle(CENTER,R,N,COLOR) draws a circle filled with COLOR that 
% has CENTER as its center and R as its radius, by using N points on the 
% periphery.
%
% Usage Examples,
%
% filledCircle([1,3],3,1000,'b'); 
% filledCircle([2,4],2,1000,'r');
%
% Sadik Hava <sadik.hava@gmail.com>
% May, 2010
%
% Inspired by: circle.m [Author: Zhenhai Wang]
%---------------------------------------------------------------------------------------------

THETA=linspace(0,2*pi,N);
RHO=ones(1,N)*r;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
h=fill(X,Y,color);
end