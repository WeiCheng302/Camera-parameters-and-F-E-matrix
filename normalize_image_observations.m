%input coordinates should be 2D image observations in homogeneous
%coordinates (x,y,1)
function [normalized_points,T]=normalize_image_observations(points)

obs_size=size(points);

%average distance x and y
mean1=mean(points);
myy_x=mean1(1);
myy_y=mean1(2);

mean_dist=sqrt(points(:,1).^2+points(:,2).^2);
%average mean distance
d1=mean(mean_dist);
scale=sqrt(2)/d1;
T=[scale 0 -scale*myy_x; 0 sqrt(2)/d1 -scale*myy_y; 0 0 1];
points=points';

for i=1:obs_size(1) 
    normalized_points(i,:)=(T*points(:,i))';
end

end