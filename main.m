
% run this m-file to calculate the volumes of intersecting 3D circles

dimensionsofanalyses = 50; % the centers of the circles will have coordinates 
                           % between x= 0 - 50m, y= 0 - 50m, z= 0 - 50m
numberofcircles =  1000; % 100000;
meanradiusofcircles = 5; % 5;

radiiofcircles = exprnd(meanradiusofcircles, numberofcircles,1); %radii of all circles
xofcircles = dimensionsofanalyses*rand(numberofcircles,1); %x-coordinates of the centres of all circles
yofcircles = dimensionsofanalyses*rand(numberofcircles,1); %y-coordinates of the centres of all circles
zofcircles = dimensionsofanalyses*rand(numberofcircles,1); %z-coordinates of the centres of all circles

xnormalvector = -1+rand(numberofcircles,1)*2; %x-vector of normal vector (normal to the circle plane)
ynormalvector = -1+rand(numberofcircles,1)*2; %y-vector of normal vector (normal to the circle plane)
znormalvector =  0+rand(numberofcircles,1)*1; %z-vector of normal vector (normal to the circle plane) 

n_ = numberofcircles;
R1 = radiiofcircles(:);
C1 = [xofcircles(:) yofcircles(:) zofcircles(:)];
N1 = [xnormalvector(:) ynormalvector(:) znormalvector(:)];
N1 = N1 ./ repmat(sqrt(dot(N1,N1,2)),1,3);

[R,idx] = sort(R1,'ascend');
% [R,idx] = sort(R1,'descend');
C = C1(idx,:);
N = N1(idx,:);

bPair = crossPointsOfCircles3D(C, R, N);

