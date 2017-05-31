
% run this m-file to calculate the volumes of intersecting 3D circles

addpath('testutils');

dimensionsofanalyses = 50; % the centers of the circles will have coordinates 
                           % between x= 0 - 50m, y= 0 - 50m, z= 0 - 50m
numberofcircles =  100; % 100000;
meanradiusofcircles = 10; % 5;

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

save('testutils\matlab100.mat');

% calculate cross points of circles
[bPair, ftriname, fptsname] = crossPointsOfCircles3D(C, R, N);

% display results
if (n_ > 1000), return; end

% read cross points info
fcrosspts = 'crosspts.dat';
fpt = fopen(fcrosspts, 'r');
CPT = fread(fpt, [6,Inf], 'float');
CPT = CPT(1:3,:)';
fseek(fpt, 16, 'bof');
IDX = fread(fpt, Inf, 'int64', 16);
IDX = IDX + 1;
fclose(fpt);

% read triples info
ft = fopen(ftriname, 'r');
fseek(ft, 0, 'eof');
flen = ftell(ft);
fseek(ft, 0, 'bof');
cnt = floor(flen/4/3);
Tr = fread(ft, [3,cnt], 'int32'); Tr = Tr';
Tr = Tr + 1;
fclose(ft);

% pairs info
Pr = [];
for j1 = 1:n_-1
    for j2 = j1 + 1 : n_
        m = mod(bPair(floor((j2-1)/8)+1, j1), 2^(mod(j2-1,8)+1));
        m = floor(double(m) / (2^mod(j2-1,8)));
        if m == 0, continue; end
        Pr = cat(1, Pr, [j1,j2]);
    end
end

AL = linspace(0,2*pi,100);
[X,Y,Z] = circlePlane3D(C, N, R, n_, AL);
figure(1); hold off
plot3(X',Y',Z','--');
labels = cell(n_,1);
for i = 1:n_, labels{i} = sprintf('%d', uint32(i)); end
legend(labels);
hold on

if isempty(Tr), return; end

disp(size(Tr,1));
disp(size(Pr,1));
  
figure(3); hold off
for i = 1: size(Pr,1)
    j1 = Pr(i,1); j2 = Pr(i,2);

    js = find((Tr(:,1)==j1 & Tr(:,2)==j2) | ...
              (Tr(:,1)==j1 & Tr(:,3)==j2) | ...
              (Tr(:,2)==j1 & Tr(:,3)==j2) );

    if length(js) < 2, continue; end

    p = [];
    for j = 1: length(js)
        cc1 = C(Tr(js(j),1),:); nn1 = N(Tr(js(j),1),:); 
        cc2 = C(Tr(js(j),2),:); nn2 = N(Tr(js(j),2),:); 
        cc3 = C(Tr(js(j),3),:); nn3 = N(Tr(js(j),3),:); 
        p1 = intersection(cc1, nn1, cc2, nn2, cc3, nn3);

        p = cat(1, p, p1(:)');
    end

    plot3(p(:,1),p(:,2),p(:,3),'-'); hold on
end

