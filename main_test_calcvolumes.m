
% run this m-file to calculate the volumes of intersecting 3D circles

addpath('testutils');

load('matlab.mat');

% calculate cross points of circles
[bPair, ftriname, fptsname] = crossPointsOfCircles3D(C, R, N);

volcnt = volumesOfIntersectingCircles3D(C, R, N, bPair, ftriname, fptsname);

% return;
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

% disp(size(Tr,1));
% disp(size(Pr,1));
  
figure(3); hold off

Trr = Tr;
cntTr = size(Tr,1);
volCnt = 0;
for ii = 1: length(IDX)
    hold off;
    
    while 1
        if isempty(Tr), break; end

        overlapCnt = zeros(size(Tr,1),1);
        i = 1;
        while (i <= size(Pr,1))
            j1 = Pr(i,1); j2 = Pr(i,2);

            js = find((Tr(:,1)==j1 & Tr(:,2)==j2) | ...
                      (Tr(:,1)==j1 & Tr(:,3)==j2) | ...
                      (Tr(:,2)==j1 & Tr(:,3)==j2) );

            if length(js) < 2
                Pr(i,:) = [];
                continue;
            end

            overlapCnt(js) = overlapCnt(js) + 1;
            i = i + 1;
        end

        if isempty(find(overlapCnt<3,1))
            break;
        end

        Tr(overlapCnt<3,:) = [];
    end

    strTr = cell(size(Tr,1),1);
    ptTr = zeros(size(Tr,1),3);
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
            strTr{js(j)} = sprintf('%d,%d,%d',Tr(js(j),1)-1,Tr(js(j),2)-1,Tr(js(j),3)-1);
            ptTr(js(j),:) = p1;
        end

%         plot3(p(:,1),p(:,2),p(:,3),'-'); hold on
    end
    text(ptTr(:,1),ptTr(:,2),ptTr(:,3),strTr);
    
    if cntTr ~= size(Tr,1)
        cntTr = size(Tr,1);
        volCnt = volCnt + 1;
%         fprintf('%d\n', ii-2);
    end
    Tr = Trr;
    Tr(IDX(1:ii),:) = [];
end

disp(' ');
disp(volCnt);

return ;

uV = [0 0 1];
rV = [1 0 0];
figure(3); hold off
for k = 1:n_
    js = find(Tr(:,1)==k | Tr(:,2)==k | Tr(:,3)==k);
    for i = 1: length(js)-1
        Ti = Tr(js(i),:);
        for j = i+1 : length(js)
            b = ~isempty(find(Ti==Tr(js(j),1),1));
            b = b + ~isempty(find(Ti==Tr(js(j),2),1));
            b = b + ~isempty(find(Ti==Tr(js(j),3),1));
            
            if b ~= 2, continue; end
            
            cc1 = C(Tr(js(j),1),:); nn1 = N(Tr(js(j),1),:); 
            cc2 = C(Tr(js(j),2),:); nn2 = N(Tr(js(j),2),:); 
            cc3 = C(Tr(js(j),3),:); nn3 = N(Tr(js(j),3),:); 
            p1 = intersection(cc1, nn1, cc2, nn2, cc3, nn3);
            cc1 = C(Tr(js(i),1),:); nn1 = N(Tr(js(i),1),:); 
            cc2 = C(Tr(js(i),2),:); nn2 = N(Tr(js(i),2),:); 
            cc3 = C(Tr(js(i),3),:); nn3 = N(Tr(js(i),3),:); 
            p2 = intersection(cc1, nn1, cc2, nn2, cc3, nn3);

            p = [p1(:)';p2(:)'];
            plot3(p(:,1),p(:,2),p(:,3),'-o'); hold on
        end
    end
    hold off
end

