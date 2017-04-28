%{
Title: calibratePhase.m
Author: M. Runyon
Date: February 21st
%}
close all
clear phaseMatrix1;
clear grayMatrix1;
clear maps1;
clear mu1;
clear fits1;
clear theta3Dmat1;
clear vecs;
clear theta3d;
clear grayVec;
clear x;
clear y;
clear yy;

% Declare data structures
numRows = length(S4d(:,1,1,1));
numCols = length(S4d(1,:,1,1));
numComp = length(S4d(1,1,:,1));
numRuns = length(S4d(1,1,1,:));
fitSteps = 200;
prec = 1e2;
startingGray = 53;
endingGray = 253;
grayVec = linspace(startingGray,endingGray, fitSteps+1);
x = linspace(startingGray,endingGray,numRuns); % Experimental x
axes = zeros(numRows,numCols,3);
theta3Dmat1 = zeros(numRows,numCols,numRuns);
total3Dphase1 = zeros(numRows,numCols,numRuns);
fits1 = zeros(numRows,numCols,8);
deltas = zeros(numRows,numCols);
mu1 = zeros(numRows,numCols,2);
maps1 = cell(numRows,numCols);
phaseMatrix1 = zeros(numRows,numCols,fitSteps+1);
step = 5;
writeFigs = 0;
wpathP = 'c:\Users\JLundeenLabo\Desktop\MattRunyon\PolarizationTomography\Calibration\Apr19\LHS\PhaseGrayFiles\';

% Determine phase-gray relatonship at every pixel (i,j) using k data sets
tic
for i = 1:numRows % For each row of CCD
    fprintf('Determining phase-gray for pixels in row %d/1236 ...\n', i);
    for j = 1:numCols % For each column in said row
        vecs = zeros(numRuns,3);
        for q = 1:numRuns
            vecs(q,:)=[S4d(i,j,1,q), S4d(i,j,2,q), S4d(i,j,3,q)];
        end
        % Find the axis/plane of rotation for pixel (i,j)
        try
            [n,V,p] = affine_fit(vecs);
        catch
            continue;
        end
        axes1(i,j,:) = n;
        % Find rotation angles such that n = S2
        phi = acos(dot([1,0],[n(2), n(3)]));
        theta = acos(dot([0,1],[n(1), n(2)]));

        % Rotate coordinate system accordingly
        R_s3 = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0;...
                0, 0, 1];
        R_s1 = [1, 0, 0; 0, cos(phi), -sin(phi); ...
                0, sin(phi), cos(phi)];

        % Project each vector onto the plane of rotation
        % Rotate coordinate system to be in S1-S3 plane
        vecsProj = zeros(size(vecs));
        for k = 1:length(vecsProj)
           vecsProj(k,:) = vecs(k,:)-dot(vecs(k,:),n').*n';
           vecsProjRot(k,:) = R_s3*R_s1*vecsProj(k,:)';
        end
        
        % For each run in calibration, determine phase delta
        for m = 1:length(vecsProj)-1
        
            % Get initial and final Stokes vectors
            p0 = vecsProjRot(1,:);
            p1 = vecsProjRot(m,:);
            p2 = vecsProjRot(m+1,:);

            % Calculate 3D phase traversed
            u = cross(p1,p2);

            theta3d = atan2(u(2),dot(p1,p2));
            theta3Dmat1(i,j,m+1) = theta3d+theta3Dmat1(i,j,m);
        end
        
        % Fit the phase-gray data with 7th order polynomial
        y = squeeze(theta3Dmat1(i,j,:))';
        [P, S, MU] = polyfit(x,y,7);
        fits1(i,j,:) = P;
        mu1(i,j,:) = MU;
        [yy, ~] = polyval(squeeze(fits1(i,j,:))',grayVec,S,MU);
        phaseMatrix1(i,j,:) = round(yy*prec)/prec;     
    end 
end

wpath = 'd:\Matlab\April27\LHS_mat\';
for i = 1:size(phaseMatrix1,1)
    maps = squeeze(phaseMatrix1(i,:,:));
    fileName = strcat(wpath, 'maps', int2str(i));
    save(fileName, 'maps')
end

% Plot to check the fit
figure(1)
r = 1;
[xp,yp,zp] = sphere(50);
lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
surf(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
grid off;
set(gca,'color','none');
set(gca,'visible','off');
plot3(squeeze(S4d(725,725,1,:))',squeeze(S4d(725,725,2,:))',squeeze(S4d(725,725,3,:))','b.');
title('Plane and Axis of Rotation -- (758,723)')
plot3([0 n(1)], [0 n(2)], [0 n(3)],'b','linewidth',2)
plot3([0 -1],[0 0],[0 0], 'k' ) % for x-axis
plot3([0 0],[0 -1],[0 0], 'k') % for y-axis
plot3([0 0],[0 0],[0 1], 'k') % for z-axis
text(-1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
text(0,-1,0,'S2','HorizontalAlignment','left','FontSize',10);
text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
[X,Y] = meshgrid(linspace(-1,1,5));
surf(X,Y, - (n(1)/n(3)*X+n(2)/n(3)*Y-dot(n,p)/n(3)),'facecolor','blue','facealpha',0.3);
axis([-1 1 -1 1 -1 1])
%view([0 0])
%{
if writeFigs
    print(strcat(wpathP,'PlaneFit_',point),'-djpeg');
end
%}
hold off;

% Plot to check the rotation
[n,V,p] = affine_fit(vecsProjRot);

figure(2)
r = 1;
[xp,yp,zp] = sphere(50);
lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
surf(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
grid off;
set(gca,'color','none');
set(gca,'visible','off');
plot3(vecsProjRot(:,1),vecsProjRot(:,2),vecsProjRot(:,3),'b.');
title('Plane and Axis of Rotation -- (758,723)')
plot3([0 n(1)], [0 n(2)], [0 n(3)],'b','linewidth',2)
plot3([0 1],[0 0],[0 0], 'k' ) % for x-axis
plot3([0 0],[0 1],[0 0], 'k') % for y-axis
plot3([0 0],[0 0],[0 1], 'k') % for z-axis
text(1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
text(0,1,0,'S2','HorizontalAlignment','left','FontSize',10);
text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
[X,Y] = meshgrid(linspace(-1,1,5));
surf(X,Y, - (n(1)/n(3)*X+n(2)/n(3)*Y-dot(n,p)/n(3)),'facecolor','blue','facealpha',0.3);
axis([-1 1 -1 1 -1 1])
view([-37.5+180 30])
%{
if writeFigs
    print(strcat(wpathP,'PlaneFit_',point),'-djpeg');
end
%}
hold off;

figure(3);
plot(grayVec,polyval(squeeze(fits1(725,725,:))',grayVec,[],squeeze(mu1(725,725,:))'),'b')
title('Fitted Data vs. Experimental Data -- (725,725)');
xlabel('Grayscale');
ylabel('Phase (rad)');
hold on;
plot(x, squeeze(theta3Dmat1(725,725,:))','r');
hold off;
legend('Polynomial Fit', 'Experiment');

figure(4);
plot(grayVec,polyval(squeeze(fits1(900,1000,:))',grayVec,[],squeeze(mu1(900,1000,:))'),'b')
title('Fitted Data vs. Experimental Data -- (900,1000)');
xlabel('Grayscale');
ylabel('Phase (rad)');
hold on;
plot(x, squeeze(theta3Dmat1(900,1000,:))','r');
hold off;
legend('Polynomial Fit', 'Experiment');

figure(5);
plot(grayVec,polyval(squeeze(fits1(600,600,:))',grayVec,[],squeeze(mu1(600,600,:))'),'b')
title('Fitted Data vs. Experimental Data -- (600,600)');
xlabel('Grayscale');
ylabel('Phase (rad)');
hold on;
plot(x, squeeze(theta3Dmat1(600,600,:))','r');
hold off;
legend('Polynomial Fit', 'Experiment');

figure(6);
plot(grayVec,polyval(squeeze(fits1(400,800,:))',grayVec,[],squeeze(mu1(400,800,:))'),'b')
title('Fitted Data vs. Experimental Data -- (400,800)');
xlabel('Grayscale');
ylabel('Phase (rad)');
hold on;
plot(x, squeeze(theta3Dmat1(400,800,:))','r');
hold off;
legend('Polynomial Fit', 'Experiment');

figure(7);
plot(grayVec,polyval(squeeze(fits1(600,700,:))',grayVec,[],squeeze(mu1(600,700,:))'),'b')
title('Fitted Data vs. Experimental Data -- (600,700)');
xlabel('Grayscale');
ylabel('Phase (rad)');
hold on;
plot(x, squeeze(theta3Dmat1(600,700,:))','r');
hold off;
legend('Polynomial Fit', 'Experiment');


figure(8)
r = 1;
[xp,yp,zp] = sphere(50);
lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
surf(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
grid off;
set(gca,'color','none');
set(gca,'visible','off');
plot3(squeeze(S4d(670, 1140,1,:))',squeeze(S4d(670, 1140,2,:))',squeeze(S4d(670, 1140,3,:))','b.');
title('Plane and Axis of Rotation -- (670, 1140)')
plot3([0 n_pix1(1)], [0 n_pix1(2)], [0 n_pix1(3)],'b','linewidth',2)
plot3([0 1],[0 0],[0 0], 'k' ) % for x-axis
plot3([0 0],[0 1],[0 0], 'k') % for y-axis
plot3([0 0],[0 0],[0 1], 'k') % for z-axis
text(1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
text(0,1,0,'S2','HorizontalAlignment','left','FontSize',10);
text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
[X,Y] = meshgrid(linspace(-1,1,5));
surf(X,Y, - ((n_pix1(1)/n_pix1(3)*X)+(n_pix1(2)/n_pix1(3))*Y - (dot(n_pix1,p_pix1)/n_pix1(3))),'facecolor','blue','facealpha',0.3);
axis([-1 1 -1 1 -1 1])
view([-37.5+180 30])

figure(9)
r = 1;
[xp,yp,zp] = sphere(50);
lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
surf(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
grid off;
set(gca,'color','none');
set(gca,'visible','off');
plot3(squeeze(S4d(480, 840,1,:))',squeeze(S4d(480, 840,2,:))',squeeze(S4d(480, 840,3,:))','b.');
title('Plane and Axis of Rotation -- (480, 840)')
plot3([0 n_pix2(1)], [0 n_pix2(2)], [0 n_pix2(3)],'b','linewidth',2)
plot3([0 1],[0 0],[0 0], 'k' ) % for x-axis
plot3([0 0],[0 1],[0 0], 'k') % for y-axis
plot3([0 0],[0 0],[0 1], 'k') % for z-axis
text(1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
text(0,1,0,'S2','HorizontalAlignment','left','FontSize',10);
text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
[X,Y] = meshgrid(linspace(-1,1,5));
surf(X,Y, - (n_pix2(1)/n_pix2(3)*X+n_pix2(2)/n_pix2(3)*Y-dot(n_pix2,p_pix2)/n_pix2(3)),'facecolor','blue','facealpha',0.3);
axis([-1 1 -1 1 -1 1])
view([-37.5+180 30])

figure(10)
r = 1;
[xp,yp,zp] = sphere(50);
lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
surf(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
grid off;
set(gca,'color','none');
set(gca,'visible','off');
plot3(squeeze(S4d(630, 570,1,:))',squeeze(S4d(630, 570,2,:))',squeeze(S4d(630, 570,3,:))','b.');
title('Plane and Axis of Rotation -- (630, 570)')
plot3([0 n_pix3(1)], [0 n_pix3(2)], [0 n_pix3(3)],'b','linewidth',2)
plot3([0 1],[0 0],[0 0], 'k' ) % for x-axis
plot3([0 0],[0 1],[0 0], 'k') % for y-axis
plot3([0 0],[0 0],[0 1], 'k') % for z-axis
text(1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
text(0,1,0,'S2','HorizontalAlignment','left','FontSize',10);
text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
[X,Y] = meshgrid(linspace(-1,1,5));
surf(X,Y, - (n_pix3(1)/n_pix3(3)*X+n_pix3(2)/n_pix3(3)*Y-dot(n_pix3,p_pix3)/n_pix3(3)),'facecolor','blue','facealpha',0.3);
axis([-1 1 -1 1 -1 1])
view([-37.5+180 30])

figure(11)
r = 1;
[xp,yp,zp] = sphere(50);
lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
surf(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
grid off;
set(gca,'color','none');
set(gca,'visible','off');
plot3(squeeze(S4d(840, 580,1,:))',squeeze(S4d(840, 580,2,:))',squeeze(S4d(840, 580,3,:))','b.');
title('Plane and Axis of Rotation -- (840, 580)')
plot3([0 n_pix4(1)], [0 n_pix4(2)], [0 n_pix4(3)],'b','linewidth',2)
plot3([0 1],[0 0],[0 0], 'k' ) % for x-axis
plot3([0 0],[0 1],[0 0], 'k') % for y-axis
plot3([0 0],[0 0],[0 1], 'k') % for z-axis
text(1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
text(0,1,0,'S2','HorizontalAlignment','left','FontSize',10);
text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
[X,Y] = meshgrid(linspace(-1,1,5));
surf(X,Y, - (n_pix4(1)/n_pix4(3)*X+n_pix4(2)/n_pix4(3)*Y-dot(n_pix4,p_pix4)/n_pix4(3)),'facecolor','blue','facealpha',0.3);
axis([-1 1 -1 1 -1 1])
view([-37.5+180 30])
