%{
Title: QST_v8_Nov9_MR.m
Author: M. Runyon
Description: This script is used to perform quantum state tomography. It
             takes in pixelated intensity data in the form of .txt files
             and astheory_matns a Stokes Vector to every pixel. It also offers
             statistics and various measures of U. It accounts for
             the non-ideal retardance of the QWP and HWP and applies a
             linear inversion technique to more accurately reconstruct the
             measured state. Stokes Parameters are defined as:
                 S1 -> D-A
                 S2 -> R-L
                 S3 -> H-V
             Note that Stokes components S1 and S3 are swapped in final
             calculations due to the absence of a HWP before the QST
             portion of the setup. Additionally, the parameters are defined
             for display purposes in plotting as:
               S1 -> H-V
               S2 -> D-A
               S3 -> R-L
File Dependencies: calcAvgI_v4.m, beamCentreX.m, beamCentreY.m,
                   calcStd.m, stokes2DensityMat.m

**Has no code that actually rotates anything, name only to suggest what
GammaZeta to use with this file.**
%}

%====================================================================%
%                    CONSTANTS                                       %
%====================================================================%
close all;
disp('Declaring constants ...')
theory_vec = [0,0,1];
theory_mat = stokes2DensityMat(theory_vec); %Set this to the conparison
theory_string = '|R>';
dir = strcat('c:\Users\JLundeenLabo\Desktop\MattRunyon\', ... 
    'PolarizationTomography\PublicationData_Nov2\NonUniform2Uniform\PreCorrection\');
if(~exist(dir))
	mkdir(dir);
end
file_name = '0804_110827'; 
rpath = strcat('c:\Users\JLundeenLabo\Desktop\MattRunyon\', ...
    'PolarizationTomography\Measurements\', file_name, '\');
ext = '.txt';       % Declare file type that stores grayscale values
xpix = 1236;        % Declare number of horizontal pixels of CCD camera
ypix = 1626;        % Declare number of vertical pixels of CCD camera
N = 6;				% Declare number of measurements of each projection
writeFigs = 0;
fprintf('Done.\n')

%====================================================================%
%                    DATA PROCESSING                                 %
%====================================================================%

disp(sprintf('Averaging intensity values at each pixel ...'));
[stdIH, IH, BIH, IHall, ~] = calcAvgI_v4(N, 'H', rpath, xpix, ypix, ext);
[stdIV, IV, BIV, IVall, ~] = calcAvgI_v4(N, 'V', rpath, xpix, ypix, ext);	
[stdIL, IL, BIL, ILall, ~] = calcAvgI_v4(N, 'LHC', rpath, xpix, ypix, ext); 
[stdIR, IR, BIR, IRall, ~] = calcAvgI_v4(N, 'RHC', rpath, xpix, ypix, ext);
[stdID, ID, BID, IDall, ~] = calcAvgI_v4(N, 'D', rpath, xpix, ypix, ext);	 
[stdIA, IA, BIA, IAall, ~] = calcAvgI_v4(N, 'A', rpath, xpix, ypix, ext);	
disp(sprintf('Done.\n'));

% Find Gaussian Parameters in each projection
[X01, stdX1] = beamCentreX(IH+IV);
[Y01, stdY1] = beamCentreY(IH+IV);
[X02, stdX2] = beamCentreX(IA+ID);
[Y02, stdY2] = beamCentreY(IA+ID);
[X03, stdX3] = beamCentreX(IR+IL);
[Y03, stdY3] = beamCentreY(IR+IL);

% Take the mean Gaussian parameters
X0vec = [X01,X02,X03];
stdXvec = [stdX1, stdX2, stdX3];
Y0vec = [Y01,Y02,Y03];
stdYvec = [stdY1, stdY2, stdY3];

X0 = round(mean(X0vec));
Y0 = round(mean(Y0vec));
stdX = min(stdXvec);
stdY = min(stdYvec);
theory_matma = min(stdX,stdY);

% Define Beam radius (in pixels) as FWTM (Gaussian)
r_fwtm2 = ceil(sqrt(2*log(10))*theory_mat);

% Construct the Stokes parameters
SS0 = (IH+IV);
SS1 = (ID-IA)./SS0;
SS3 = (IR-IL)./SS0;
SS2 = (IH-IV)./SS0;

% Calculate the average Stokes vector for the whole light field (FWTM)
avg_stokes_vec = [0, 0, 0];
pixels_in_fwtm = 0;
for i = 1:xpix
    for j = 1:ypix
        if (j-Y0)^2 + (i-X0)^2 <= r_fwtm2^2
            if ~isnan(SS1(i,j)) && ~isnan(SS2(i,j)) && ~isnan(SS3(i,j))
                if SS1(i,j) < 10^2 && SS2(i,j) < 10^2 && SS3(i,j) < 10^2
                    avg_stokes_vec(1) = avg_stokes_vec(1) + SS1(i,j);
                    avg_stokes_vec(2) = avg_stokes_vec(2) + SS2(i,j);
                    avg_stokes_vec(3) = avg_stokes_vec(3) + SS3(i,j);
                    pixels_in_fwtm = pixels_in_fwtm + 1;
                end
            end
        end
    end
end
avg_stokes_vec = avg_stokes_vec./pixels_in_fwtm;

% FOR QUADRANT/TETRA BEAM ONLY  -------------------------------------------
%{
Quadrant1 -- lower left
Quadrant2 -- upper left
Quadrant3 -- upper right
Quadrant4 -- lower right
LL : 'Lower Left'
%}
%{
quad1_sum = [0, 0, 0];
quad2_sum = [0, 0, 0];
quad3_sum = [0, 0, 0];
quad4_sum = [0, 0, 0];
quad1_count = 0;
quad2_count = 0;
quad3_count = 0;
quad4_count = 0;
quad_offset = 40;

for i = 1:xpix
    for j = 1:ypix
        if ((j-Y0)^2 + (i-X0)^2 <= r_fwtm2^2) && ... 
                j<(Y0-quad_offset) && i<(X0-quad_offset)
            quad1_sum(1) = quad1_sum(1) + SS1(i,j);
            quad1_sum(2) = quad1_sum(2) + SS2(i,j);
            quad1_sum(3) = quad1_sum(3) + SS3(i,j);
            quad1_count = quad1_count + 1;
        end
        
        if ((j-Y0)^2 + (i-X0)^2 <= r_fwtm2^2) && ... 
                j<(Y0-quad_offset) && i>(X0+quad_offset)
            quad2_sum(1) = quad2_sum(1) + SS1(i,j);
            quad2_sum(2) = quad2_sum(2) + SS2(i,j);
            quad2_sum(3) = quad2_sum(3) + SS3(i,j);
            quad2_count = quad2_count + 1;
        end

        if ((j-Y0)^2 + (i-X0)^2 <= r_fwtm2^2) && ...
                j>(Y0+quad_offset) && i>(X0+quad_offset)
            quad3_sum(1) = quad3_sum(1) + SS1(i,j);
            quad3_sum(2) = quad3_sum(2) + SS2(i,j);
            quad3_sum(3) = quad3_sum(3) + SS3(i,j);
            quad3_count = quad3_count + 1;
        end
        
        if ((j-Y0)^2 + (i-X0)^2 <= r_fwtm2^2) && ...
                j>(Y0+quad_offset) && i<(X0-quad_offset)
            quad4_sum(1) = quad2_sum(1) + SS1(i,j);
            quad4_sum(2) = quad2_sum(2) + SS2(i,j);
            quad4_sum(3) = quad2_sum(3) + SS3(i,j);
            quad4_count = quad2_count + 1;
        end
    end
end
Sq1_theory = 1/sqrt(3).*[-1, 1, -1];
Sq2_theory = 1/sqrt(3).*[-1, -1, 1];
Sq3_theory = 1/sqrt(3).*[1, -1, -1];
Sq4_theory = 1/sqrt(3).*[1, 1, 1];

Sq1 = quad1_sum./quad1_count;
Sq2 = quad2_sum./quad2_count;
Sq3 = quad3_sum./quad3_count;
Sq4 = quad4_sum./quad4_count;

Sq1(3) = -Sq1(3);
Sq2(3) = -Sq2(3);
Sq3(3) = -Sq3(3);
Sq4(3) = -Sq4(3);
%}
% -------------------------------------------------------------------------
%====================================================================%
%                    CALCULATIONS                                    %
%====================================================================%

disp(sprintf('Calculating density matrix ...\n'))
rho = stokes2DensityMat([avg_stokes_vec(1), avg_stokes_vec(2), ...
    avg_stokes_vec(3)]);
disp(sprintf('Done.\n'))

disp(sprintf('Calculating Fidelity, F ...'))
F_dens = calcFidelity(rho,theory_mat);
F = 0.5*(1 + dot(theory_vec,[avg_stokes_vec(1), avg_stokes_vec(2), ... 
    avg_stokes_vec(3)]));
disp(sprintf('Done.\n'))

disp(sprintf('Calculating purity parameter, P ...'))
P = trace(rho^2);
disp(sprintf('Done.\n'))

disp(sprintf('Calculating Uniformity, U ...'))
U = (avg_stokes_vec(1)^2 + avg_stokes_vec(2)^2 + avg_stokes_vec(3)^2)^(0.5)
disp(sprintf('Done.\n'))

% FOR QUADRANT/TETRA BEAM ONLY  -------------------------------------------
%{
fprintf('Calculating density matrices for each quadrant ...\n')
rho1 = stokes2DensityMat(Sq1);
rho2 = stokes2DensityMat(Sq2);
rho3 = stokes2DensityMat(Sq3);
rho4 = stokes2DensityMat(Sq4);
rho1_theory = stokes2DensityMat(Sq1_theory);
rho2_theory = stokes2DensityMat(Sq2_theory);
rho3_theory = stokes2DensityMat(Sq3_theory);
rho4_theory = stokes2DensityMat(Sq4_theory);
disp(sprintf('Done.\n'))

fprintf('Calculating fidelities (density matrix) for each quadrant ...\n')
F1_dens = calcFidelity(rho1_theory,rho1);
F2_dens = calcFidelity(rho2_theory,rho2);
F3_dens = calcFidelity(rho3_theory,rho3);
F4_dens = calcFidelity(rho4_theory,rho4);
disp(sprintf('Done.\n'))

fprintf('Calculating fidelities (dot product) for each quadrant ...\n')
F1 = 0.5*(1+dot(Sq1,Sq1_theory));
F2 = 0.5*(1+dot(Sq2,Sq2_theory));
F3 = 0.5*(1+dot(Sq3,Sq3_theory));
F4 = 0.5*(1+dot(Sq4,Sq4_theory));
disp(sprintf('Done.\n'))

disp(sprintf('Calculating purities in each quadrant ...'))
P1 = trace(rho1^2);
P2 = trace(rho2^2);
P3 = trace(rho3^2);
P4 = trace(rho4^2);
disp(sprintf('Done.\n'))

disp(sprintf('Calculating uniformities in each quadrant ...'))
U1 = (Sq1(1)^2 + Sq1(2)^2 + Sq1(3)^2)^(0.5);
U2 = (Sq2(1)^2 + Sq2(2)^2 + Sq2(3)^2)^(0.5);
U3 = (Sq3(1)^2 + Sq3(2)^2 + Sq3(3)^2)^(0.5);
U4 = (Sq4(1)^2 + Sq4(2)^2 + Sq4(3)^2)^(0.5);
disp(sprintf('Done.\n'))
%}
if writeFigs
    diary(strcat(dir,'QST_Summary.txt'));
end

%====================================================================%
%                    OUTPUT                                          %
%====================================================================%
fprintf(strcat('Data from file ', file_name,'\n'));

fprintf('PRINTING RESULTS WITHIN A FULL WIDTH TENTH MAX BEAM WAIST:\n\n');
fprintf('The average FWTM Stokes vector is [%.3f %.3f %.3f].\n\n', ...
    avg_stokes_vec(1), avg_stokes_vec(2), avg_stokes_vec(3))
fprintf('The uniformity U of the FWTM light field is %.4f\n', U)
fprintf('The density matrix of the FWTM light field is:\n');disp(rho);
fprintf('The quantum purity P of the FWTM light field is: %.3f.\n',P);
fprintf('The fidelity F (trace) of this state with %s', theory_string)
fprintf(' is %0.4f.\n', F_dens);
fprintf('The fidelity F (dot product) of this state with %s', theory_string)
fprintf(' is %0.4f.\n', F);
%}
%{
fprintf('PRINTING RESULTS WITHIN EACH QUADRANT:\n\n');
fprintf('Weighted Scheme ...\n\n');
fprintf('Expected average Stokes vector in Quadrant 1: ')
fprintf('[%.3f %.3f %.3f]\n', Sq1_theory(1), Sq1_theory(2), Sq1_theory(3));
fprintf('Actual average Stokes vector in Quadrant 1: ')
fprintf('[%.3f %.3f %.3f]\n', Sq1(1), Sq1(2), Sq1(3));
fprintf('Fidelity (density matrix) of Quadrant 1: %.3f\n', F1_dens);
fprintf('Fidelity (dot product) of Quadrant 1: %.3f\n', F1);
fprintf('Purity of Quadrant 1: %.3f\n', P1);
fprintf('Uniformity of Quadrant 1: %.3f\n\n',U1);

fprintf('Expected average Stokes vector in Quadrant 2: ')
fprintf('[%.3f %.3f %.3f]\n', Sq2_theory(1), Sq2_theory(2), Sq2_theory(3));
fprintf('Actual average Stokes vector in Quadrant 2: ')
fprintf('[%.3f %.3f %.3f]\n', Sq2(1), Sq2(2), Sq2(3));
fprintf('Fidelity (density matrix) of Quadrant 2: %.3f\n', F2_dens);
fprintf('Fidelity (dot product) of Quadrant 2: %.3f\n', F2);
fprintf('Purity of Quadrant 2: %.3f\n', P2);
fprintf('Uniformity of Quadrant 2: %.3f\n\n',U2);

fprintf('Expected average Stokes vector in Quadrant 3: ')
fprintf('[%.3f %.3f %.3f]\n', Sq3_theory(1), Sq3_theory(2), Sq3_theory(3));
fprintf('Actual average Stokes vector in Quadrant 3: ')
fprintf('[%.3f %.3f %.3f]\n', Sq3(1), Sq3(2), Sq3(3));
fprintf('Fidelity (density matrix) of Quadrant 3: %.3f\n', F3_dens);
fprintf('Fidelity (dot product) of Quadrant 3: %.3f\n', F3);
fprintf('Purity of Quadrant 3: %.3f\n', P3);
fprintf('Uniformity of Quadrant 3: %.3f\n\n',U3);

fprintf('Expected average Stokes vector in Quadrant 4: ')
fprintf('[%.3f %.3f %.3f]\n', Sq4_theory(1), Sq4_theory(2), Sq4_theory(3));
fprintf('Actual average Stokes vector in Quadrant 4: ')
fprintf('[%.3f %.3f %.3f]\n', Sq4(1), Sq4(2), Sq4(3));
fprintf('Fidelity (density matrix) of Quadrant 4: %.3f\n', F4_dens);
fprintf('Fidelity (dot product) of Quadrant 4: %.3f\n', F4);
fprintf('Purity of Quadrant 4: %.3f\n', P4);
fprintf('Uniformity of Quadrant 4: %.3f\n\n',U4);
%}
if writeFigs
    diary off;
end
%====================================================================%
%                    PLOTTING                                        %
%====================================================================%

angle_deg = -1.7;
disp(sprintf('Plotting ...'))
t = linspace(0,2*pi);
figure(1);
h7=pcolor((IH+IV));
colormap(jet);
colorbar
set(h7,'LineStyle','none')
set(gcf,'color','w');
axis equal
axis tight
hold on;
plot(r_fwtm2.*cos(t)+Y0,r_fwtm2.*sin(t)+X0, 'k-.', 'LineWidth', 2);
hold off;
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
if writeFigs
    saveas(h7, strcat(dir,'Intensity_FWTM'),'epsc');
    saveas(h7, strcat(dir,'Intensity_FWTM'),'svg');
    saveas(h7, strcat(dir,'Intensity_FWTM'),'pdf');
    saveas(h7, strcat(dir,'Intensity_FWTM.jpg'));
end

figure(2);
h8=pcolor(((imrotate(SS1,angle_deg,'crop'))));
colormap(jet);
colorbar
caxis([-1,1]);
set(h8,'LineStyle','none')
set(gcf,'color','w');
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
hold on;
plot(r_fwtm2.*cos(t)+Y0,r_fwtm2.*sin(t)+X0, 'k-.', 'LineWidth', 2);
hold off;
if writeFigs
    saveas(h8, strcat(dir,'S1'),'epsc');
    saveas(h8, strcat(dir,'S1'),'svg');
    saveas(h8, strcat(dir,'S1'),'pdf');
    saveas(h8, strcat(dir,'S1.jpg'));
end

figure(3)
h9=pcolor(((imrotate(SS2,angle_deg,'crop'))));
colormap(jet);
colorbar
caxis([-1,1]);
set(h9,'LineStyle','none')
set(gcf,'color','w');
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
hold on;
plot(r_fwtm2.*cos(t)+Y0,r_fwtm2.*sin(t)+X0, 'k-.', 'LineWidth', 2);
hold off;
if writeFigs
    saveas(h9, strcat(dir,'S2'),'epsc');
    saveas(h9, strcat(dir,'S2'),'svg');
    saveas(h9, strcat(dir,'S2'),'pdf');
    saveas(h9, strcat(dir,'S2.jpg'));
end

figure(4)
h10=pcolor(((imrotate(SS3,angle_deg,'crop'))));
colormap(jet);
colorbar
caxis([-1,1]);
set(h10,'LineStyle','none')
set(gcf,'color','w');
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
%title('S3 Gradient - FWTM waist, unweighted');
hold on;
plot(r_fwtm2.*cos(t)+Y0,r_fwtm2.*sin(t)+X0, 'k-.', 'LineWidth', 2);
hold off;
if writeFigs
    saveas(h10, strcat(dir,'S3'),'epsc');
    saveas(h10, strcat(dir,'S3'),'svg');
    saveas(h10, strcat(dir,'S3'),'pdf');
    saveas(h10, strcat(dir,'S3.jpg'));
end
disp(sprintf('Done.'));