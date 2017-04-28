%{
Title: QST_v6_getS4d.m
Author: M. Runyon
Description: This script is used or calibration purposes. It has the same
             function as QST_v6.m, but this functionality is placed in a
             loop that iterates over many measurements for calibration. It
             is the basis in providing a phase-gray relationship for each
             pixel of the SLM. It also offers statistics on any 'bad'
             polarizations at any pixels. Assumes the retardances of wave
             plates have been measured.

File Dependencies: calcAvgI_v4.m,
                   stokes2DensityMat.m
%}
%{
% LHS Calibration Files (April 26th, 2017)
trials = cellstr(['0426_133035';'0426_133407';'0426_133739'; ...
    '0426_134110';'0426_134441';'0426_134813';'0426_135144';'0426_135516'; ...
    '0426_135847';'0426_140220';'0426_140551';'0426_140923';'0426_141255'; ...
    '0426_141628';'0426_141959';'0426_142331';'0426_142702';'0426_143034'; ...
    '0426_143405';'0426_143737';'0426_144110';'0426_144444';'0426_144816'; ...
    '0426_145147';'0426_145658';'0426_150025';'0426_150357';'0426_150728'; ...
    '0426_151100';'0426_151431';'0426_151803';'0426_152135';'0426_152507'; ...
    '0426_152839';'0426_153211';'0426_153545';'0426_153921';'0426_154252'; ...
    '0426_154629';'0426_155008';'0426_155339';]);
%}
% RHS Calibration Files (April 27th, 2017)
trials = cellstr(['0427_102329';'0427_102705';'0427_103046'; ...
    '0427_103436';'0427_103813';'0427_104202';'0427_104544';'0427_104917'; ...
    '0427_105251';'0427_105643';'0427_110014';'0427_110359';'0427_110732'; ...
    '0427_111110';'0427_111446';'0427_111818';'0427_112152';'0427_112528'; ...
    '0427_112901';'0427_113238';'0427_113613';'0427_113945';'0427_114320'; ...
    '0427_115022';'0427_120016';'0427_120352';'0427_120723';'0427_121055'; ...
    '0427_121446';'0427_121820';'0427_122156';'0427_122529';'0427_122904'; ...
    '0427_123236';'0427_123608';'0427_123941';'0427_124316';'0427_124652'; ...
    '0427_125024';'0427_125405';'0427_125737']);

% Relevant data structures to populate
disp('Declaring constants ...')
S4d = zeros(1236,1626,3);
numBadPix = zeros(1,length(trials));
numNansS0 = zeros(1,length(trials));

% Constants
%rpath = 'c:\Users\JLundeenLabo\Desktop\MattRunyon\PolarizationTomography\Measurements\1123_174219\'; % Path for the folder housing tomography images
wpath = 'c:\Users\Jlundeenlabo\PolTom\Meaurements\';
ext = '.txt';       % Declare file type that stores grayscale values
xpix = 1236;        % Declare number of horizontal pixels of CCD camera
ypix = 1626;        % Declare number of vertical pixels of CCD camera
N = 6;				% Declare number of measurements of each projection for averaging
writeFigs =0;
S_field = zeros(length(trials),3);
S_pix1 = zeros(length(trials),3);
S_pix2 = zeros(length(trials),3);
S_pix3 = zeros(length(trials),3);
S_pix4 = zeros(length(trials),3);
fprintf('Done.\n')

for k = 1:length(trials)
    diary off;
    close all;

    file = (char(trials(k)));
    rpath = strcat('c:\Users\Jlundeenlabo\PolTom\Measurements\',file,'\');
    
    fprintf('Averaging intensity values at each pixel ...');
    [stdIH, IH, BIH, IHall, ~] = calcAvgI_v4(N, 'H', rpath, xpix, ypix, ext);	% Average Horizontal Projection - xpix x ypix double matrix
    [stdIV, IV, BIV, IVall, ~] = calcAvgI_v4(N, 'V', rpath, xpix, ypix, ext);	% Average Vertical Projection - xpix x ypix double matrix
    [stdIL, IL, BIL, ILall, ~] = calcAvgI_v4(N, 'LHC', rpath, xpix, ypix, ext);   % Average Left Hand Circular Projection - xpix x ypix double matrix
    [stdIR, IR, BIR, IRall, ~] = calcAvgI_v4(N, 'RHC', rpath, xpix, ypix, ext);% Average Right Hand Circular Projection - xpix x ypix double matrix
    [stdID, ID, BID, IDall, ~] = calcAvgI_v4(N, 'D', rpath, xpix, ypix, ext);	% Average Diagonal Projection - xpix x ypix double matrix
    [stdIA, IA, BIA, IAall, ~] = calcAvgI_v4(N, 'A', rpath, xpix, ypix, ext);	% Average Antidiagonal Projection - xpix x ypix double matrix
    fprintf('Done.\n');

    fprintf('Total Intensity in |H>: %d\n', sum(sum(IH)));
    fprintf('Total Intensity in |V>: %d\n', sum(sum(IV)));
    fprintf('Total Intensity in |D>: %d\n', sum(sum(ID)));
    fprintf('Total Intensity in |A>: %d\n', sum(sum(IA)));
    fprintf('Total Intensity in |R>: %d\n', sum(sum(IR)));
    fprintf('Total Intensity in |L>: %d\n', sum(sum(IL)));
    fprintf('Total Background Intensity: %d\n', sum(sum(BIH)));

    %{
    fprintf('Finding Gaussian Parameters ...')
    % Find Gaussian Parameters in each projection
    [X01, stdX1] = beamCentreX(IH+IV);
    [Y01, stdY1] = beamCentreY(IH+IV);
    [X02, stdX2] = beamCentreX(IA+ID);
    [Y02, stdY2] = beamCentreY(IA+ID);
    [X03, stdX3] = beamCentreX(IR+IL);
    [Y03, stdY3] = beamCentreY(IR+IL);

    % Take the mean Gaussian parameters
    sigma = min(min([stdX1, stdX2, stdX3]),min([stdY1, stdY2, stdY3]));

    % 2D Gaussian fit
    [X,Y] = meshgrid(1:ypix,1:xpix);
    g2dm = double(max(max(IR+IL))).*exp(-((X-double(round(mean([X01,X02,X03]))).^2)./(2.*double(min([stdX1, stdX2, stdX3])).^2) ...
            -((Y-double(round(mean([Y01,Y02,Y03])))).^2)./(2.*double(min([stdY1, stdY2, stdY3])).^2)));
    g2d = @ (x,y) double(max(max(IR+IL))).*exp(-((x-double(round(mean([X01,X02,X03])))).^2)./(2.*double(min([stdX1, stdX2, stdX3])).^2) ...
            -((y-double(round(mean([Y01,Y02,Y03])))).^2)./(2.*double(min([stdY1, stdY2, stdY3])).^2));
    fprintf('Done.\n')
    %}
    fprintf('Begin waveplate compensation:\n');
    fprintf('Constructing inversion matrix ...')
    Tmat0 = [1.00000000000000022204, -0.035235598019462792407, 0.0019364022625764736301, 0.00046617334881802119026];
    Tmat1 = [7.2807584862345766142E-17, -0.99948794683696196017, -0.0029837764449056765857, 0.049895544452960322035];
    Tmat2 = [-2.1718837704922462108E-16, -0.00079293476769616024445, 1.00092015145929269515, -0.022493306973316081615];
    Tmat3 = [-2.2444507888817004976E-16, 0.035235598019462792407, -0.0019364022625764736301, 0.99953382665118195582];
    T = [Tmat0;Tmat1;Tmat2;Tmat3];
    fprintf('Done.\n')

    fprintf('Constructing T-parameters for each pixel ...')
    T0_p = (IH+IV);
    T1_p = (ID-IA);
    T2_p = (IR-IL);
    T3_p = (IH-IV);
    T_p = cat(3,T0_p, T1_p, T2_p, T3_p);
    S_temp = (zeros(size(T_p)));
    S = (zeros(size(T_p)));
    for i = 1:xpix
        for j=1:ypix
            t = [T_p(i,j,1);T_p(i,j,2);T_p(i,j,3);T_p(i,j,4)];
            S_temp(i,j,:) = T*t;
            S(i,j,:) = S_temp(i,j,:)./S_temp(i,j,1);
        end
    end
    SS0 = S(:,:,1);SS1 = S(:,:,2);SS2=S(:,:,4);SS3=S(:,:,3);
    fprintf('Done.\n')

    fprintf('Constructing T-parameters for the entire light field ...')
    T0 = sum(sum(IH))+sum(sum(IV));
    T1 = sum(sum(ID))-sum(sum(IA));
    T2 = sum(sum(IR))-sum(sum(IL));
    T3 = sum(sum(IH))-sum(sum(IV));
    fprintf('Done.\n')

    fprintf('Converting T-Parameters to Stokes vectors ...')
    s_temp = T*[T0;T1;T2;T3];
    sp = s_temp./s_temp(1);
    s = [sp(1), sp(4), sp(3), sp(2)]; % Stokes vec for light field.
    fprintf('Done.\n')

    fprintf('Calculating density matrix ...\n')
    rho = stokes2DensityMat([s(2), s(3), s(4)]);
    fprintf('Done.\n')

    fprintf('Calculating purity parameter ...')
    P = trace(rho^2);
    fprintf('Done.\n')

    fprintf('The Stokes vector of the light field is [%.6f %.6f %.6f %.6f]\n', s(4), s(2), s(3));
    fprintf('The density matrix (rho) is given by:\n');
    disp(rho);
    fprintf('The trace of the density matrix is given by: %.2f\n', trace(rho));
    fprintf('The purity parameter of the light field (average intensities) is %.4f.\n\n', P);

    S4d(:,:,1,k) = SS1;
    S4d(:,:,2,k) = SS2;
    S4d(:,:,3,k) = SS3;
    S_field(k,:) = [s(2), s(3), s(4)];
    S_pix1(k,:) = [SS1(670, 1140),SS2(670, 1140),SS3(670, 1140)];
    S_pix2(k,:) = [SS1(480, 840),SS2(480, 840),SS3(480, 840)];
    S_pix3(k,:) = [SS1(630, 570),SS2(630, 570),SS3(630, 570)];
    S_pix4(k,:) = [SS1(840, 580),SS2(840, 580),SS3(840, 580)];
    
end
[n,V,p] = affine_fit(S_field);
[n_pix1,V_pix1,p_pix1] = affine_fit(S_pix1);
[n_pix2,V_pix2,p_pix2] = affine_fit(S_pix2);
[n_pix3,V_pix3,p_pix3] = affine_fit(S_pix3);
[n_pix4,V_pix4,p_pix4] = affine_fit(S_pix4);
theta = acos(dot(n,[0,1,0])/(norm([0,1,0])*norm(n)));
fprintf('The angle between the axis of rotation and S2 is %d\n', theta);