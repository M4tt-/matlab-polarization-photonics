%{
Title: StokesPolarimetry.m
Author: M. Runyon
Description: This script is used to perform Stokes Polarimetry to
             characterize the polarization state of a classical light
             field. If the field is completely depolarized, a quantum 
             purity of 0.5 is assigned to the light field. A completely
             polarized light field has a quantum purity of 1.

             Note: Quantum purity, in this context, is the same as the 
                   degree of polarization.

             The following convention is used for Stokes vectors [S1 S2 S3]
             
                 |H> (horizontal) = [1 0 0]
                 |V> (vertical) = [-1 0 0]
                 |D> (diagonal) = [0 1 0]
                 |A> (antidiagonal) = [0 -1 0]
                 |R> (right-hand circular) = [0 0 1]
                 |L> (left-hand circular) = [0 0 -1]

            Also,

                |D> = 1/sqrt(2)*(|H> + |V>)
                |R> = 1/sqrt(2)*(|H> + i|V>)

If it is desired to compensate for non-ideal waveplate retardance, see 
pages 6-8 of this paper:
http://research.physics.illinois.edu/QI/Photonics/
       tomography-files/amo_tomo_chapter.pdf
Note that the Stokes vector convention is different here.

Input: The only input required by this script is intensity data for each
       projection of the light field of interest onto each axis of the
       Poincare Sphere. This data should be in the form of a tab delimited
       .txt file, where each entry corresponds to the intensity measured
       at a pixel of detector (CCD camera). Each text file should be named
       so that the first letter corresponds to the projection its for
       followed by an underscore and number, where the number is ordered
       from 1 to N (see defn of N below). For example, if one wishes to
       save 3 files for intensity data corresponding to a horizontal
       projection, they would be named H_1.txt, H_2.txt, H_3.txt, etc. It
       is also assumed that N data captures of only the ambient light has
       been made. These should be named 'BACK_1.txt, BACK_2.txt' etc.

       xpix, ypix = dimensions of the CCD camera's pixelated grid
       N = # of CCD captures performed per projection
       rpath = directory storing the .txt files

        
%}

disp('Declaring constants ...')
rpath = 'c:\SCHOOL\Graduate\LundeenLab\PolTom\Measurements\0217_171045\';
xpix = 1236;        % Declare number of horizontal pixels of CCD camera
ypix = 1626;        % Declare number of vertical pixels of CCD camera
N = 6;				% Declare number of measurements of each projection for averaging
I = [1 0;0 1];        % Identity matrix
sigz = [0 1;1 0];     % Pauli Matrix 1
sigy = [0 -1i; 1i 0]; % Pauli Matrix 2
sigx = [1 0;0 -1];    % Pauli Matrix 3
fprintf('Done.\n\n');

%--------------------------------------------------------------------%
%               Get Intensity data                                   %
%--------------------------------------------------------------------%
fprintf('Reading intensity data ...');
for i = 1 : N
    
   filenameH=strcat(rpath,'H_',int2str(i),'.txt');
   IH(:,:,i)=importdata(filenameH);
   
   filenameV=strcat(rpath,'V_',int2str(i),'.txt');
   IV(:,:,i)=importdata(filenameV);
   
   filenameD=strcat(rpath,'D_',int2str(i),'.txt');
   ID(:,:,i)=importdata(filenameD);
  
   filenameA=strcat(rpath,'A_',int2str(i),'.txt');
   IA(:,:,i)=importdata(filenameA);
   
   filenameR=strcat(rpath,'R_',int2str(i),'.txt');
   IR(:,:,i)=importdata(filenameR);
   
   filenameL=strcat(rpath,'L_',int2str(i),'.txt');
   IL(:,:,i)=importdata(filenameL);
   
   filenameB=strcat(rpath, 'BACK_',int2str(i),'.txt');
   IB(:,:,i)=importdata(filenameB);
   
end
fprintf('Done.\n\n');

%--------------------------------------------------------------------%
%               Average intensity data                               %
%--------------------------------------------------------------------%
fprintf('Averaging intensity data ...');
IHavg=mean(IH,3);
IVavg=mean(IV,3);
IDavg=mean(ID,3);
IAavg=mean(IA,3);
IRavg=mean(IR,3);
ILavg=mean(IL,3);
IBavg=mean(IB,3);
fprintf('Done.\n\n');

%--------------------------------------------------------------------%
%               Construct reduced Stokes vectors                     %
%--------------------------------------------------------------------%
fprintf('Constructing Stokes vectors for each pixel ...');
S1top=IHavg-IVavg;
S1bottom=IHavg+IVavg-2*IBavg;
S1=S1top./S1bottom;

S2top=IDavg-IAavg;
S2bottom=IHavg+IVavg-2*IBavg;
S2=S2top./S2bottom;

S3top=IRavg-ILavg;
S3bottom=IHavg+IVavg-2*IBavg;
S3=S3top./S3bottom;

S0=sqrt(S1.^2+S2.^2+S3.^2);
fprintf('Done.\n\n');

fprintf('Constructing Stokes vectors for entire light field ...');
s0 = sum(sum(IHavg)) + sum(sum(IVavg));
s1 = (sum(sum(IHavg)) - sum(sum(IVavg)))/s0;
s2 = (sum(sum(IDavg)) - sum(sum(IAavg)))/s0;
s3 = (sum(sum(IRavg)) - sum(sum(ILavg)))/s0;
fprintf('Done.\n\n')

fprintf('Calculating density matrix ...\n\n')
rho = 1/2.*(I + s3*sigz+s2*sigy+s1*sigx);
fprintf('Done.\n')

fprintf('Calculating purity parameter ...\n')
P = trace(rho^2);
fprintf('Done.\n\n')

%--------------------------------------------------------------------%
%               OUTPUT                                               %
%--------------------------------------------------------------------%

fprintf('The Stokes vector of thelight field is given by: [%.4f %.4f %.4f].\n', s1,s2,s3);
fprintf('The density matrix (rho) of the light field is given by:\n\n');disp(rho);
fprintf('The trace of rho is given by: %.2f\n', trace(rho));
fprintf('The purity parameter of the light field is %.4f.\n\n', P);
fprintf('Plotting ...\n')

figure(1);
h1=pcolor(IHavg+IVavg);
colormap(jet);
colorbar
set(h1,'LineStyle','none')
axis equal
axis tight
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
title('Intensity Profile');

figure(2)
h3 = surf(IHavg+IVavg);
set(h3,'LineStyle','none')
title('Intensity Surface Plot')
xlabel('Pixels(X)')
ylabel('Pixels(Y)')


figure(3)
h4=pcolor(S1);
colormap(jet);
colorbar
caxis([-1,1]);
set(h4,'LineStyle','none')
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
title('S1 Gradient');

figure(4)
h5=pcolor(S2);
colormap(jet);
colorbar
caxis([-1,1]);
set(h5,'LineStyle','none')
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
title('S2 Gradient');

figure(5)
h6=pcolor(S3);
colormap(jet);
colorbar
caxis([-1,1]);
set(h6,'LineStyle','none')
xlabel('Pixels (X)');
ylabel('Pixels (Y)');
title('S3 Gradient');
fprintf('Done.\n\n');