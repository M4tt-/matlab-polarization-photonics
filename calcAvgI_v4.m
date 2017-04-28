%{
Title: calcAvgI.m
Author: M. Runyon
Description: This is a function file that will average matrix entries N times.

@param N: Integer, number of files to be averaged
@param e: Character, specifies which basis the intensity profile was measured in
@param path: String, the full path to the folder housing the intensity profile data
@param m: Integer, number of pixel rows that CCD camera has
@param n: Integer, number of pixel columns that CCD camera has
@param ext: String, file type of intensity file
@param numbit: Integer, number of bits in images

@return: Double mxn Matrix, averaged intensity matrix
@return: Double, standard deviation between N runs
%}

function [Istd, Iavg, Bp, I, B] = calcAvgI_v4(N, e, path, m, n, ext)
	
    currDir = pwd;		% Get cwd
    cd(path);			% Change directory to where the images are stored
    B = zeros(m,n,N);
    I = zeros(m,n,N);	% Declare matrix to house intensity profiles of all N runs
    for k = 1:N         % For each exposure ...
        f1 = strcat(upper(e),'_', int2str(k),ext); % Get filename of intensity profile
        f2 = strcat('BACK','_0_', int2str(k),ext); % Get filename of background intensity profile
        if ext == '.txt'
            Ipp = importdata(f1);           % Attempt to read in intensity data
            tempB = importdata(f2);         % Attempt to read in background data
        else
            Ipp = imread(f1);
            tempB = imread(f2);
        end
        Bp = (im2double(tempB));	% Convert into decimal grayscale values
        Ip = (im2double(Ipp));   % Convert into decimal grayscale values
        I(:,:,k) = Ip;
        B(:,:,k) = Bp;
    end
    Iavg1 = mean(I,3);						% Average matrix entries over N runs
    Bavg = mean(B,3);
    Iavg = Iavg1-Bavg;
    Iavg(Iavg<0)=0;                         % Ensure no negative numbers
    Istd = std(I,1,3);                      % Standard Deviation between N runs
    cd(currDir);							% Go back to working directory
end