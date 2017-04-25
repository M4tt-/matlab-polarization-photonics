
%{
Title: calcSNR_CCD.m
Author: M. Runyon
Description: This is a function file that determines the SNR in dB of a CCD
             detector. The input detector data should be in the form of a
             .txt file, tab delimited if applicable, where each entry
             corresponds to the intensity measured at a pixel of detector
             (CCD camera). Each text file should be named
             so that the first letter corresponds to the projection its for,
          followed by an underscore and number, where the number is ordered
          from 1 to N (see defn of N below). For example, if one wishes to
          save 3 files for intensity data corresponding to a horizontal
         projection, they would be named H_1.txt, H_2.txt, H_3.txt, etc. It
         is also assumed that N data captures of only the ambient light has
            been made. These should be named 'BACK_1.txt, BACK_2.txt' etc.
            
            The acceptable projection prefixes (before _X.txt) are:
            H - Horizontal
            V - Vertical
            D - Diagonal
            A - Antidiagonal
            R - Right hand circular
            L - Left hand circular
            
@param path: String, directory housing signal/background intensity files
@param N: Number of files to average over
@return SNR: Signal to noise ratio of detector

%}

function SNR = calcSNR_CCD(path, N)

    for i=1:N
        background_path=strcat(path,'BACK','_',int2str(i),'.txt'); %Background 
        signal_path1=strcat(path,'H','_',int2str(i),'.txt'); %Laser
        signal_path2=strcat(path,'V','_',int2str(i),'.txt'); %Laser

        background=importdata(background_path);
        signal1=importdata(signal_path1);
        signal2=importdata(signal_path2);
        
        background=im2double(background);
        signal=im2double(signal1)+im2double(signal2);

        totalBackground(i)=sum(sum(background));
        totalNoisySignal(i)=sum(sum(signal));
        totalPureSignal(i)=totalNoisySignal(i)-totalBackground(i);
    end

    meanNoisySignal = mean(totalNoisySignal);
    meanBackground = mean(totalBackground);
    meanPureSignal = meanNoisySignal - meanBackground;
    SNR = 10*log10(meanPureSignal/meanBackground);

end
