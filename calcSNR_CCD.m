%{
Title: calcSNR_CCD.m
Author: M. Runyon
Description: This is a function file that determines the SNR in dB of a 
             detector. The input detector data should be in the form of a
             .txt file, tab delimited if applicable, where each entry
             corresponds to the intensity measured at a pixel of detector
             (CCD camera). Each text file should be named
             so that the first letter corresponds to the projection its for
          followed by an underscore and number, where the number is ordered
          from 1 to N (see defn of N below). For example, if one wishes to
          save 3 files for intensity data corresponding to a horizontal
         projection, they would be named H_1.txt, H_2.txt, H_3.txt, etc. It
        is also assumed that N data captures of only the ambient light has
            been made. These should be named 'BACK_1.txt, BACK_2.txt' etc.

@param path: String, directory housing signal/background intensity files
@param e: The alphabetical prefix describing which intensity projection
           is to be selected.
@param N: Number of files to average over

@return SNR: Signal to noise ratio of detector
%}

function SNR = calcSNR_CCD(path,e, N)

    for i=1:N
        background_path=strcat(path,'BACK','_',int2str(i),'.txt'); %Background 
        signal_path=strcat(path,e,'_',int2str(i),'.txt'); %Laser

        background=importdata(background_path);
        signal=importdata(signal_path);
        
        background=im2double(background);
        signal=im2double(signal);

        totalBackground(i)=sum(sum(background));
        totalNoisySignal(i)=sum(sum(signal));
        totalPureSignal(i)=totalNoisySignal(i)-totalBackground(i);
    end

    meanNoisySignal = mean(totalNoisySignal);
    meanBackground = mean(totalBackground);
    meanPureSignal = meanNoisySignal - meanBackground;
    SNR = 10*log10(meanPureSignal/meanBackground);

end