%{
Title: jVec2sVec.m
Author: M. Runyon
Description: This is a function file that converts a Jones vector to a
             Stokes vector. Assumes Jones vector is written in H-V basis.

%}

function [s0 s1 s2 s3] = jVec2sVec(j)

    s0 = j(1)*conj(j(1)) + j(2)*conj(j(2));
    s1 = j(1)*conj(j(1)) - j(2)*conj(j(2));
    s2 = 2*real(j(1)*conj(j(2)));
    s3 = -2*imag(j(1)*conj(j(2)));
    
end