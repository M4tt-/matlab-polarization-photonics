%{
Title: gaussFit.m
Author: M. Runyon
Description: This script is a function file that finds the gaussian
             parameters mu (average) and sigma (standard deviation) of a
             one dimensional normal distribution.
%}

function [mu sigma] = gaussFit(A)
    x0 = 0;         % initialize x0 to zero
    x02 = 0;        % initialize x02 to zero
    for k = 1:length(A)
        x0 = x0 + k.*A(k)./sum(A);       % expectation value of x
        x02 = x02 + (k.^2).*A(k)./sum(A); % expectation value of x2
    end
    
    sigma = sqrt(x02-x0.^2);
    mu = x0;
end
