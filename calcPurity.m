%{
Title: calcPurity.m
Author: M. Runyon
Description: This script is a function file calculates the quantum
             purity of an input density matrix.
%}

function P = calcPurity(rho)
    
    P = trace(rho^2);

end