%{
Title: calcFidelity
Author M. Runyon
Description: This script is  function file that returns the quantum
             fidelity of two density matrices, rho and sig.
%}

function F = calcFidelity(rho,sig)

   % sig = 1/sqrt(2).*[1 0;0 1];  %D
   % sig = 1/2.*[1 0;0 -1]; %A
   % sig = [1 0;0 0]; %H
   % sig = [0 0;0 1]; %V
    F = trace(sqrtm(sqrtm(rho)*sig*sqrtm(rho)));

end