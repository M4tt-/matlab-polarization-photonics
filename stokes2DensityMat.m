%{
Title: Stokes2DensityMat
Author M. Runyon
Description: This script is  function file that takes a Stokes vectorand
             outputs a density matrix
%}

function rho = stokes2DensityMat(S)

    I = [1 0;0 1];
    sig3 = [0 1;1 0];
    sig2 = [0 -1i; 1i 0];
    sig1 = [1 0;0 -1];

    rho = 1/2.*(I + (S(3)*sig3)+(S(2)*sig2)+(S(1)*sig1));

end