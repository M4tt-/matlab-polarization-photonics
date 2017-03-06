%{
Title: sVec2ell.m
Author: M. Runyon
Description: This is a function file that converts a Stokes vectors to the
             a polarization ellipse characterized by parameters psi, chi:
             https://en.wikipedia.org/wiki/Stokes_parameters

             psi: orientation of ellipse (polar angle in xy plane) [0,pi]
             chi: ellipticity/eccentricity of ellipse [-pi/4,pi/4]
             a: semi-major axis
             b: semi-minor axis
%}

function [psi, chi, a, b] = sVec2ell(s)

    I = s(1);
    p = sqrt(s(2)^2+s(3)^2+s(4)^2);
    psi = mod(0.5*atan2(s(3),s(2)),pi);
    chi = 0.5*atan2(s(4),sqrt(s(2)^2)+s(3)^2);
    theta = atan2(sec(chi),tan(2*psi));
    a = sqrt((1+sqrt(1-(sin(2*psi)^2)*(sin(chi))^2))/2);
    b = sqrt((1-sqrt(1-(sin(2*psi)^2)*(sin(chi))^2))/2);
    
end