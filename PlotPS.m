%{
Title: PlotPS.m
Author M. Runyon
Description: This script is function file that plots a 3D vector [s1,s2,s3]
             on the Poincare sphere.

    az: view point - azimuthal angle in degrees
    polx: view point - polar angle in degrees
%}

function PlotPS(s1,s2,s3, az, polx)
    r = 1;
    [x,y,z] = sphere(50);
    figure
    lightGrey = 0.85*[1 1 1]; % It looks better if the lines are lighter
    surf(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey);
    hold on
    plot3([0 s1],[0 s2],[0 s3],'b','linewidth',2);
    plot3([0 1],[0 0],[0 0], 'k' ) % for x-axis
    plot3([0 0],[0 1],[0 0], 'k') % for y-axis
    plot3([0 0],[0 0],[0 1], 'k') % for z-axis
    text(1.05,0,0,'S1','HorizontalAlignment','left','FontSize',10);
    text(0,1,0,'S2','HorizontalAlignment','left','FontSize',10);
    text(0,0,1.1,'S3','HorizontalAlignment','left','FontSize',10);
    hold off;
    xlabel('S1');ylabel('S2'),zlabel('S3');title('Poincare Sphere')
    view([az polx])
end