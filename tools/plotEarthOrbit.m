function plotEarthOrbit(ts,sats,plot_frame,name)
%PLOTEARTHORBIT Plots the trajectory of an Earth orbit, given the 3- or
%6-state output of ODE45 (or similar propagators).
%   Inputs:
%    - ts; time steps
%    - sats; satellite positions over time
%    - plot_frame; frame that sats is represented in for SPICE
%    - name; plot title
arguments
    ts          (1,:)   double
    sats        (:,:,:) double
    plot_frame  (1,:)   char
    name        (1,1)   string
end

figure();
% Display moon in trajectory plot
R_e = 6378.137;         % km, earth equatorial radius
R_p = 6356.752;         % km, earth polar radius
[Imoon, ~] = imread("ModifiedBlueMarble.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, R_e, R_e, R_p);

% Rotate moon from ITRF frame to plot_frame
T = cspice_pxform('ITRF93', plot_frame, ts(end));
for j=1:size(xx,1)
    for k=1:size(xx,2)
        % -z to flip image
        tmp = T * [xx(j,k); yy(j,k); -zz(j,k)];
        xx(j,k) = tmp(1); yy(j,k) = tmp(2); zz(j,k) = tmp(3);
    end
end

globe = surf(xx, yy, zz);
set(globe, 'FaceColor', 'texturemap', 'CData', Imoon, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');
hold on;

% plot user trajectory for same time frame
for k=1:size(sats,3)
    plot3(sats(:,1,k), sats(:,2,k), sats(:,3,k), "LineWidth", 1.5, ...
        "Marker", "diamond", "MarkerIndices", length(ts));
end

grid on; axis equal;
SUB = strsplit(plot_frame,"_");
SUB = SUB(end);
xlabel("x_{"+SUB+"} (km)");
ylabel("y_{"+SUB+"} (km)");
zlabel("z_{"+SUB+"} (km)");
title(name);
end

