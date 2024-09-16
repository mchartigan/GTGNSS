function plotLunarOrbit(ts,sats,plot_frame,name)
%PLOTLUNARORBIT Plots the trajectory of a lunar orbit, given the 3- or
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
R_me = 1738.1;          % km, moon equatorial radius
R_mp = 1736;            % km, moon polar radius
[Imoon, ~] = imread("Moon_HermesCelestiaMotherlode.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, R_me, R_me, R_mp);

% Rotate moon from MOON_ME frame to plot_frame
T = cspice_pxform('MOON_ME', plot_frame, ts(end));
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

% put a little star on the south pole
lsp = T * [0;0;-R_mp];
scatter3(lsp(1), lsp(2), lsp(3), 100, "red", "filled", "pentagram");

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

