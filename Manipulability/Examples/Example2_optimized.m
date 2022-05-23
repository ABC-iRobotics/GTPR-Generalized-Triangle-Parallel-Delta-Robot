% GTPR Test

% clear; clc;



% ABB
% a1 = 0.251;
% a2 = 0.251;   
% a3 = 0.251;
% e1 = 0.0451;
% e2 = 0.0451;
% e3 = 0.0451;
% L1 = 0.3;
% L2 = 0.3;
% L3 = 0.3;
% l1 = 0.8;
% l2 = 0.8;
% l3 = 0.8; 
% alpha12 = deg2rad(120);
% alpha13 = deg2rad(240);

% Case 2 - GTPR Triangle Changed

a1 = 0.5;
a2 = 0.2; 
a3 = 0.4;
h = 0.1797;
e1 = a1*h;
e2 = a2*h;
e3 = a3*h;
L1 = 0.3;
L2 = 0.4;
L3 = 0.3;
l1 = 0.95;
l2 = 0.8;
l3 = 0.95; 
alpha12 = deg2rad(90);
alpha13 = deg2rad(-90);

% a1 = 0.5;
% a2 = 0.2; 
% a3 = 0.4;
% h = 0.1797;
% e1 = a1*h;
% e2 = a2*h;
% e3 = a3*h;
% L1 = 0.3;
% L2 = 0.4;
% L3 = 0.3;
% l1 = 0.95;
% l2 = 0.8;
% l3 = 0.95; 
% alpha12 = deg2rad(90);
% alpha13 = deg2rad(-90);

% a1 = 0.6;
% a2 = 0.5; 
% a3 = 0.251;
% h = 0.1797;
% e1 = a1*h;
% e2 = a2*h;
% e3 = a3*h;
% L1 = 0.4;
% L2 = 0.3;
% L3 = 0.4;
% l1 = 0.8;
% l2 = 0.9;
% l3 = 0.7; 
% alpha12 = deg2rad(90);
% alpha13 = deg2rad(-90);

% initial angle
thetas0 = deg2rad([0;0;0]); % [rad]

% how to plot the gtpr mechanism
depiction = 'real';
%depiction = 'theoretical';

gtpr2 = GTPR(a1, a2, a3, L1, L2, L3, l1, l2, l3, e1, e2, e3, alpha12, alpha13);

TCP2 = gtpr2.directGeometry(thetas0);


% resolution of the workspace
% resolution = 0.05; %[m]
% resolution = 0.01; %[m]
resolution = 0.005; %[m]

figure(7);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

gtpr2.plotMechanism(depiction);

view(30,30);


% give the arm angles in which to plot sample the workspace
thetas_map2 = rand(1000,3)*deg2rad(141.5) - deg2rad(42);

%get the approximated workspace cuboid:
% ws_cube = gtpr.getCubeAboutRobot(resolution, 'thetas', thetas_map);%, 'edgePlot');

% map the workspace completely, based on the arm angles
[ws_cube2, PtsMapped2, thetas_map2, coordsMapped2] = gtpr2.MapWorkspace(resolution, thetas_map2);

%%
% Skip these points, lines in plotting, see gtpr.bits
skipIdx_plotting = 2:49;

figure(8);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;
gtpr2.plotWSCube(ws_cube2, resolution, 'skipIds', skipIdx_plotting);
view(30,30);





%% after the workspace is done, get the manipulability etc, for a given slice parallel to one of the planes

xAxis2 = ws_cube2.grid.xAxis.points;
yAxis2 = ws_cube2.grid.yAxis.points;
zAxis2 = ws_cube2.grid.zAxis.points;

[xind2,yind2,zind2] = ind2sub(size(ws_cube2.grid.idxMap), 1:length(ws_cube2.manipulability.mu1));

% Select the plane in which to show the workspace, from the previous values:

xVal2 = 0;
yVal2 = 0;
zVal2 = -0.7;

% Give the plane to look at: can be: XY, XZ, YZ.
plane = "xz";

% get the points associated to the plane:
switch lower(plane)

case 'xy'

    [zVal2, zidx2, zdistances2, zminDist2] = findClosestValue(zVal2, zAxis2);

    ctr = 1;
    for x2 = 1:length(xAxis2)
        for y2 = 1:length(zAxis2)
            TCPCoords2(:,ctr) = [xAxis2(x2mdate);yAxis2(y2);zAxis2(zidx2)];
            q2(:,ctr) = gtpr2.inverseGeometry(TCPCoords2(:,ctr));
            [mu12(ctr), mu22(ctr), mu32(ctr)] = gtpr2.getManipulabilityMeasures();
    
            if ~withinValueRange(TCPCoords2(1,ctr), xAxis2(x2), resolution) || ~withinValueRange(TCPCoords2(2,ctr), yAxis2(y2), resolution) || ~withinValueRange(TCPCoords2(3,ctr), zAxis2(zidx2), resolution)
                cplxIdx2(ctr) = ctr;
            end

            ctr = ctr + 1;
        end
    end

    az = 0;
    el = 90;

case 'yz'

    [xVal2, xidx2, xdistances2, xminDist2] = findClosestValue(xVal2, xAxis2);

    ctr = 1;
    for y2 = 1:length(yAxis2)
        for z2 = 1:length(zAxis2)
            TCPCoords2(:,ctr) = [xAxis2(xidx2);yAxis2(y2);zAxis2(z2)];
            q2(:,ctr) = gtpr2.inverseGeometry(TCPCoords2(:,ctr));
            [mu12(ctr), mu22(ctr), mu32(ctr)] = gtpr2.getManipulabilityMeasures();
    
            if ~withinValueRange(TCPCoords2(1,ctr), xAxis2(xidx2), resolution) || ~withinValueRange(TCPCoords2(2,ctr), yAxis2(y2), resolution) || ~withinValueRange(TCPCoords2(3,ctr), zAxis2(z2), resolution)
                cplxIdx2(ctr) = ctr;
            end

            ctr = ctr + 1;
        end
    end

    az = 90;
    el = 0;

case 'xz'
	[yVal2, yidx2, ydistances2, yminDist2] = findClosestValue(yVal2, yAxis2);

    ctr = 1;
    for x2 = 1:length(xAxis2)
        for z2 = 1:length(zAxis2)
            TCPCoords2(:,ctr) = [xAxis2(x2);yAxis2(yidx2);zAxis2(z2)];
            q2(:,ctr) = gtpr2.inverseGeometry(TCPCoords2(:,ctr));
            [mu12(ctr), mu22(ctr), mu32(ctr)] = gtpr2.getManipulabilityMeasures();
    
            if ~withinValueRange(TCPCoords2(1,ctr), xAxis2(x2), resolution) || ~withinValueRange(TCPCoords2(2,ctr), yAxis2(yidx2), resolution) || ~withinValueRange(TCPCoords2(3,ctr), zAxis2(z2), resolution)
                cplxIdx2(ctr) = ctr;
            end

            ctr = ctr + 1;
        end
    end

    az = 0;
    el = 0;

end

[r2,c2] = find(imag(q2) ~= 0);

TCPCoords2(:,c2) = [];
q2(:,c2) = [];
mu12(:,c2) = [];
mu22(:,c2) = [];
mu32(:,c2) = [];

[ws_cube2, PtsMapped2, thetas_map2, coordsMapped2] = gtpr2.MapWorkspace(resolution, q2');

% xAxisPlot2 = ws_cube2.grid.xAxis.points(sort(unique(ws_cube2.grid.mappedIndeces(:,1))));
% yAxisPlot2 = ws_cube2.grid.yAxis.points(sort(unique(ws_cube2.grid.mappedIndeces(:,2))));
% zAxisPlot2 = ws_cube2.grid.zAxis.points(sort(unique(ws_cube2.grid.mappedIndeces(:,3))));


xAxisPlot2 = round(min(xAxis2), 1) - 0.2 : 0.1: round(max(xAxis2), 1) + 0.2;
yAxisPlot2 = round(min(yAxis2), 1) - 0.2 : 0.1: round(max(yAxis2), 1) + 0.2;
zAxisPlot2 = round(min(zAxis2), 1) - 0.2 : 0.1: round(max(zAxis2), 1) + 0.2;

labelsx2 = xAxisPlot2;
labelsy2 = yAxisPlot2;
labelsz2 = zAxisPlot2;

xAxis2 = xAxisPlot2;
yAxis2 = yAxisPlot2;
zAxis2 = zAxisPlot2;

% labelsx2 = xAxisPlot2(1:2:end); % extract
% labelsy2 = yAxisPlot2(1:2:end); % extract
% labelsz2 = zAxisPlot2(1:2:end); % extract

% r = find(ws_cube.grid.mappedIndeces(:,2) ~= yidx);
% idxs = ws_cube.grid.mappedIndeces(r,:);

% xAxis2 = xAxis2(1:2:end);
% yAxis2 = yAxis2(1:2:end);
% zAxis2 = zAxis2(1:2:end);

% ws_cube.grid.idxMap(idxs(:,1), idxs(:,2),idxs(:,3)) = 0;
figure(7);
xticks(xAxis2);
yticks(yAxis2);
zticks(zAxis2);
set(gca, 'XTickLabel', string(xAxis2));
set(gca, 'YTickLabel', string(yAxis2));
set(gca, 'ZTickLabel', string(zAxis2));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig1.eps', 'epsc');
saveas(gcf, 'case1_fig1.png');

figure(8);
xticks(xAxis2);
yticks(yAxis2);
zticks(zAxis2);
set(gca, 'XTickLabel', string(xAxis2));
set(gca, 'YTickLabel', string(yAxis2));
set(gca, 'ZTickLabel', string(zAxis2));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig2.eps', 'epsc');
saveas(gcf, 'case1_fig2.png');


figure(9);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;
gtpr2.plotWSCube(ws_cube2, resolution, 'skipIds', skipIdx_plotting);
view(az,el);
view(az,el);
saveas(gcf, 'case1_fig3.eps', 'epsc');
saveas(gcf, 'case1_fig3.png');

set(gcf,'Color','w');



xticks(labelsx2);
yticks(labelsy2);
zticks(labelsz2);
set(gca, 'XTickLabel', string(labelsx2));
set(gca, 'YTickLabel', string(labelsy2));
set(gca, 'ZTickLabel', string(labelsz2));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');

% plot the manipulability stuff:

mu12 = ws_cube2.manipulability.mu1;
mu22 = ws_cube2.manipulability.mu2;
mu32 = ws_cube2.manipulability.mu3;

x2 = ws_cube2.grid.xAxis.points(ws_cube2.grid.mappedIndeces(:,1));
y2 = ws_cube2.grid.yAxis.points(ws_cube2.grid.mappedIndeces(:,2));
z2 = ws_cube2.grid.zAxis.points(ws_cube2.grid.mappedIndeces(:,3));

x2(r2) = [];
y2(r2) = [];
z2(r2) = [];

mu12(r2) = [];
mu22(r2) = [];
mu32(r2) = [];

mu1HeatmapRange2 = range(mu12);
mu2HeatmapRange2 = range(mu22);
mu3HeatmapRange2 = range(mu32);

c2 = colormap('parula');
% c = colormap('gray');

d2(:,1) = resample(c2(:,1), length(mu12), length(c2(:,1)));
d2(:,2) = resample(c2(:,2), length(mu12), length(c2(:,2)));
d2(:,3) = resample(c2(:,3), length(mu12), length(c2(:,3)));

c2 = d2;

c2 = rescale(c2,0,1);

mu1Map2 = rescale(movmean(sort(mu12), floor(length(c2)/10)), 0, 1);
mu2Map2 = rescale(movmean(sort(mu22), floor(length(c2)/10)), 0, 1);
mu3Map2 = rescale(movmean(sort(mu32), floor(length(c2)/10)), 0, 1);

mu1r = rescale(mu12, 0, 1);
mu2r = rescale(mu22, 0, 1);
mu3r = rescale(mu32, 0, 1);

figure(10);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

title("\mu_1 & \mu_2");

for i = 1:length(mu12)
    [~, row_mu1] = min(abs(mu1Map2-mu1r(i)));
    markerColors_mu12(i,:) = c2(row_mu1,:);
%     plot3(questionableTCPs(1,i), questionableTCPs(2,i), questionableTCPs(3,i), 'Color', markerColors_mu1(i,:), 'Marker', 'o');
end

scatter3(x2,y2,z2, 14, markerColors_mu12, 'filled');% 'k', 'filled');%, markerColors_mu1, 'filled');
set(gcf,'Color','w');
set(gca, 'clim', [min(mu12)/max(mu12) max(mu12)/max(mu12)]);
colorbar
view(az,el);
% set(gca, 'XTickLabel', h3.XTickLabel);
% set(gca, 'YTickLabel', h3.YTickLabel);
% set(gca, 'ZTickLabel', h3.ZTickLabel);

xticks(labelsx2);
yticks(labelsy2);
zticks(labelsz2);
set(gca, 'XTickLabel', string(labelsx2));
set(gca, 'YTickLabel', string(labelsy2));
set(gca, 'ZTickLabel', string(labelsz2));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case2_fig4.eps', 'epsc');
saveas(gcf, 'case2_fig4.png');

figure(11);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

title("\mu_2");


for i = 1:length(mu22)
    [~, row_mu2] = min(abs(mu2Map2-mu2r(i)));
    markerColors_mu22(i,:) = c2(row_mu2,:);
%     plot3(questionableTCPs(1,i), questionableTCPs(2,i), questionableTCPs(3,i), 'Color', markerColors_mu2(i,:), 'Marker', 'o');
end

scatter3(x2,y2,z2, 14, markerColors_mu22, 'filled');%'k', 'filled');%, markerColors_mu2, 'filled');

set(gca, 'clim', [min(mu22)/max(mu22) max(mu22)/max(mu22)]);
colorbar
view(az,el);
% set(gca, 'XTickLabel', h3.XTickLabel);
% set(gca, 'YTickLabel', h3.YTickLabel);
% set(gca, 'ZTickLabel', h3.ZTickLabel);

xticks(labelsx2);
yticks(labelsy2);
zticks(labelsz2);
set(gca, 'XTickLabel', string(labelsx2));
set(gca, 'YTickLabel', string(labelsy2));
set(gca, 'ZTickLabel', string(labelsz2));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case2_fig5.eps', 'epsc');
saveas(gcf, 'case2_fig5.png');


figure(12);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

title("\mu_3");

% title("Manipulability Measure");

for i = 1:length(mu32)
    [~, row_mu3] = min(abs(mu3Map2-mu3r(i)));
    markerColors_mu32(i,:) = c2(row_mu3,:);
%     plot3(questionableTCPs(1,i), questionableTCPs(2,i), questionableTCPs(3,i), 'Color', markerColors_mu3(i,:), 'Marker', 'o');
end

scatter3(x2,y2,z2, 14, markerColors_mu32, 'filled');

set(gca, 'clim', [min(mu32)/max(mu32) max(mu32)/max(mu32)]);
colorbar
view(az,el);
% set(gca, 'XTickLabel', h3.XTickLabel);
% set(gca, 'YTickLabel', h3.YTickLabel);
% set(gca, 'ZTickLabel', h3.ZTickLabel);


xticks(labelsx2);
yticks(labelsy2);
zticks(labelsz2);
set(gca, 'XTickLabel', string(labelsx2));
set(gca, 'YTickLabel', string(labelsy2));
set(gca, 'ZTickLabel', string(labelsz2));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');


set(gcf,'Color','w');

saveas(gcf, 'case2_fig6.eps', 'epsc');
saveas(gcf, 'case2_fig6.png');


save("case2.mat");

































































