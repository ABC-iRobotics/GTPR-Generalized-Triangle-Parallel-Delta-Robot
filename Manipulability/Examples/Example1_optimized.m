% GTPR Test

clear; clc;



% ABB
a1 = 0.251;
a2 = 0.251;   
a3 = 0.251;
e1 = 0.0451;
e2 = 0.0451;
e3 = 0.0451;
L1 = 0.3;
L2 = 0.3;
L3 = 0.3;
l1 = 0.8;
l2 = 0.8;
l3 = 0.8; 
alpha12 = deg2rad(120);
alpha13 = deg2rad(240);

% Case 2 - GTPR Triangle Changed
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

gtpr1 = GTPR(a1, a2, a3, L1, L2, L3, l1, l2, l3, e1, e2, e3, alpha12, alpha13);

TCP1 = gtpr1.directGeometry(thetas0);


% resolution of the workspace
% resolution = 0.05; %[m]
% resolution = 0.01; %[m]
resolution = 0.005; %[m]


figure(1);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

gtpr1.plotMechanism(depiction);

view(30,30);


% give the arm angles in which to plot sample the workspace
thetas_map1 = rand(1000,3)*deg2rad(141.5) - deg2rad(42);

%get the approximated workspace cuboid:
% ws_cube = gtpr.getCubeAboutRobot(resolution, 'thetas', thetas_map);%, 'edgePlot');

% map the workspace completely, based on the arm angles
[ws_cube1, PtsMapped, thetas_map1, coordsMapped1] = gtpr1.MapWorkspace(resolution, thetas_map1);


%%
% Skip these points, lines in plotting, see gtpr.bits
skipIdx_plotting = 2:49;

figure(2);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;
gtpr1.plotWSCube(ws_cube1, resolution, 'skipIds', skipIdx_plotting);
view(30,30);





%% after the workspace is done, get the manipulability etc, for a given slice parallel to one of the planes

xAxis1 = ws_cube1.grid.xAxis.points;
yAxis1 = ws_cube1.grid.yAxis.points;
zAxis1 = ws_cube1.grid.zAxis.points;

[xind1,yind1,zind1] = ind2sub(size(ws_cube1.grid.idxMap), 1:length(ws_cube1.manipulability.mu1));

% Select the plane in which to show the workspace, from the previous values:

xVal1 = 0;
yVal1 = 0;
zVal1 = -0.7;

% Give the plane to look at: can be: XY, XZ, YZ.
plane = "xz";

% get the points associated to the plane:
switch lower(plane)

case 'xy'

    [zVal1, zidx, zdistances, zminDist] = findClosestValue(zVal1, zAxis1);

    ctr = 1;
    for x1 = 1:length(xAxis1)
        for y1 = 1:length(zAxis1)
            TCPCoords1(:,ctr) = [xAxis1(x1);yAxis1(y1);zAxis1(zidx)];
            q1(:,ctr) = gtpr1.inverseGeometry(TCPCoords1(:,ctr));
            [mu11(ctr), mu21(ctr), mu31(ctr)] = gtpr1.getManipulabilityMeasures();
    
            if ~withinValueRange(TCPCoords1(1,ctr), xAxis1(x1), resolution) || ~withinValueRange(TCPCoords1(2,ctr), yAxis1(y1), resolution) || ~withinValueRange(TCPCoords1(3,ctr), zAxis1(zidx), resolution)
                cplxIdx1(ctr) = ctr;
            end

            ctr = ctr + 1;
        end
    end

    az = 0;
    el = 90;

case 'yz'

    [xVal1, xidx1, xdistances1, xminDist1] = findClosestValue(xVal1, xAxis1);

    ctr = 1;
    for y1 = 1:length(yAxis1)
        for z1 = 1:length(zAxis1)
            TCPCoords1(:,ctr) = [xAxis1(xidx1);yAxis1(y1);zAxis1(z1)];
            q1(:,ctr) = gtpr1.inverseGeometry(TCPCoords1(:,ctr));
            [mu11(ctr), mu21(ctr), mu31(ctr)] = gtpr1.getManipulabilityMeasures();
    
            if ~withinValueRange(TCPCoords1(1,ctr), xAxis1(xidx1), resolution) || ~withinValueRange(TCPCoords1(2,ctr), yAxis1(y1), resolution) || ~withinValueRange(TCPCoords1(3,ctr), zAxis1(z1), resolution)
                cplxIdx1(ctr) = ctr;
            end

            ctr = ctr + 1;
        end
    end

    az = 90;
    el = 0;

case 'xz'
	[yVal1, yidx1, ydistances1, yminDist1] = findClosestValue(yVal1, yAxis1);

    ctr = 1;
    for x1 = 1:length(xAxis1)
        for z1 = 1:length(zAxis1)
            TCPCoords1(:,ctr) = [xAxis1(x1);yAxis1(yidx1);zAxis1(z1)];
            q1(:,ctr) = gtpr1.inverseGeometry(TCPCoords1(:,ctr));
            [mu11(ctr), mu21(ctr), mu31(ctr)] = gtpr1.getManipulabilityMeasures();
    
            if ~withinValueRange(TCPCoords1(1,ctr), xAxis1(x1), resolution) || ~withinValueRange(TCPCoords1(2,ctr), yAxis1(yidx1), resolution) || ~withinValueRange(TCPCoords1(3,ctr), zAxis1(z1), resolution)
                cplxIdx1(ctr) = ctr;
            end

            ctr = ctr + 1;
        end
    end

    az = 0;
    el = 0;

end

[r1,c1] = find(imag(q1) ~= 0);

TCPCoords1(:,c1) = [];
q1(:,c1) = [];
mu11(:,c1) = [];
mu21(:,c1) = [];
mu31(:,c1) = [];

[ws_cube1, PtsMapped, thetas_map1, coordsMapped1] = gtpr1.MapWorkspace(resolution, q1');


% xAxisPlot1 = ws_cube1.grid.xAxis.points(sort(unique(ws_cube1.grid.mappedIndeces(:,1))));
% yAxisPlot1 = ws_cube1.grid.yAxis.points(sort(unique(ws_cube1.grid.mappedIndeces(:,2))));
% zAxisPlot1 = ws_cube1.grid.zAxis.points(sort(unique(ws_cube1.grid.mappedIndeces(:,3))));

xAxisPlot1 = round(min(xAxis1), 1) - 0.2 : 0.1: round(max(xAxis1), 1) + 0.2;
yAxisPlot1 = round(min(yAxis1), 1) - 0.2 : 0.1: round(max(yAxis1), 1) + 0.2;
zAxisPlot1 = round(min(zAxis1), 1) - 0.2 : 0.1: round(max(zAxis1), 1) + 0.2;

labelsx1 = xAxisPlot1;
labelsy1 = yAxisPlot1;
labelsz1 = zAxisPlot1;

xAxis1 = xAxisPlot1;
yAxis1 = yAxisPlot1;
zAxis1 = zAxisPlot1;

% labelsx1 = xAxisPlot1(1:2:end); % extract
% labelsy1 = yAxisPlot1(1:2:end); % extract
% labelsz1 = zAxisPlot1(1:2:end); % extract

% r = find(ws_cube.grid.mappedIndeces(:,2) ~= yidx);
% idxs = ws_cube.grid.mappedIndeces(r,:);

% xAxis1 = xAxis1(1:2:end);
% yAxis1 = yAxis1(1:2:end);
% zAxis1 = zAxis1(1:2:end);

% ws_cube.grid.idxMap(idxs(:,1), idxs(:,2),idxs(:,3)) = 0;
figure(1);
xticks(xAxis1);
yticks(yAxis1);
zticks(zAxis1);
set(gca, 'XTickLabel', string(xAxis1));
set(gca, 'YTickLabel', string(yAxis1));
set(gca, 'ZTickLabel', string(zAxis1));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig1.eps', 'epsc');
saveas(gcf, 'case1_fig1.png');

figure(2);
xticks(xAxis1);
yticks(yAxis1);
zticks(zAxis1);
set(gca, 'XTickLabel', string(xAxis1));
set(gca, 'YTickLabel', string(yAxis1));
set(gca, 'ZTickLabel', string(zAxis1));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig2.eps', 'epsc');
saveas(gcf, 'case1_fig2.png');


figure(3);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;
gtpr1.plotWSCube(ws_cube1, resolution, 'skipIds', skipIdx_plotting);
view(az,el);
view(az,el);
saveas(gcf, 'case1_fig3.eps', 'epsc');
saveas(gcf, 'case1_fig3.png');

set(gcf,'Color','w');



xticks(labelsx1);
yticks(labelsy1);
zticks(labelsz1);
set(gca, 'XTickLabel', string(labelsx1));
set(gca, 'YTickLabel', string(labelsy1));
set(gca, 'ZTickLabel', string(labelsz1));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');

% plot the manipulability stuff:

mu11 = ws_cube1.manipulability.mu1;
mu21 = ws_cube1.manipulability.mu2;
mu31 = ws_cube1.manipulability.mu3;

x1 = ws_cube1.grid.xAxis.points(ws_cube1.grid.mappedIndeces(:,1));
y1 = ws_cube1.grid.yAxis.points(ws_cube1.grid.mappedIndeces(:,2));
z1 = ws_cube1.grid.zAxis.points(ws_cube1.grid.mappedIndeces(:,3));

x1(r1) = [];
y1(r1) = [];
z1(r1) = [];

mu11(r1) = [];
mu21(r1) = [];
mu31(r1) = [];

mu1HeatmapRange1 = range(mu11);
mu2HeatmapRange1 = range(mu21);
mu3HeatmapRange1 = range(mu31);

c1 = colormap('parula');
% c = colormap('gray');

d1(:,1) = resample(c1(:,1), length(mu11), length(c1(:,1)));
d1(:,2) = resample(c1(:,2), length(mu11), length(c1(:,2)));
d1(:,3) = resample(c1(:,3), length(mu11), length(c1(:,3)));

c1 = d1;

c1 = rescale(c1,0,1);

mu1Map1 = rescale(movmean(sort(mu11), floor(length(c1)/10)), 0, 1);
mu2Map1 = rescale(movmean(sort(mu21), floor(length(c1)/10)), 0, 1);
mu3Map1 = rescale(movmean(sort(mu31), floor(length(c1)/10)), 0, 1);

mu1r = rescale(mu11, 0, 1);
mu2r = rescale(mu21, 0, 1);
mu3r = rescale(mu31, 0, 1);

figure(4);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

title("\mu_1 & \mu_2");

for i = 1:length(mu11)
    [~, row_mu1] = min(abs(mu1Map1-mu1r(i)));
    markerColors_mu1(i,:) = c1(row_mu1,:);
%     plot3(questionableTCPs(1,i), questionableTCPs(2,i), questionableTCPs(3,i), 'Color', markerColors_mu1(i,:), 'Marker', 'o');
end

scatter3(x1,y1,z1, 6, markerColors_mu1, 'filled');% 'k', 'filled');%, markerColors_mu1, 'filled');
set(gcf,'Color','w');
set(gca, 'clim', [min(mu11)/max(mu11) max(mu11)/max(mu11)]);
colorbar
view(az,el);
% set(gca, 'XTickLabel', h3.XTickLabel);
% set(gca, 'YTickLabel', h3.YTickLabel);
% set(gca, 'ZTickLabel', h3.ZTickLabel);

xticks(labelsx1);
yticks(labelsy1);
zticks(labelsz1);
set(gca, 'XTickLabel', string(labelsx1));
set(gca, 'YTickLabel', string(labelsy1));
set(gca, 'ZTickLabel', string(labelsz1));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig4.eps', 'epsc');
saveas(gcf, 'case1_fig4.png');

figure(5);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

title("\mu_2");


for i = 1:length(mu21)
    [~, row_mu2] = min(abs(mu2Map1-mu2r(i)));
    markerColors_mu2(i,:) = c1(row_mu2,:);
%     plot3(questionableTCPs(1,i), questionableTCPs(2,i), questionableTCPs(3,i), 'Color', markerColors_mu2(i,:), 'Marker', 'o');
end

scatter3(x1,y1,z1, 14, markerColors_mu2, 'filled');%'k', 'filled');%, markerColors_mu2, 'filled');

set(gca, 'clim', [min(mu21)/max(mu21) max(mu21)/max(mu21)]);
colorbar
view(az,el);
% set(gca, 'XTickLabel', h3.XTickLabel);
% set(gca, 'YTickLabel', h3.YTickLabel);
% set(gca, 'ZTickLabel', h3.ZTickLabel);

xticks(labelsx1);
yticks(labelsy1);
zticks(labelsz1);
set(gca, 'XTickLabel', string(labelsx1));
set(gca, 'YTickLabel', string(labelsy1));
set(gca, 'ZTickLabel', string(labelsz1));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig5.eps', 'epsc');
saveas(gcf, 'case1_fig5.png');


figure(6);
clf;
hold on;
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
axis equal;

title("\mu_3");

% title("Manipulability Measure");

for i = 1:length(mu31)
    [~, row_mu3] = min(abs(mu3Map1-mu3r(i)));
    markerColors_mu3(i,:) = c1(row_mu3,:);
%     plot3(questionableTCPs(1,i), questionableTCPs(2,i), questionableTCPs(3,i), 'Color', markerColors_mu3(i,:), 'Marker', 'o');
end

scatter3(x1,y1,z1, 14, markerColors_mu3, 'filled');

set(gca, 'clim', [min(mu31)/max(mu31) max(mu31)/max(mu31)]);
colorbar
view(az,el);
% set(gca, 'XTickLabel', h3.XTickLabel);
% set(gca, 'YTickLabel', h3.YTickLabel);
% set(gca, 'ZTickLabel', h3.ZTickLabel);


xticks(labelsx1);
yticks(labelsy1);
zticks(labelsz1);
set(gca, 'XTickLabel', string(labelsx1));
set(gca, 'YTickLabel', string(labelsy1));
set(gca, 'ZTickLabel', string(labelsz1));
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');


set(gcf,'Color','w');

saveas(gcf, 'case1_fig6.eps', 'epsc');
saveas(gcf, 'case1_fig6.png');



save("case1.mat");

































































