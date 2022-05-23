clear; clc;

[X,Y,Z] = ellipsoid(0,0,0,0.8,0.6,1);

figure(1);
clf;
grid off;
axis equal;

plot3(0,0,0,'ko');

s = surf(X,Y,Z, 'FaceAlpha', 0.1, 'EdgeColor',[0 0 0], 'FaceColor', [0 0 0], 'LineWidth', 0.0001, 'LineStyle','none');


hold on;

cx = ellipse3D(0.8, 0.6, 0,0,0);

cy = ellipse3D(1, 0.6, 0,0,0, 300, 0,pi/2,0);

cz = ellipse3D(1, 0.8, 0,0,0, 300, pi/2,pi/2,0);

plot3(cx(1,:), cx(2,:), cx(3,:), 'k');
plot3(cy(1,:), cy(2,:), cy(3,:), 'k');
plot3(cz(1,:), cz(2,:), cz(3,:), 'k');


dirx = [0.8*1.5; 0; 0];
diry = [0; 0.6*1.5; 0];
dirz = [0; 0; 1*1.5];

quiver3(0,0,0, dirx(1), dirx(2), dirx(3), 'r:');
quiver3(0,0,0, diry(1), diry(2), diry(3), 'm:');
quiver3(0,0,0, dirz(1), dirz(2), dirz(3), 'b:');

% quiver3(norm(dirx),0,0, 0.8/2, 0, 0, 'r');
% quiver3(0,norm(diry),0, 0, 0.6/2, 0, 'm');
% quiver3(0,0,norm(dirz), 0, 0, 1/2, 'b');


quiver3D([norm(dirx),0,0], [0.2, 0, 0], 'r', 0.8);
quiver3D([0,norm(diry),0], [0, 0.2, 0], 'm', 0.8);
quiver3D([0,0,norm(dirz)], [0, 0, 0.2], 'b', 0.8);

text(norm(dirx),0.15,0, '\Lambda_1', 'FontSize', 20, 'Color', 'r');
text(0.1,norm(diry),0, '\Lambda_2', 'FontSize', 20, 'Color', 'm');
text(0.05,0,norm(dirz)+0.1, '\Lambda_3', 'FontSize', 20, 'Color', 'b');


text(norm(dirx)/2-0.1,0.1,0, '$\sqrt{\lambda_1}$', 'FontSize', 15, 'Color', 'r', 'Interpreter', 'latex');
text(0.1,norm(diry)/2-0.2,0, '$\sqrt{\lambda_2}$', 'FontSize', 15, 'Color', 'm', 'Interpreter', 'latex');
text(0,0,norm(dirz)/2, '$\sqrt{\lambda_3}$', 'FontSize', 15, 'Color', 'b', 'Interpreter', 'latex');

set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
set(gca,'Visible','off');
set(gcf,'color','w');

view(30,30);