% Define helix parameters
r_helix = 25;
h = 50;
n_turns = 3;
n_points_per_turn = 100/n_turns; % number of points per turn

figure;
hold on;
axis equal;
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
t = linspace(0, n_turns*2*pi, n_turns*n_points_per_turn);
z_helix = - h*(t/(2*pi*n_turns));
x_helix = r_helix.*cos(t);
y_helix = r_helix.*sin(t);
plot3(x_helix, y_helix, z_helix, 'LineWidth', 2);
curve = [x_helix; y_helix; z_helix];
hold off
save('new_helix.mat', 'curve')