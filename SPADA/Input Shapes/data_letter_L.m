%% ordered points of letter L
curve_1_x = zeros(1,50);
curve_1_y = 100:-80/49:20;
center_x = 20; %mm
center_y = 20;
angle = pi:pi/2/49:pi*3/2;
curve_2_x = center_x + 20*cos(angle);
curve_2_y = center_y + 20*sin(angle);
curve_3_x = 20:30/49:50;
curve_3_y = zeros(1,50);
l_x = [curve_1_x curve_2_x curve_3_x];
l_y = [curve_1_y curve_2_y curve_3_y];
l_z = zeros(1,length(l_x));
figure(1)
hold on
axis equal
plot3(l_x,l_y,l_z,'kx')
hold off
curve = [l_x;l_y;l_z];
save('new_letter_L.mat', 'curve')