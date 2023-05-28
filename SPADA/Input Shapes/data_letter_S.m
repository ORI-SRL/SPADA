%% ordered points of letter S
upper_center_x = 25; %mm
upper_center_y = 25*3;
upper_angle = 0:pi*3/2/89:pi*3/2;
upper_x = upper_center_x + 25*cos(upper_angle);
upper_y = upper_center_y + 25*sin(upper_angle);
lower_center_x = 25;
lower_center_y = 25;
lower_angle = pi/2:-pi*3/2/89:-pi;
lower_x = lower_center_x + 25*cos(lower_angle);
lower_y = lower_center_y + 25*sin(lower_angle);
s_x = [upper_x lower_x];
s_y = [upper_y lower_y];
s_z = zeros(1,length(s_x));
figure(1)
hold on
axis equal
plot3(s_x,s_y,s_z,'kx')
hold off
curve = [s_x;s_y;s_z];
save('new_letter_S.mat', 'curve')