curve_1_x = zeros(1,50);
curve_1_y = 0:50*1.5/49:50*1.5;

R_2 = 20;
curve_2_center_x = R_2;
curve_2_center_y = 50*1.5;
curve_2_angle = pi:(0-pi)/49:0;
curve_2_x = curve_2_center_x + R_2*cos(curve_2_angle);
curve_2_y = curve_2_center_y + R_2*sin(curve_2_angle);



curve_3_center_x = 0;
curve_3_center_y = 50*1.5;
R_3 = 40;
curve_3_angle = 0:(-pi/4)/49:-pi/4;
curve_3_x = curve_3_center_x + R_3*cos(curve_3_angle);
curve_3_y = curve_3_center_y + R_3*sin(curve_3_angle);


R_4 = sqrt((50/4*3 - curve_3_x(end))^2 + (50/4*3 - curve_3_y(end))^2);
curve_4_center_x = curve_3_x(end) + R_4*sqrt(2)/2;
curve_4_center_y = curve_3_y(end) - R_4*sqrt(2)/2;
curve_4_angle = pi*3/4:(pi*1/4+pi/180*25)/49:pi+pi/180*25;
curve_4_x = curve_4_center_x + R_4*cos(curve_4_angle);
curve_4_y = curve_4_center_y + R_4*sin(curve_4_angle);

curve_5_x = curve_4_x(end):(curve_4_center_y*tan(pi/180*25))/49:curve_4_x(end)+curve_4_center_y*tan(pi/180*25);
curve_5_y = curve_4_y(end):(0-curve_4_y(end))/49:0;

R_x = [curve_1_x curve_2_x curve_3_x curve_4_x curve_5_x];
R_y = [curve_1_y curve_2_y curve_3_y curve_4_y curve_5_y];
R_z = zeros(1,length(R_x));
curve = [R_x; R_y; R_z];

figure
hold on
axis equal
plot3(R_x, R_y, R_z)
hold off

save('new_letter_R.mat', 'curve')