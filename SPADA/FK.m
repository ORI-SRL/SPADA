ri = 5;
t = 1.5;
rm = 8;
l = 10;
N = 8;

angle_FEM_half_hyperelastic = 2*[0.00145080000000000,0.0116840000000000,0.0212660000000000,0.0302700000000000,0.0387600000000000,0.0467900000000000,0.0544050000000000,0.0616470000000000,0.0685490000000000,0.0751410000000000,0.0814510000000000,0.0875020000000000,0.0933140000000000,0.0989060000000000,0.104300000000000,0.109500000000000,0.114520000000000,0.119380000000000,0.124090000000000,0.128660000000000,0.133100000000000];
phi = pi/4;
Tz = [cos(phi), -sin(phi), 0, 0; sin(phi), cos(phi), 0, 0 ; 0, 0, 1, 0; 0, 0, 0, 1];
theta = angle_FEM_half_hyperelastic(end);
r = l/theta;

Ty = [cos(theta), 0, sin(theta), -r*(1-cos(theta)); 0, 1, 0, 0; -sin(theta), 0, cos(theta), -r*sin(theta); 0, 0, 0, 1];
figure(1)
hold on
cm = jet;
for i = 1:N
    center_x = -r;
    center_y = 0;
    center_z = 0;
    angle = 0:-theta/9:-theta;
    x = center_x + r*cos(angle);
    y = center_y + zeros(1,length(x));
    z = center_z + r*sin(angle);
    
    T = (Ty*Tz)^(i-1);
    arc = T*[x;y;z;ones(1,length(x))];
    plot3(arc(1,:), arc(2,:), arc(3,:), '-', 'color', cm(1+(i-1)*30,:), 'LineWidth', 2);
end
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on
box on
axis equal
hold off