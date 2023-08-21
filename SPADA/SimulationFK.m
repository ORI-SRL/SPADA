function FK = SimulationFK(assemblyData, repeatIndex, uniqueUnitResult)
% Output the position of the designed actuator by using FEM to simulate unique modules and 
% applying forward kinematics to get the whole actuator position
%
% Input arguments:
% assemblyData:
%    n-by-5 matrix of the modular design parameters: phi, r, t, R, l 
% repeatIndex:
%   list of arcs with its index of the unique arcs e.g.: [1 2 3 2 2 3]
% uniqueUnitResult:
%   list of the angular deflection of unique modules
% 
% Output arguments:
% FK:
%   figure of the deformed actuator's position 

FK = figure('Name', 'Simulated Bellow SPA', 'Visible','on', "NumberTitle", "off");

% Get screen size
screenSize = get(groot, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
% Set the figure position
figWidth = 560;  % Set the width of the figure
figHeight = 420; % Set the height of the figure
FK.Position =  [(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight];


hold on

xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
view(30,30)
plot3(0,0,0,'ko','MarkerSize',5);
cm = jet;
T = eye(4);
X = [];
Y = [];
Z = [];
for i = 1:size(assemblyData, 1)
    % plot original position
    ori_x = zeros(1,10);
    ori_y = zeros(1,10);
    ori_z = -(sum(assemblyData(1:i,end)) - assemblyData(i,end)):-assemblyData(i,end)/9:-sum(assemblyData(1:i,end));
    ori_position = plot3(ori_x, ori_y, ori_z, '-','color', [0,0,0,0.5], 'LineWidth', 2);
    plot3(ori_x(end), ori_y(end), ori_z(end), 'ko', 'MarkerSize',5);
    
    if i == 1
        phi = assemblyData(i,1);
    else
        phi = assemblyData(i,1) - assemblyData(i-1,1);
    end
    theta = uniqueUnitResult(repeatIndex(i));
    new_l = (assemblyData(i,end)/theta + assemblyData(i,2) + assemblyData(i,3)/2)*theta;
    r = new_l/theta;
    Tz = [cos(phi), -sin(phi), 0, 0; sin(phi), cos(phi), 0, 0 ; 0, 0, 1, 0; 0, 0, 0, 1];
    Ty = [cos(theta), 0, sin(theta), -r*(1-cos(theta)); 0, 1, 0, 0; -sin(theta), 0, cos(theta), -r*sin(theta); 0, 0, 0, 1];
    nTz = [cos(-phi), -sin(-phi), 0, 0; sin(-phi), cos(-phi), 0, 0 ; 0, 0, 1, 0; 0, 0, 0, 1];
    
    center_x = -r;
    center_y = 0;
    center_z = 0;
    angle = 0:-theta/9:-theta;
    x = center_x + r*cos(angle);
    y = center_y + zeros(1,length(x));
    z = center_z + r*sin(angle);
    
    
    T = T*Tz;
    arc = T*[x;y;z;ones(1,length(x))];
    X = [X arc(1,:)];
    Y = [Y arc(2,:)];
    Z = [Z arc(3,:)];
    T = T*Ty;
    
    deform_position = plot3(arc(1,:), arc(2,:), arc(3,:), '-', 'color', cm(1+((i-1)*50)-250*floor((i-1)*50/250),:), 'LineWidth', 2);
    plot3(arc(1,end), arc(2,end), arc(3,end), 'ko','MarkerSize',5);
    
end

grid on
box on
axis equal
hold off


end

