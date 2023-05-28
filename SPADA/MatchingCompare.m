function MatchingCompare(curve, PCC_result, uniqueArc, repeatIndex, optVariable)
% plot the new_curve to compare with the original curve
% 
% Input arguments:
% curve:
%   3-by-m matrix of the original curve data
% PCC_result:
%   structure array
%   struct('type', {1}, 'center', {[0,0,0]}, 'axis1', {[0,0,0]}, 'axis2', {[0,0,0]}, 'radius', {0}, 'angle', {0}, 'length', {0}, 'rotation', {0});
% uniqueArc: 
%   numberOfUniqueArc-by-2 matrix of the arclength and curvature of the unique arcs
% repeatIndex:
%   list of arcs with its index of the unique arcs
% optVariable:
%   optimal design parameters for the unique arcs

% 0bjective value: 0.0156: optVariable = [2.52037209777459,1.03624595483529,4186.01721039465,5.02719882050958,4.65658381049268,5.03645724753291,4.80413318584171]
ri = optVariable(1);
t = optVariable(2);
P = optVariable(3); % Pa
rm = [];
l = [];
optimal_angle = [];
N = [];

load NewANN.mat
for i = 1:(length(optVariable)-3)/2
    rm = [rm optVariable(3+(i-1)*2+1)];
    l = [l optVariable(3+(i-1)*2+2)];
    optimal_angle = [optimal_angle 2*sim(net,[ri, t, rm(end), l(end), P]')];
    N = [N round(uniqueArc(i,1)/l(i))];
end

new_curve = [];
j = 1; % arc_index
for i = 1:size(PCC_result,2)

    if isempty(new_curve)
        % if it is the first segment
        start_point = curve(:,1); % 3-by-1
    else
        start_point = new_curve(:,end); 
    end
    
    if PCC_result(i).type == 1
        % line segment
        new_curve = [new_curve start_point start_point - 2*PCC_result(i).axis1'];
    else
        % arc segment
        arc_angle = optimal_angle(repeatIndex(j)) * N(repeatIndex(j));
        arc_length = l(repeatIndex(j)) * N(repeatIndex(j));
        arc_radius = arc_length / arc_angle;
        
        arc_center = start_point - PCC_result(i).axis1';
        arc = [];
        theta = linspace(0, arc_angle, 20);
        for k = 1:length(theta)
            arc = [arc arc_center + PCC_result(i).axis1'*cos(theta(k)) + PCC_result(i).axis2'*sin(theta(k))];
        end
        new_curve = [new_curve arc];
        j = j+1;

    end
end

f = figure('Name', 'Matching Comparison',"NumberTitle", "off", "Position", [480, 250, 560, 420]);
% Get screen size
screenSize = get(groot, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
% Set the figure position
figWidth = 560;  % Set the width of the figure
figHeight = 420; % Set the height of the figure
f.Position =  [(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight];

hold on
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
origin = plot3(curve(1,:), curve(2,:), curve(3,:), 'k-', 'LineWidth', 2);
optimal = plot3(new_curve(1,:), new_curve(2,:), new_curve(3,:),  'g-', 'LineWidth', 5);
optimal.Color(4) = 0.5;
legend([origin,optimal],'Original Curve', 'Optimal Result')
hold off

end