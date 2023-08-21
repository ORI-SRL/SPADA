function [x, fval] = GAOpt(ArcSegment, BD)
% genetic algorithm to find the optimal design parameters for shape matching
% 
% Input argument:
% ArcSegment:
%    N-by-2 matrix of the unique arc length and curvature
% BD:
%    2-by-2 matrix of the outer radius and pressure boundary conditions
% 
% output argument:
% x:
%    1-by-(3+2N) matrix of the optimal design parameters
% fval:
%    the optimal value for the objective function

Segment_num = size(ArcSegment,1);
Segment_length = ArcSegment(:,1); %mm
Segment_curvature = ArcSegment(:,2); % 1/mm
[max_curvature, ~] = max(Segment_curvature);
[min_length, ~] = min(Segment_length);

% x = [ri, t, P, rm1, l1, rm2, l2, ...];

ro_lb = BD(1,1); 
ro_ub = BD(1,2);
P_lb = max(0, BD(2,1));
P_ub = min(100000, BD(2,2));

A = zeros(8+9*Segment_num, 3+2*Segment_num);
b = zeros(8+9*Segment_num, 1);
lb = zeros(1,3+2*Segment_num);
ub = zeros(1,3+2*Segment_num);
% bounds and linear constraint
% ri >= 0 is -x(1) <= 0
% ri <= min(1/max_curvature, ro_ub) is x(1) <= min(1/max_curvature, ro_ub)
% t >= 0 is -x(2) <= 0
% t <= min_length/4 is x(2) <= min_length/4
% P >= max(0, BD(2,1)) is -x(3) <= -max(0, BD(2,1))
% P <= min(100000, BD(2,2)) is x(3) <= min(100000, BD(2,2))
% t >= ri/4 is x(1) - 4x(2) <= 0
% t < ri/2 is -x(1) + 2x(2) < 0
A(1,1) = -1;
A(2,1) = 1;
b(2) = min(1/max_curvature, ro_ub);
A(3,2) = -1; 
b(3) = -1.5; % to make the wall thickness larger than 1.5mm
A(4,2) = 1;
b(4) = min_length/4;
A(5,3) = -1;
b(5) = -max(0, BD(2,1));
A(6,3) = 1;
b(6) = min(100000, BD(2,2));
A(7,1) = 1;
A(7,2) = -4;
A(8,1) = -1;
A(8,2) = 2;
lb(2) = 1;
ub(1) = min(1/max_curvature, ro_ub);
ub(2) = min_length/4;
ub(3) = min(100000,P_ub);


for i = 1:Segment_num
% rm_lb = max(ri + t, (ri + ro_lb)/2);
% rm_ub = min([(1/Segment_curvature(i)+ ri)/2, ri*2, (ri + ro_ub)/2]);
% l_lb = 4*t;
% l_ub = min([Segment_length(i), 4*ri,4*(rm-ri)]);

% rm >= ri + t is x(1) + x(2) - x(3+2i-1) <= 0
% rm >= (ri + ro_lb)/2 is x(1) - 2x(3+2i-1) <= -ro_lb
% rm <= (1/Segment_curvature(i)+ ri)/2 is -x(1) + 2x(3+2i-1) <=
% 1/Segment_curvature(i)
% rm <= ri*2 is -2x(1) + x(3+2i-1) <= 0
% rm <= (ri + ro_ub)/2 is -x(1) + 2x(3+2i-1) < ro_ub
% l > 2t is 3x(2) - x(3+2i) <= 0
% l <= Segment_length(i) is x(3+2i) <= Segment_length(i)
% l <= 4*ri is -4x(1) + x(3+2i) <= 0
% l <= 4*(rm-ri) is 4x(1) - 4x(3+2i-1) + x(3+2i) <= 0
    A(8+9*(i-1)+1,1) = 1;
    A(8+9*(i-1)+1,2) = 1;
    A(8+9*(i-1)+1,3+2*(i-1)+1) = -1;
    A(8+9*(i-1)+2,1) = 1;
    A(8+9*(i-1)+2,3+2*(i-1)+1) = -2;
    b(8+9*(i-1)+2) = -ro_lb;
    A(8+9*(i-1)+3,1) = -1;
    A(8+9*(i-1)+3,3+2*(i-1)+1) = 2;
    b(8+9*(i-1)+3) = 1/Segment_curvature(i);
    A(8+9*(i-1)+4,1) = -2;
    A(8+9*(i-1)+4,3+2*(i-1)+1) = 1;
    A(8+9*(i-1)+5,1) = -1;
    A(8+9*(i-1)+5,3+2*(i-1)+1) = 2;
    b(8+9*(i-1)+5) = ro_ub;
    A(8+9*(i-1)+6,2) = 3; % l > 2t is 3x(2) - x(3+2i) <= 0
    A(8+9*(i-1)+6,3+2*(i-1)+2) = -1;
    A(8+9*(i-1)+7,3+2*(i-1)+2) = 1;
    b(8+9*(i-1)+7) = Segment_length(i);
    A(8+9*(i-1)+8,1) = -4;
    A(8+9*(i-1)+8,3+2*(i-1)+2) = 1;
    A(8+9*(i-1)+9,1) = 4;
    A(8+9*(i-1)+9,3+2*(i-1)+1) = -4;
    A(8+9*(i-1)+9,3+2*(i-1)+2) = 1;
    lb(3+2*(i-1)+1)= ro_lb/2;
    ub(3+2*(i-1)+1)= min(1/Segment_curvature(i), ro_ub);
    ub(3+2*(i-1)+2)= min([Segment_length(i), Segment_length(i)/(Segment_length(i)*Segment_curvature(i)/(pi/4)), 2*ro_ub]);
end

Aeq = [];
beq = [];
nonlcon = [];
intcon = [];

fun = @(x)GAObj(x, Segment_length, Segment_curvature);
options= optimoptions('ga', 'MaxGenerations', 500, 'OutputFcn', @(options,state,flag)GAoutputfun(options, state, flag, Segment_length, Segment_curvature), 'PlotFcn', @(options, state, flag)modifiedGaplotBestF(options, state, flag));
[x, fval] = ga(fun, 3+2*Segment_num, A, b, Aeq, beq, lb, ub, nonlcon, intcon, options);
end

function objective = GAObj(x, Length, Curvature)
obj = 0;
for i = 1:length(Length)
    load NewANN.mat
    unit_angle_ANN = 2*sim(net,[x(1), x(2), x(3+2*(i-1)+1), x(3+2*(i-1)+2), x(3)]');
    new_l = (x(3+2*(i-1)+2)/unit_angle_ANN + x(1) + x(2)/2)*unit_angle_ANN;
    bellow_num = round(Length(i)/new_l);
    obj = obj + abs(Length(i) - bellow_num*new_l)/Length(i) + abs(Curvature(i) - unit_angle_ANN/new_l)/Curvature(i);
end
objective = obj/length(Length);

end

function [state,options,optchanged] = GAoutputfun(options,state,flag, Length, Curvature)

optchanged = false;
switch flag
    case 'init'

    case 'iter'

        ibest = state.Best(end);
        ibest = find(state.Score == ibest,1,'last');
        bestx = state.Population(ibest,:);
        bestf = GAObj(bestx, Length, Curvature);
        if bestf < 0.0255
            state.StopFlag = 'y';
            disp(bestx);
            disp(['Best function value is ', bestf]);
        end
end
end

function state = modifiedGaplotBestF(options, state, flag)
    % Call the original gaplotbestf function
    state = gaplotbestf(options, state, flag);
    
    % Adjust the window position after the plot is generated
%     if strcmp(flag, 'done')
        % Retrieve the handle of the GA plot figure
        hFig = gcf;
        
        % Get screen size
        screenSize = get(groot, 'ScreenSize');
        screenWidth = screenSize(3);
        screenHeight = screenSize(4);
        % Set the figure position
        figWidth = 560;  % Set the width of the figure
        figHeight = 420; % Set the height of the figure
        
        % Set the window position using the specified coordinates
        set(hFig, 'Position',[(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight]);  % Adjust width and height as desired
%     end
end




