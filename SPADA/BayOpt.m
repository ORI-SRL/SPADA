function output = BayOpt(ArcSegment, BD)
% Upper Level Optimisation
global OptVariable OptResult
OptVariable = [];
OptResult = [];

Segment_num = size(ArcSegment,1);
Segment_length = ArcSegment(:,1); %mm
Segment_curvature = ArcSegment(:,2); % 1/mm

[max_curvature, ~] = max(Segment_curvature);
[min_length, ~] = min(Segment_length);


ro_lb = BD(1,1); 
ro_ub = BD(1,2);
P_lb = BD(2,1);
P_ub = BD(2,2);

optimVar_inner_radius = optimizableVariable('ri',[0, min(1/max_curvature, ro_ub)]);

% optimVar_thickness = optimizableVariable('t',[0,min(1/max_curvature, Segment_length(index_max_curvature))/4]);
optimVar_thickness = optimizableVariable('t',[0,min_length/4]);


pressure = optimizableVariable('Pf',[max(0,P_lb),min(10000,P_ub)]);

optimVars = [optimVar_inner_radius optimVar_thickness pressure];


fun = @(x) BayUpper(x, Segment_num, Segment_length, Segment_curvature, [ro_lb, ro_ub]);

results = bayesopt(fun, optimVars, 'XConstraintFcn', @xconstraint, 'OutputFcn',@outputfun, 'MaxObjectiveEvaluations',1000, 'NumSeedPoints', 1, 'InitialX', table(5.5, 1.3, 3000));
% results = bayesopt(fun, optimVars, 'XConstraintFcn', @xconstraint, 'OutputFcn',@outputfun, 'MaxObjectiveEvaluations',1000);
h = gcf;
h.CurrentAxes.Title.String = 'Bayesian Optimisation';
if isnan(results.MinObjective)
    output = [];
else
    output = OptVariable(end,:);
end

end

function tf = xconstraint(X)
% geometric boundary constraints

tf1 = X.t >= X.ri/4; % t >= ri/4
tf2 = X.t < X.ri/2; % t < ri/2

tf = tf1 & tf2;
end

function stop = outputfun(results,state)
stop = false;
switch state
    case 'initial'

    case 'iteration'
        if results.MinObjective < 0.02
            stop = true;
        end

end
end



% Upper Level Objective and Constraint
function objective = BayUpper(x, Segment_num, Segment_length, Segment_curvature, ro_BD)
global OptVariable OptResult
disp('---------------------');
disp([x.ri, x.t, x.Pf]);

ro_lb = ro_BD(1);
ro_ub = ro_BD(2);

CurrentOptVariable = [x.ri, x.t, x.Pf];
obj_value = zeros(1,Segment_num);

for i = 1:Segment_num
    disp('Segment Number:');
    disp(i);
    
    rm_lb = max(x.ri + x.t, (x.ri + ro_lb)/2);
    rm_ub = min([(1/Segment_curvature(i)+x.ri)/2, x.ri*2, (x.ri + ro_ub)/2]);
    l_lb = 4*x.t;
%     l_ub = min([1/Segment_curvature(i), Segment_length(i), 4*x.ri]);
    l_ub = min([Segment_length(i), 4*x.ri]);
    
    if 4*(rm_ub - x.ri) >= l_lb
        feasible_point = true;
    else
        feasible_point = false;
    end
    
    if (rm_lb < rm_ub)&&(l_lb < l_ub)&&(feasible_point)
        optimVar_aver_radius = optimizableVariable('rm',[rm_lb, rm_ub]);
        optimVar_unit_length = optimizableVariable('l',[l_lb, l_ub]);
        
        fun = @(y) BayLower(y, x.ri, x.t, x.Pf, Segment_length(i), Segment_curvature(i));
        lowerconst = @(y) yconstraint(y, x.ri);
        
        results = bayesopt(fun, [optimVar_aver_radius optimVar_unit_length], 'XConstraintFcn', lowerconst, 'OutputFcn', @outputfun, 'MaxObjectiveEvaluations',15, 'PlotFcn', []);
        obj_value(i) = results.MinObjective;
%         constraint(i) = results.MinObjective - 0.05^2*(Segment_length(i)*Segment_curvature(i))^2*Segment_length(i)^2;
%         constraint(i) = results.MinObjective - 0.02^2;
        
        CurrentOptVariable = [CurrentOptVariable results.XAtMinObjective.Variables];
        
    else
        obj_value(i) = 800;
%         constraint(i) = 800;
        CurrentOptVariable = [CurrentOptVariable [-1 -1]];
    end
end

% objective  = sum(sqrt(obj_value/((Segment_length(i)*Segment_curvature(i))^2*Segment_length(i)^2)));
objective = sum(obj_value)/Segment_num;
OptVariable = [OptVariable; CurrentOptVariable];
OptResult = [OptResult objective];


end

function tf = yconstraint(Y, ri)
% geometric boundary constraints
tf1 = Y.l <= 4*(Y.rm-ri); % l <= min(4*(rm-ri), L)
tf = tf1;
end

function objective = BayLower(y, Inner_radius, Thickness, Pressure, Length, Curvature)

load NewANN.mat
unit_angle_ANN = 2*sim(net,[Inner_radius, Thickness, y.rm, y.l, Pressure]');

bellow_num = round(Length/y.l);

% objective = (Length*Curvature)^2*(Length - bellow_num*y.l)^2 + (Length)^2*(Length*Curvature - bellow_num*unit_angle_ANN)^2;
% objective = sqrt(((Length*Curvature)^2*(Length - bellow_num*y.l)^2 + (Length)^2*(Length*Curvature - bellow_num*unit_angle_ANN)^2)/(Length^2*Curvature)^2);
objective = sqrt(((Curvature)^2*(Length - bellow_num*y.l)^2 + Length^2*(Curvature - unit_angle_ANN/y.l)^2)/(Length*Curvature)^2);


end