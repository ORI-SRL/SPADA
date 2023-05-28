function [sample_no, material_name, material_value, variable_name, variable_value, result_name, result_value] =  GenerateDataset(E,nu,rho)
% only one material is considered for the dataset

% Material properties: Young's Modulus, Poisson's Ratio, Density
material_name = {'Youngs Modulus'; 'Poissons Ratio'; 'Density'};
material_value = zeros(3,1);

% Design Variable: ri, t, l, rm
variable_name = {'Inner Radius'; 'Wall Thickness'; 'Average Radius'; 'Unit Length'};
variable_value = zeros(4, 1);

% Results: pressure, bending angle, principal strain, elastic strain energy
result_name = {'Pressure'; 'Angular Deflection'; 'Principal Strain'; 'Elastic Strain Energy'};
result_value = zeros(4, 21);

N_sample = 1;

pause(2);
S = [];
for ri = 2:2:10
    for t = ri/4:(ri/3-ri/4)/2:ri/3
        for rm = ri+t:(2*ri - (ri+t))/4:2*ri
            for l = max(4*t,ri):(4*(rm-ri)-max(4*t,ri))/4:4*(rm-ri)
                disp('---------------------------------');
                disp('N_sample: ' + string(N_sample));
                
                material_value(:,:,N_sample) = [E; nu; rho*1000]; 
                variable_value(:,:,N_sample) = [ri; t; rm; l];
                output = DataSetSimulation(material_value(:,:,N_sample),variable_value(:,:,N_sample), 10000);
%                 output = zeros(4, 21);
                if output ~= 800
                    if length(output) < 21
                        NaNset = NaN(4, 21 - length(output));
                        output = [output NaNset];
                    end
                    result_value(:,:,N_sample) = output;
                    N_sample = N_sample + 1;
                else
                    continue;   
                end
            end
        end
    end
end
N_sample = N_sample - 1;
sample_no = 1:N_sample;
sample_no = transpose(sample_no);
end

function output = DataSetSimulation(material_value, variable_value, final_pressure)
disp('-----------------------------------');
tic
import com.comsol.model.*
import com.comsol.model.util.*
flank = variable_value(3)*2 - 2*variable_value(1) - variable_value(4)/2;

if flank >= 0
    % geometry checked
    try

        Ps = 500; % pressure step [Pa]

        model = mphopen('Half_Unit_DataSet.mph', 'Model');
        % set material properties
        model.param.set('Es', string(material_value(1))+'[kPa]', 'soft material stiffness');
        model.param.set('miu', string(material_value(2)), "Possion's Ratio");
        model.param.set('density', string(material_value(3))+'[kg/m^3]', 'density');

        % set design variables

        model.param.set('ri', string(variable_value(1))+'[mm]', 'average inner radius');
        model.param.set('t', string(variable_value(2))+'[mm]', 'thickness');
        model.param.set('rm', string(variable_value(3))+'[mm]', 'average average radius'); % check rank
        model.param.set('l', string(variable_value(4))+'[mm]', 'single unit length'); % check rank
        
        % set pressure
        model.param.set('Pf', string(final_pressure)+'[Pa]', 'final pressure'); % check for unit


        model.geom('geom1').runAll;
        model.geom('geom1').run;

        model.sol('sol1').runAll;

        model.result.table.create('tbl2', 'Table');
        model.result.numerical('int1').set('table', 'tbl2');
        model.result.numerical('int1').setResult;

        tab1 = mphtable(model,'tbl1');
        P = tab1.data(:,1); % Pa
        ep1 = tab1.data(:,2); % 1
        angle = tab1.data(:,3); % rad
        tab2 = mphtable(model,'tbl2');
        energy = tab2.data(:,2);

        results = [transpose(P); transpose(angle); transpose(ep1); transpose(energy)];
        disp(results);

        model.result.table.remove('tbl1');
        model.result.table.remove('tbl2');

        ModelUtil.remove('Model');

        output = results;

    catch msg
        disp(msg);
        output = 800;
    end
else
    output = 800;
    disp('Negative flank length.');
end
toc;
end
