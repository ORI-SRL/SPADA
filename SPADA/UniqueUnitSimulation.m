function output = UniqueUnitSimulation(uniqueUnit, pressure, material)
% Simulate half modules' deformation using COMSOL
%
% Input arguments:
% uniqueUnit:
%    n-by-4 matrix of the unique modules' design parameters: r, t, R, l
% pressure:
%    actuation pressure [Pa]
% material:
%    1-by-3 array of material properties: Young's Modulus [kPa], Possion's Ratio, Density [kg/m^3]
%
% Output arguments:
% output:
%    list of the angular deflection of unique modules or error messages

if isempty(uniqueUnit)
    output = [];
else
    uniqueUnitResult = [];
    for i = 1:size(uniqueUnit,1)

        flank = uniqueUnit(i,3)*2 - 2*uniqueUnit(i,1) - uniqueUnit(i,4)/2;
        
        if flank >= 0
            % geometry checked
            try

                Ps = 500; % pressure step [Pa]

                model = mphopen('Half_Unit_Simulation.mph', 'Model');
                % set material properties
                model.param.set('Es', string(material(1))+'[kPa]', 'soft material stiffness');
                model.param.set('miu', string(material(2)), "Possion's Ratio");
                model.param.set('density', string(material(3))+'[kg/m^3]', 'density');

                % set design variables

                model.param.set('ri', string(uniqueUnit(i,1))+'[mm]', 'average inner radius');
                model.param.set('t', string(uniqueUnit(i,2))+'[mm]', 'thickness');
                model.param.set('rm', string(uniqueUnit(i,3))+'[mm]', 'average average radius'); 
                model.param.set('l', string(uniqueUnit(i,4))+'[mm]', 'single unit length');

                % set pressure
                model.param.set('Pf', string(pressure)+'[Pa]', 'final pressure'); % check for unit


                model.geom('geom1').runAll;
                model.geom('geom1').run;
                
                if (round(pressure/Ps)*Ps == pressure) 
                    model.sol('sol1').feature('v1').set('clist', {'range(0,Ps,Pf)'});
                    model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(0,Ps,Pf)'});
                    end_step = pressure/Ps + 1;
                    
                else
                    model.sol('sol1').feature('v1').set('clist', {'range(0,Ps,Pf) Pf'});
                    model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(0,Ps,Pf) Pf'});
                    end_step = floor(pressure/Ps) + 2;
                end

                model.sol('sol1').runAll;
                
                model.result('pg1').set('solnum', end_step);
                model.result('pg1').run;

                model.result('pg3').set('solnum', end_step);
                model.result('pg3').run;

    
                tab1 = mphtable(model,'tbl1');
                angle = 2*tab1.data(end,end); % rad
                uniqueUnitResult = [uniqueUnitResult angle];
                model.result.table.remove('tbl1');
                
                bellow_geometry = sprintf('( ri = %s, t = %s, rm = %s, l = %s)', string(uniqueUnit(i,1)), string(uniqueUnit(i,2)), string(uniqueUnit(i,3)), string(uniqueUnit(i,4)));
                
                % Get screen size
                screenSize = get(groot, 'ScreenSize');
                screenWidth = screenSize(3);
                screenHeight = screenSize(4);
                % Set the figure position
                figWidth = 560;  % Set the width of the figure
                figHeight = 420; % Set the height of the figure
                
                stress = figure("Name", "Stress Plot of Bellow Module "+ bellow_geometry, "NumberTitle", "off", 'visible', 'on');
                stress.Position =  [(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight];

                hold on
                mphplot(model,'pg1','rangenum',1) 
                view(-30,30);
                xlabel('X [mm]');
                ylabel('Y [mm]');
                zlabel('Z [mm]');
                hold off

                strain = figure("Name", "Strain Plot of Bellow Module "+ bellow_geometry, "NumberTitle", "off", 'visible', 'on');
                strain.Position =  [(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight];
                hold on
                mphplot(model,'pg3','rangenum',1)
                view(-30,30);
                xlabel('X [mm]');
                ylabel('Y [mm]');
                zlabel('Z [mm]');
                hold off

                
            catch msg
                output = [uniqueUnitResult -1];
                break;
            end
        else
            output = [uniqueUnitResult -2];
            break;
        end
    end
%     ModelUtil.remove('Model');
    output = uniqueUnitResult;
end