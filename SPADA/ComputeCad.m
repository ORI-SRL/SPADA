function out = ComputeCad(filepath, filename, PCC_result, uniqueArc, repeatIndex, optVariable)
% create a geometry of the designed actuator in COMSOL and export a STL file
% 
% Input arguments:
% PCC_result: 
%   structure array
%   struct('type', {1}, 'center', {[0,0,0]}, 'axis1', {[0,0,0]}, 'axis2', {[0,0,0]}, 'radius', {0}, 'angle', {0}, 'length', {0}, 'rotation', {0});
% uniqueArc: 
%   numberOfUniqueArc-by-2 matrix of the arclength and curvature of the unique arcs
% repeatIndex:
%   list of arcs with its index of the unique arcs
% optVariable:
%   optimal design parameters for the unique arcs


r = optVariable(1); % inner radius in mm
t = optVariable(2); % wall thickness in mm
R = []; % outer radius in mm
l = []; % module length in mm
N = []; % number of bellow in each segment

for i = 1:(length(optVariable)-3)/2
    R = [R optVariable(3+(i-1)*2+1)];
    l = [l optVariable(3+(i-1)*2+2)]/4;
    l = l*4;
    N = [N round(uniqueArc(i,1)/l(i))];
end

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath(filepath);

model.label('CAD.mph');


model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').geomRep('comsol');

actuatorLength = 0;
rotationAngle = 0;
m = 1; % line_index
j = 1; % arc_index
k = 1; % wp_index

for i = 1:size(PCC_result,2) % segment_index
    
    if PCC_result(i).type == 1 
        % line segment
        wp = "wp" + num2str(k);
        rev = "rev" + num2str(k);
        k = k + 1; % wp_index
        arr = "arr" + num2str(i);
        lineLabel = "line" + num2str(m);
        m = m + 1; % line index
        
        
        R_line = max(R);
        N_line = ceil(PCC_result(i).length/4/(R_line-r));
        l_line = PCC_result(i).length/N_line;
        
        model.component('comp1').geom('geom1').create(wp, 'WorkPlane');
        model.component('comp1').geom('geom1').feature(wp).label(lineLabel);
        model.component('comp1').geom('geom1').feature(wp).set('unite', true);
        model.component('comp1').geom('geom1').feature(wp).geom.create('ca1', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca1').set('center', [0 r+l_line/4]); % [0 r+l_line/4]
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca1').set('r', l_line/4-t/2); % l_line/4-t/2
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca1').set('angle1', 270);
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca1').set('angle2', 0);
        model.component('comp1').geom('geom1').feature(wp).geom.create('ls1', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls1').set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls1').set('coord2', [l_line/4-t/2 R_line*2-r-l_line/4]); % [l_line/4-t/2 R_line*2-r-l_line/4]
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls1').selection('vertex1').set('ca1(1)', 2);
        model.component('comp1').geom('geom1').feature(wp).geom.create('ca2', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca2').set('center', [l_line/2 R_line*2-r-l_line/4]); %[l_line/2 R_line*2-r-l_line/4]
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca2').set('r', l_line/4+t/2); % l_line/4+t/2
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca2').set('angle1', 90);
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ca2').set('angle2', 180);
        model.component('comp1').geom('geom1').feature(wp).geom.create('ls3', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls3').set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls3').set('coord2', [0 t]); % [0 t]
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls3').selection('vertex1').set('ca1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp).geom.create('ccur1', 'ConvertToCurve');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ccur1').selection('input').set({'ca1' 'ca2' 'ls1' 'ls3'});
        model.component('comp1').geom('geom1').feature(wp).geom.create('mir1', 'Mirror');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('mir1').set('keep', true);
        model.component('comp1').geom('geom1').feature(wp).geom.feature('mir1').set('pos', [l_line/2 0]); %[l_line/2 0]
        model.component('comp1').geom('geom1').feature(wp).geom.feature('mir1').selection('input').set({'ccur1'});
        model.component('comp1').geom('geom1').feature(wp).geom.create('ls4', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls4').selection('vertex1').set('ccur1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp).geom.feature('ls4').selection('vertex2').set('mir1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp).geom.create('csol1', 'ConvertToSolid');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('csol1').selection('input').set({'ccur1' 'ls4' 'mir1'});
        model.component('comp1').geom('geom1').feature(wp).geom.create('mov1', 'Move');
        model.component('comp1').geom('geom1').feature(wp).geom.feature('mov1').set('displx', actuatorLength); % previous length
        model.component('comp1').geom('geom1').feature(wp).geom.feature('mov1').selection('input').set({'csol1'});
        model.component('comp1').geom('geom1').create(rev, 'Revolve');
        model.component('comp1').geom('geom1').feature(rev).set('angtype', 'full');
        model.component('comp1').geom('geom1').feature(rev).set('axis', [1 0]);
        model.component('comp1').geom('geom1').feature(rev).selection('input').set(wp);
        model.component('comp1').geom('geom1').create(arr, 'Array');
        model.component('comp1').geom('geom1').feature(arr).set('fullsize', [N_line 1 1]);
        model.component('comp1').geom('geom1').feature(arr).set('displ', [l_line 0 0]);
        model.component('comp1').geom('geom1').feature(arr).selection('input').set(rev);
        actuatorLength = actuatorLength + PCC_result(i).length;
        
    else
        % arc segment
        wp1 = "wp" + num2str(k);
        wp2 = "wp" + num2str(k+1);
        rev1 = "rev" + num2str(k);
        rev2 = "rev" + num2str(k+1);
        arr = "arr" + num2str(i);
        k = k + 2; % wp_index
        bellowLabel = "arc" + num2str(j) + "-bellow";
        constraintLabel = "arc" + num2str(j) + "-constraint";
        rotationAngle = rotationAngle + PCC_result(i).rotation*180/pi;
        
        model.component('comp1').geom('geom1').create(wp1, 'WorkPlane');
        model.component('comp1').geom('geom1').feature(wp1).label(bellowLabel);
        model.component('comp1').geom('geom1').feature(wp1).set('unite', true);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ca1', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca1').set('center', [0 r+l(repeatIndex(j))/4]); % [0 r+l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca1').set('r', l(repeatIndex(j))/4 + t/2); % l(repeatIndex(j))/4 + t/2
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca1').set('angle1', 270); 
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca1').set('angle2', 0);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ls1', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls1').set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls1').set('coord2', [l(repeatIndex(j))/4+t/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]); % [l(repeatIndex(j))/4+t/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls1').selection('vertex1').set('ca1(1)', 2);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ca2', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca2').set('center', [l(repeatIndex(j))/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]); % [l(repeatIndex(j))/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca2').set('r', l(repeatIndex(j))/4-t/2); % l(repeatIndex(j))/4-t/2
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca2').set('angle1', 90);
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca2').set('angle2', 180);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ca3', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca3').set('center', [0 r+l(repeatIndex(j))/4]); % r+l(repeatIndex(j))/4
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca3').set('r', l(repeatIndex(j))/4-t/2); % l(repeatIndex(j))/4-t/2
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca3').set('angle1', 270);
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca3').set('angle2', 0);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ls2', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls2').set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls2').set('coord2', [l(repeatIndex(j))/4-t/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]); % [l(repeatIndex(j))/4-t/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls2').selection('vertex1').set('ca3(1)', 2);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ca4', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca4').set('center', [l(repeatIndex(j))/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]); % [l(repeatIndex(j))/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca4').set('r', l(repeatIndex(j))/4+t/2); % l(repeatIndex(j))/4+t/2
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca4').set('angle1', 90);
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ca4').set('angle2', 180);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ls3', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls3').selection('vertex1').set('ca1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ls3').selection('vertex2').set('ca3(1)', 1);
        model.component('comp1').geom('geom1').feature(wp1).geom.create('ccur1', 'ConvertToCurve');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('ccur1').selection('input').set({'ca1' 'ca2' 'ca3' 'ca4' 'ls1' 'ls2' 'ls3'});
        model.component('comp1').geom('geom1').feature(wp1).geom.create('mir1', 'Mirror');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('mir1').set('keep', true);
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('mir1').set('pos', [l(repeatIndex(j))/2 0]); %[l(repeatIndex(j))/2 0]
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('mir1').selection('input').set({'ccur1'});
        model.component('comp1').geom('geom1').feature(wp1).geom.create('csol1', 'ConvertToSolid');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('csol1').selection('input').set({'ccur1' 'mir1'});
        model.component('comp1').geom('geom1').feature(wp1).geom.create('mov1', 'Move');
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('mov1').set('displx', actuatorLength); % previous length
        model.component('comp1').geom('geom1').feature(wp1).geom.feature('mov1').selection('input').set({'csol1'});
        model.component('comp1').geom('geom1').create(rev1, 'Revolve');
        model.component('comp1').geom('geom1').feature(rev1).set('angtype', 'full');
        model.component('comp1').geom('geom1').feature(rev1).set('axis', [1 0]);
        model.component('comp1').geom('geom1').feature(rev1).selection('input').set(wp1);
        
        % create fan-shaped constraint
        model.component('comp1').geom('geom1').create(wp2, 'WorkPlane');
        model.component('comp1').geom('geom1').feature(wp2).label(constraintLabel);
        model.component('comp1').geom('geom1').feature(wp2).set('unite', true);
        model.component('comp1').geom('geom1').feature(wp2).geom.create('ca1', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca1').set('center', [0 r+l(repeatIndex(j))/4]); % [0 r+l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca1').set('r', l(repeatIndex(j))/4+r/2); % l(repeatIndex(j))/4+r/2
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca1').set('angle1', 270);
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca1').set('angle2', 0);
        model.component('comp1').geom('geom1').feature(wp2).geom.create('ls1', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls1').set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls1').set('coord2', [l(repeatIndex(j))/4+t/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]); % [l(repeatIndex(j))/4+t/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls1').selection('vertex1').set('ca1(1)', 2);
        model.component('comp1').geom('geom1').feature(wp2).geom.create('ca2', 'CircularArc');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca2').set('center', [l(repeatIndex(j))/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]); % [l(repeatIndex(j))/2 R(repeatIndex(j))*2-r-l(repeatIndex(j))/4]
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca2').set('r', l(repeatIndex(j))/4-t/2); % l(repeatIndex(j))/4-t/2
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca2').set('angle1', 90);
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ca2').set('angle2', 180);
        model.component('comp1').geom('geom1').feature(wp2).geom.create('ls2', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls2').set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls2').set('coord2', [0 t]); % [0 t]
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls2').selection('vertex1').set('ca1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp2).geom.create('ccur1', 'ConvertToCurve');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ccur1').selection('input').set({'ca1' 'ca2' 'ls1' 'ls2'});
        model.component('comp1').geom('geom1').feature(wp2).geom.create('mir1', 'Mirror');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('mir1').set('keep', true);
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('mir1').set('pos', [l(repeatIndex(j))/2 0]); % [l(repeatIndex(j))/2 0]
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('mir1').selection('input').set({'ccur1'});
        model.component('comp1').geom('geom1').feature(wp2).geom.create('ls3', 'LineSegment');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls3').selection('vertex1').set('ccur1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('ls3').selection('vertex2').set('mir1(1)', 1);
        model.component('comp1').geom('geom1').feature(wp2).geom.create('csol1', 'ConvertToSolid');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('csol1').selection('input').set({'ccur1' 'ls3' 'mir1'});
        model.component('comp1').geom('geom1').feature(wp2).geom.create('mov1', 'Move');
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('mov1').set('displx', actuatorLength); % previous length
        model.component('comp1').geom('geom1').feature(wp2).geom.feature('mov1').selection('input').set({'csol1'});
        model.component('comp1').geom('geom1').create(rev2, 'Revolve');
        model.component('comp1').geom('geom1').feature(rev2).set('angle1', -45+rotationAngle);
        model.component('comp1').geom('geom1').feature(rev2).set('angle2', 45+rotationAngle);
        model.component('comp1').geom('geom1').feature(rev2).set('axis', [1 0]);
        model.component('comp1').geom('geom1').feature(rev2).selection('input').set(wp2);
        model.component('comp1').geom('geom1').create(arr, 'Array');
        model.component('comp1').geom('geom1').feature(arr).set('fullsize', [N(repeatIndex(j)) 1 1]); % [N(repeatIndex(j)) 1 1]
        model.component('comp1').geom('geom1').feature(arr).set('displ', [l(repeatIndex(j)) 0 0]); % [l(repeatIndex(j)) 0 0]
        model.component('comp1').geom('geom1').feature(arr).selection('input').set({char(rev1), char(rev2)});

        actuatorLength = actuatorLength + N(repeatIndex(j))*l(repeatIndex(j));
        j = j + 1; % arc_index
    end
end
% add top and bottom parts
wpt = "wp" + num2str(k);
wpb = "wp" + num2str(k+1);
model.component('comp1').geom('geom1').create(wpt, 'WorkPlane');
model.component('comp1').geom('geom1').feature(wpt).label('top');
model.component('comp1').geom('geom1').feature(wpt).set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature(wpt).set('unite', true);
model.component('comp1').geom('geom1').feature(wpt).geom.create('c1', 'Circle');
model.component('comp1').geom('geom1').feature(wpt).geom.feature('c1').set('r', t);
model.component('comp1').geom('geom1').feature(wpt).geom.create('c2', 'Circle');
model.component('comp1').geom('geom1').feature(wpt).geom.feature('c2').set('r', r+t/2);
model.component('comp1').geom('geom1').feature(wpt).geom.create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature(wpt).geom.feature('dif1').selection('input').set({'c2'});
model.component('comp1').geom('geom1').feature(wpt).geom.feature('dif1').selection('input2').set({'c1'});
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', '5', 0);
model.component('comp1').geom('geom1').feature('ext1').set('reverse', true);
model.component('comp1').geom('geom1').feature('ext1').selection('input').set(wpt);
model.component('comp1').geom('geom1').create(wpb, 'WorkPlane');
model.component('comp1').geom('geom1').feature(wpb).set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature(wpb).set('quickx', actuatorLength);
model.component('comp1').geom('geom1').feature(wpb).set('unite', true);
model.component('comp1').geom('geom1').feature(wpb).geom.create('c1', 'Circle');
model.component('comp1').geom('geom1').feature(wpb).geom.feature('c1').set('r', t);
model.component('comp1').geom('geom1').feature(wpb).geom.create('c2', 'Circle');
model.component('comp1').geom('geom1').feature(wpb).geom.feature('c2').set('r', r+t/2);
model.component('comp1').geom('geom1').feature(wpb).geom.create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature(wpb).geom.feature('dif1').selection('input').set({'c2'});
model.component('comp1').geom('geom1').feature(wpb).geom.feature('dif1').selection('input2').set({'c1'});
model.component('comp1').geom('geom1').create('ext2', 'Extrude');
model.component('comp1').geom('geom1').feature('ext2').setIndex('distance', '5', 0);
model.component('comp1').geom('geom1').feature('ext2').selection('input').set(wpb);
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').geom('geom1').export.setSTLFormat('text');
model.component('comp1').geom('geom1').export(filepath + "\" + filename + ".stl");

out = model;
end