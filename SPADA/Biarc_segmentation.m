function [PCC_figure, PCC_result] = Biarc_segmentation(bspline, tangential_vectors, accuracy)
% Biarc segmentation
%
% Input arguments
% bspline:
%    3-by-m matrix of the projected points of the original data on the approximated b-spline
% tangential_vectors:
%    3-by-m matrix of the normalized tangential vectors of the projected points on the b-spline
% accuracy: 
%    the maximum error allowed between the segmented arc and b-spline curve
% 
% Output arguments:
% allArc:
%    segmented constant curvature segments
% arcIndex:
%    list of the points that separate biarcs

tol = 10e-2;
M = bspline;
T = tangential_vectors;

iIndex = 1;
eIndex = 3;
j = 1;
allArc = [];
    
while eIndex <= size(M,2)
    % D =[D ga(@(d1)biarcError (M, T, iIndex, eIndex, d1, tol),1)];
    error = biarcError(M, T, iIndex, eIndex, tol);
    if error <= accuracy
        eIndex = eIndex + 1;
        continue;
    else
        eIndex = eIndex - 1;
        Arc = biarc(M, T, iIndex, eIndex, tol);
        allArc = [allArc Arc];
        error = biarcError(M, T, iIndex, eIndex, tol);
        
        j = j+1;
        iIndex = eIndex;
        eIndex = iIndex+2;
    end
end
eIndex = eIndex - 1;
allArc = [allArc biarc(M, T, iIndex, eIndex, tol)];

% do the classification of the constant curvature segments

segment = struct('type', {1}, 'center', {[0,0,0]}, 'axis1', {[0,0,0]}, 'axis2', {[0,0,0]}, 'radius', {0}, 'angle', {0}, 'length', {0}, 'rotation', {0});

PCC_result = [];

for i = 1: size(allArc,2)
    if isempty(PCC_result)
        PCC_result = [PCC_result BiarcToSegment(allArc(i), segment, i)];
    else
        segment = BiarcToSegment(allArc(i), segment, i);
        if segment.type == 1 && PCC_result(end).type == 1
            % two lines can be combined
            comb_segment = struct('type', {1}, 'center', {[0,0,0]}, 'axis1', {[0,0,0]}, 'axis2', {[0,0,0]}, 'radius', {0}, 'angle', {0}, 'length', {0}, 'rotation', {0});
            comb_segment.center = (PCC_result(end).center + PCC_result(end).axis1 + segment.center - segment.axis1)/2;
            comb_segment.axis1 = PCC_result(end).center + PCC_result(end).axis1 - comb_segment.center;
            comb_segment.axis2 = [0,0,0];
            comb_segment.radius = 0;
            comb_segment.angle = 0;
            comb_segment.length = PCC_result(end).length + segment.length;
            comb_segment.rotation = 0;
            PCC_result(end) = comb_segment; % replace the segment
            
        elseif segment.type == 2 && PCC_result(end).type == 2
            n1 = cross(PCC_result(end).axis1, PCC_result(end).axis2);
            n1 = n1/norm(n1);
            n2 = cross(segment.axis1, segment.axis2);
            n2 = n2/norm(n2);
            
            if (norm(cross(n1,n2)) <= tol) && (dot(n1,n2) > 0) && (abs(PCC_result(end).radius - segment.radius)/max(PCC_result(end).radius, segment.radius) <= tol)
                % two arcs can be combined
                comb_segment = struct('type', {2}, 'center', {[0,0,0]}, 'axis1', {[0,0,0]}, 'axis2', {[0,0,0]}, 'radius', {0}, 'angle', {0}, 'length', {0}, 'rotation', {0});
                comb_segment.center = PCC_result(end).center;
                comb_segment.axis1 = PCC_result(end).axis1;
                comb_segment.axis2 = PCC_result(end).axis2;
                comb_segment.radius = PCC_result(end).radius;
                comb_segment.angle = PCC_result(end).angle + segment.angle;
                comb_segment.length = PCC_result(end).length + segment.length;
                comb_segment.rotation = PCC_result(end).rotation;
                PCC_result(end) = comb_segment; % replace the segment
            else
                % two arcs are separate
                v1 = -(PCC_result(end).axis1*cos(PCC_result(end).angle) + PCC_result(end).axis2*sin(PCC_result(end).angle));
                v2 = -segment.axis1;
                cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
                sin_angle = norm(cross(v1, v2)) / (norm(v1) * norm(v2));
                segment.rotation = atan2(norm(cross(v1, v2)), dot(v1, v2));% positive (counter-clockwise) or negative (clockwise).
                PCC_result = [PCC_result segment];
            end
        else
            PCC_result = [PCC_result BiarcToSegment(allArc(i), segment, i)];
        end
    end
        
end
    function segment = BiarcToSegment(arc, segment, i)
        if arc.angle <= 0.11 || arc.arcLen <= 1 || arc.angle/arc.arcLen <= 0.011
            % line segment
            segment.type = 1; 
            segment.axis2 = [0,0,0];
            segment.radius = 0;
            segment.angle = 0;
            segment.rotation = 0;
            
            if arc.radius == 0
                segment.center = arc.center;
                if rem(i,2) == 1
                    segment.axis1 = arc.axis1;
                else
                    segment.axis1 = -arc.axis1;
                end
            else
                startPoint = arc.center + arc.axis1;
                endPoint = arc.center + arc.axis1*cos(arc.angle) + arc.axis2*sin(arc.angle);
                segment.center = (startPoint + endPoint)/2;
                if rem(i,2) == 1
                    segment.axis1 = arc.center + arc.axis1 - segment.center;
                else
                    segment.axis1 = arc.center + arc.axis1*cos(arc.angle) + arc.axis2*sin(arc.angle) - segment.center;
                end
            end
            segment.length = 2*norm(segment.axis1);
                
        else
            % arc segment
            segment.type = 2;
            segment.center = arc.center;
            segment.radius = arc.radius;
            segment.angle = arc.angle;
            segment.length = arc.arcLen;
            segment.rotation = 0; % will be determined later
            if rem(i,2) == 1
                segment.axis1 = arc.axis1;
                segment.axis2 = arc.axis2;
            else
                segment.axis1 = arc.midPoint - arc.center;
                segment.axis2 = cross(segment.axis1, cross(arc.axis1, arc.axis2)/(arc.radius)^2);
            end
        end
    end

% for i = 1:2:size(allArc,2)-1 
%     if (allArc(i).angle <= 0.05) && (allArc(i+1).angle <= 0.05)
%         % line segment
%         segment.type = 1; 
%         segment.center = allArc(i).midPoint;
%         segment.axis1 = allArc(i).center + allArc(i).axis1 - segment.center;
%         segment.axis2 = [0,0,0];
%         segment.radius = 0;
%         segment.angle = 0;
%         segment.length = 2*norm(segment.axis1);
%         segment.rotation = 0;
%        
%         PCC_result = [PCC_result segment];
%     else
%         % arc segment
%         segment.type = 2; 
%         n1 = cross(allArc(i).axis1, allArc(i).axis2);
%         n1 = n1/norm(n1);
%         n2 = cross(allArc(i+1).axis2, allArc(i+1).axis1);
%         n2 = n2/norm(n2);
%         if (norm(cross(n1,n2)) <= tol) && (dot(n1,n2) > 0) && (abs(allArc(i).radius - allArc(i+1).radius)/max(allArc(i).radius,allArc(i+1).radius) <= tol)
%             % two arcs can be combined
%             segment.center = allArc(i).center;
%             segment.axis1 = allArc(i).axis1;
%             segment.axis2 = allArc(i).axis2;
%             segment.radius = allArc(i).radius;
%             segment.angle = allArc(i).angle + allArc(i+1).angle;
%             segment.length = allArc(i).arcLen + allArc(i+1).arcLen;
%             
%             if isempty(PCC_result) || PCC_result(end).type == 1
%                 segment.rotation = 0;
%             else
%                 v1 = -(PCC_result(end).axis1*cos(PCC_result(end).angle) + PCC_result(end).axis2*sin(PCC_result(end).angle));
%                 v2 = -segment.axis1;
%                 cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
%                 sin_angle = norm(cross(v1, v2)) / (norm(v1) * norm(v2));
%                 segment.rotation = atan2(norm(cross(v1, v2)), dot(v1, v2));% positive (counter-clockwise) or negative (clockwise).
%             end
%             PCC_result = [PCC_result segment];
%         else
%             % two separate arcs
%             % i-th 
%             segment.center = allArc(i).center;
%             segment.axis1 = allArc(i).axis1;
%             segment.axis2 = allArc(i).axis2;
%             segment.radius = allArc(i).radius;
%             segment.angle = allArc(i).angle;
%             segment.length = allArc(i).arcLen;
%             
%             if isempty(PCC_result) || PCC_result(end).type == 1
%                 segment.rotation = 0;
%             else
%                 v1 = -(PCC_result(end).axis1*cos(PCC_result(end).angle) + PCC_result(end).axis2*sin(PCC_result(end).angle));
%                 v2 = -segment.axis1;
%                 cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
%                 sin_angle = norm(cross(v1, v2)) / (norm(v1) * norm(v2));
%                 segment.rotation = atan2(norm(cross(v1, v2)), dot(v1, v2));% positive (counter-clockwise) or negative (clockwise).
%             end
%             PCC_result = [PCC_result segment];
%             
%             % i+1-th
%             segment.center = allArc(i+1).center;
%             segment.axis1 = allArc(i+1).midPoint - segment.center;
%             segment.axis2 = cross(segment.axis1, cross(allArc(i+1).axis1, allArc(i+1).axis2)/(allArc(i+1).radius)^2);
%             segment.radius = allArc(i+1).radius;
%             segment.angle = allArc(i+1).angle;
%             segment.length = allArc(i+1).arcLen;
%             
%             v1 = -(PCC_result(end).axis1*cos(PCC_result(end).angle) + PCC_result(end).axis2*sin(PCC_result(end).angle));
%             v2 = -segment.axis1;
%             cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
%             sin_angle = norm(cross(v1, v2)) / (norm(v1) * norm(v2));
%             segment.rotation = atan2(norm(cross(v1, v2)), dot(v1, v2));
%             
%             PCC_result = [PCC_result segment];   
%         end
%     end  
% end

PCC_figure = figure('Name', 'Curve Segmentation', 'Visible','off', "NumberTitle", "off");
screenSize = get(groot, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
% Set the figure position
figWidth = 560;  % Set the width of the figure
figHeight = 420; % Set the height of the figure
PCC_figure.Position =  [(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight];

hold on
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
c = jet;
for i = 1:size(PCC_result,2)
    if PCC_result(i).type == 1
        line = [PCC_result(i).center + PCC_result(i).axis1; PCC_result(i).center - PCC_result(i).axis1];
        line_plot = plot3(line(:,1), line(:,2), line(:,3), '-', 'Color', c(rem(i*50,256),:), 'LineWidth',10);
        line_plot.Color(4) = 0.5;
    else
        theta = linspace(0, PCC_result(i).angle, 20);
        arc = [];
        for n_theta = 1:length(theta)
            arc = [arc; PCC_result(i).center + PCC_result(i).axis1*cos(theta(n_theta)) + PCC_result(i).axis2*sin(theta(n_theta))];
        end
        arc_plot = plot3(arc(:,1), arc(:,2), arc(:,3), '-', 'Color', c(rem(i*50,256),:), 'LineWidth',10);
        arc_plot.Color(4) = 0.5;
    end
end
bspline_curve = plot3(M(1,:), M(2,:), M(3,:), 'k-', 'LineWidth', 3);

% add legend
if ~isempty(PCC_result)
    str = strings(1, size(PCC_result, 2)+1);
    for i = 1:size(PCC_result, 2)
        str(i) = 'Segment ' + string(i);
    end
    str(end) = 'Original Curve';
    legend(str)
end

hold off
end


function error = biarcError(M, T, iIndex, eIndex, tol)
Arc = biarc(M, T, iIndex, eIndex, tol);
E = [];
midE = [];
for i = 1:eIndex-iIndex-1
    E = [E [Point2ArcError(Arc(1), M(:,iIndex+i)', tol); Point2ArcError(Arc(2), M(:,iIndex+i)', tol)]];
%     midE = [midE norm(Arc(2).midPoint - M(:,iIndex+j)')];
end
E = min(E,[],1);
% error = E;
error = mean(E);
end

function Arc = biarc(M, T, iIndex, eIndex, tol)
% get biarc with equal d between the end points and control points

p1 = M(:,iIndex)';
p2 = M(:,eIndex)';
t1 = T(:,iIndex)';
t2 = T(:,eIndex)';
Arc = struct('center',{[0,0,0], [0,0,0]},'axis1',{[0,0,0], [0,0,0]},'axis2',{[0,0,0], [0,0,0]},'radius',{0,0},'angle',{0,0},'arcLen',{0,0}, 'midPoint', {[0,0,0], [0,0,0]});

v = p2 - p1;
% if control points are equal, we don't need to interpolate
if (dot(v,v) < tol)
    Arc(1).center = p1;
    Arc(1).radius = 0;
    Arc(1).axis1 = v;
    Arc(1).axis2 = v;
    Arc(1).angle = 0;
    Arc(1).arcLen = 0;
    Arc(1).midPoint = p1;
    
    Arc(2).center = p1;
    Arc(2).radius = 0;
    Arc(2).axis1 = v;
    Arc(2).axis2 = v;
    Arc(2).angle = 0;
    Arc(2).arcLen = 0;
    Arc(2).midPoint = p1;
    return 
end


% d1 = d2
t = t1 + t2;
DEN = 2*(1-dot(t1,t2));
if DEN < tol
    if (abs(dot(v,t2)) < tol)
        % d being infinity, arc1,2 be semicircles
        % first arc
        pm = (p1 + p2)/2;

        p1ToPm = pm - p1;
        planeNormal1 = cross(p1ToPm,t1);
        perpAxis1 = cross(planeNormal1,p1ToPm);
        Arc(1).center = (pm+p1)/2;
        Arc(1).radius = norm(p1ToPm)/2;
        Arc(1).axis1 = p1 - Arc(1).center;
        Arc(1).axis2 = perpAxis1*(-Arc(1).radius)/dot(p1ToPm,p1ToPm);
        Arc(1).angle = pi;
        Arc(1).arcLen = pi*norm(p1ToPm)/2;
        Arc(1).midPoint = pm;
        
        % second arc 
        p2ToPm = pm - p2;
        planeNormal2 = cross(p2ToPm,t2);
        perpAxis2 = cross(planeNormal2,p2ToPm);
        Arc(2).center = (pm+p2)/2;
        Arc(2).radius = norm(p2ToPm)/2;
        Arc(2).axis1 = p2 - Arc(2).center;
        Arc(2).axis2 = perpAxis2*(-Arc(2).radius)/dot(p2ToPm,p2ToPm);
        Arc(2).angle = pi;
        Arc(2).arcLen = pi*norm(p2ToPm)/2;
        Arc(2).midPoint = pm;
        return;
    else
        d = dot(v,v)/(4*dot(v,t2));
    end
else
    disc = dot(v,t)^2 + DEN*dot(v,v);
    d = (-dot(v,t) + sqrt(disc))/DEN;
end

pm = ((p1+d*t1) + (p2-d*t2))/2;
p1ToPm = pm - p1;
p2ToPm = pm - p2;
[Arc(1).center, Arc(1).radius, angle1] = BiarcInterp_ComputeArc(p1, t1, p1ToPm, tol);
[Arc(2).center, Arc(2).radius, angle2] = BiarcInterp_ComputeArc(p2, t2, p2ToPm, tol);

if d < 0
    angle1 = sign(angle1)*2*pi - angle1;
    angle2 = sign(angle2)*2*pi - angle2;
end

Arc(1).angle = angle1;
Arc(1).axis1 = p1 - Arc(1).center;
Arc(1).axis2 = t1*Arc(1).radius;
if Arc(1).radius == 0
    Arc(1).arcLen = norm(p1ToPm);
else
    Arc(1).arcLen = angle1*Arc(1).radius;
end
Arc(1).midPoint = pm;

Arc(2).angle = angle2;
Arc(2).axis1 = p2 - Arc(2).center;
Arc(2).axis2 = -t2*Arc(2).radius;
if Arc(2).radius == 0
    Arc(2).arcLen = norm(p2ToPm);
else
    Arc(2).arcLen = angle2*Arc(2).radius;
end
Arc(2).midPoint = pm;

end


function [center, radius, angle] = BiarcInterp_ComputeArc(point, tangent, pointToMid, tol)
% compute an arc
    
    normal = cross(pointToMid, tangent); % the normal to the arc plane
    perpAxis = cross(tangent, normal)/norm(cross(tangent, normal));
    % axis within the arc plane that is perpendicular to the tangent,
    % colinear with the vector from the center to the end point
    
    if abs(2*dot(perpAxis, pointToMid)) < tol
        % radius is infinite, so use a straight line
        center = point + 0.5*pointToMid;
        radius = 0;
        angle = 0;
    else
        centerDist = dot(pointToMid,pointToMid)/(2*dot(perpAxis, pointToMid)); % the distance to the center along perpAxis
        center = point + centerDist*perpAxis;
        radius = abs(centerDist*norm(perpAxis));
        if radius < tol
            angle = 0; % compact arc
        else
            centerToEndDir = (point - center)/radius;
            centerToMidDir = (point - center + pointToMid)/radius;
            twist = dot(perpAxis, pointToMid);
            angle = acos(dot(centerToEndDir, centerToMidDir))*sign(twist);
        end
    end
end

function pResult = ArcIntep(arc, frac, tol)
if arc.arcLen < tol
    % just output the end point
    pResult = arc.center + arc.axis1;
else

    if arc.radius == 0
        % interpolate along the line from c+axis1 to c-axis1
        pResult = arc.center + (1-2*frac)*arc.axis1;
    else
        angle = arc.angle*frac;
        sinRot = sin(angle);
        cosRot = cos(angle);
        pResult = arc.center + arc.axis1*cosRot + arc.axis2*sinRot;
    end
end
end

function pResult = BiarcIntep(arc1, arc2, frac, tol)

    totalDist = arc1.arcLen + arc2.arcLen;
    fracDist = frac * totalDist;
    if (fracDist < arc1.arcLen)
        % on arc1
        if arc1.arcLen < tol
            % just output the end point
            pResult = arc1.center + arc1.axis1;
        else
            arcFrac = fracDist / arc1.arcLen;
            if arc1.radius == 0
                % interpolate along the line from c+axis1 to c-axis1
                pResult = arc1.center + (1-2*arcFrac)*arc1.axis1;
            else
                angle = arc1.angle*arcFrac;
                sinRot = sin(angle);
                cosRot = cos(angle);
                pResult = arc1.center + arc1.axis1*cosRot + arc1.axis2*sinRot;
            end
        end
    else
        if arc2.arcLen < tol
            pResult = arc2.center + arc2.axis1;
        else
            arcFrac = (fracDist - arc1.arcLen) / arc2.arcLen;
            if arc2.radius == 0
                pResult = arc2.center + (arcFrac*2-1)*arc2.axis1;
            else
                angle = arc2.angle * (1-arcFrac);
                sinRot = sin(angle);
                cosRot = cos(angle);
                pResult = arc2.center + arc2.axis1*cosRot + arc2.axis2*sinRot;
            end
        end
    end
end

function E = Point2ArcError(arc, point, tol)

    if arc.arcLen < tol
        % compact arc
        E = norm(point - arc.center);
    else
        if arc.radius == 0
            % a stright line
            iniPoint = arc.center + arc.axis1;
            endPoint = arc.center - arc.axis1;
            E = min([norm(point-iniPoint), norm(point-endPoint), Point2LineDist(iniPoint, endPoint, point)]);
        else
            % a normal arc
            iniPoint = arc.center + arc.axis1;
            endPoint = arc.center + arc.axis1*cos(arc.angle) + arc.axis2*sin(arc.angle);
            projPoint2ArcPlane = ProjPoint2Plane(arc.axis1, arc.axis2, iniPoint, point);
            sameDir1 = dot(cross(endPoint-arc.center, arc.axis1),cross(endPoint-arc.center, projPoint2ArcPlane-arc.center));
            sameDir2 = dot(cross(arc.axis1, endPoint-arc.center),cross(arc.axis1,projPoint2ArcPlane-arc.center));
            angle1 = acos(dot(endPoint-arc.center,projPoint2ArcPlane-arc.center)/norm(endPoint-arc.center)/norm(projPoint2ArcPlane-arc.center));
            angle2 = acos(dot(iniPoint-arc.center,projPoint2ArcPlane-arc.center)/norm(iniPoint-arc.center)/norm(projPoint2ArcPlane-arc.center));
            sameAngle = abs(angle1 + angle2 - arc.angle);
            if (sameDir1 >= 0) && (sameDir2 >= 0) && (sameAngle < tol)
                E = sqrt((norm(point-projPoint2ArcPlane))^2 + (arc.radius - norm(projPoint2ArcPlane - arc.center))^2);
            else
                E = min(norm(point-iniPoint), norm(point-iniPoint));
            end    
        end
    end
end

function d = Point2LineDist(A,B,P)
    num = 0;
    den = 0;
    for i = 1:length(A)
        num = num - (B(i)-A(i))*(A(i)-P(i));
        den = den + (B(i)-A(i))^2;
    end
    t = num/den;
    ProjectionPoint = A + t*(B-A);
    d = norm(ProjectionPoint-P);
end

function ProjectionPoint = ProjPoint2Plane(V1,V2,A,P)
    n = cross(V1,V2);
    n = n/norm(n);
    ProjectionPoint = P - dot(P - A, n) * n;
end
