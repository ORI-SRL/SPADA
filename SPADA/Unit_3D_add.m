function [AssemblyFigure, SurfHandle, LineHandle] = Unit_3D_add(x,y,z,new_unit_data,AssemblyFigure,SurfHandle,LineHandle)


phi = -new_unit_data(1);
ri = new_unit_data(2);
t = new_unit_data(3);
rm = new_unit_data(4);
l = new_unit_data(5);

Unit3D = figure(AssemblyFigure);
set(AssemblyFigure,'visible','off')


hold on
axis equal
view(30,30);
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');

ro = rm*2 - ri;
r1 = l/4;
r2 = l/4;
f = ro - ri - r1 - r2;

% outer bellow
Ro = 0:(pi*2-0)/100:pi*2;

LeftLowerArcCenterX = x + (ri+r2)*cos(Ro);
LeftLowerArcCenterY = y + (ri+r2)*sin(Ro);
LeftLowerArcCenterZ = z;
UpperArcCenterX = x + (ro-r1)*cos(Ro);
UpperArcCenterY = y + (ro-r1)*sin(Ro);
UpperArcCenterZ = z - l/2;
RightLowerArcCenterX = x + (ri+r2)*cos(Ro);
RightLowerArcCenterY = y + (ri+r2)*sin(Ro);
RightLowerArcCenterZ = z - l;

ThetaLeftLower = pi: (pi*3/2-pi)/100: pi*3/2;
ThetaUpper = pi/2: (-pi/2-pi/2)/100: -pi/2;
ThetaRightLower = pi/2:(pi-pi/2)/100: pi;


obx_1 = LeftLowerArcCenterX + (r2-t/2)*transpose(cos(ThetaLeftLower))*cos(Ro);
oby_1 = LeftLowerArcCenterY + (r2-t/2)*transpose(cos(ThetaLeftLower))*sin(Ro);
obz_1 = LeftLowerArcCenterZ + (r2-t/2)*repmat(transpose(sin(ThetaLeftLower)),1,length(sin(Ro)));
obx_2 = [obx_1(end,:); obx_1(end,:)+f*(cos(Ro))];
oby_2 = [oby_1(end,:); oby_1(end,:)+f*(sin(Ro))];
obz_2 = [obz_1(end,:); obz_1(end,:)];
obx_3 = UpperArcCenterX + (r1+t/2)*transpose(cos(ThetaUpper))*cos(Ro);
oby_3 = UpperArcCenterY + (r1+t/2)*transpose(cos(ThetaUpper))*sin(Ro);
obz_3 = UpperArcCenterZ + (r1+t/2)*repmat(transpose(sin(ThetaUpper)),1,length(sin(Ro)));
obx_4 = [obx_3(end,:); obx_3(end,:)-f*(cos(Ro))];
oby_4 = [oby_3(end,:); oby_3(end,:)-f*(sin(Ro))];
obz_4 = [obz_3(end,:); obz_3(end,:)];
obx_5 = RightLowerArcCenterX + (r2-t/2)*transpose(cos(ThetaRightLower))*cos(Ro);
oby_5 = RightLowerArcCenterY + (r2-t/2)*transpose(cos(ThetaRightLower))*sin(Ro);
obz_5 = RightLowerArcCenterZ + (r2-t/2)*repmat(transpose(sin(ThetaRightLower)),1,length(sin(Ro)));
obx = [obx_1; obx_2; obx_3; obx_4; obx_5];
oby = [oby_1; oby_2; oby_3; oby_4; oby_5];
obz = [obz_1; obz_2; obz_3; obz_4; obz_5];

ob_surf = surf(obx,oby,obz,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.5,'EdgeColor','none');

% inner bellow
Ri = (phi+pi*3/4):(-pi*3/2)/100:(phi-pi*3/4);

LeftLowerArcCenterX = x + (ri+r2)*cos(Ri);
LeftLowerArcCenterY = y + (ri+r2)*sin(Ri);
LeftLowerArcCenterZ = z;
UpperArcCenterX = x + (ro-r1)*cos(Ri);
UpperArcCenterY = y + (ro-r1)*sin(Ri);
UpperArcCenterZ = z - l/2;
RightLowerArcCenterX = x + (ri+r2)*cos(Ri);
RightLowerArcCenterY = y + (ri+r2)*sin(Ri);
RightLowerArcCenterZ = z - l;

ThetaLeftLower = pi: (pi*3/2-pi)/100: pi*3/2;
ThetaUpper = pi/2: (-pi/2-pi/2)/100: -pi/2;
ThetaRightLower = pi/2:(pi-pi/2)/100: pi;

ibx_1 = LeftLowerArcCenterX + (r2+t/2)*transpose(cos(ThetaLeftLower))*cos(Ri);
iby_1 = LeftLowerArcCenterY + (r2+t/2)*transpose(cos(ThetaLeftLower))*sin(Ri);
ibz_1 = LeftLowerArcCenterZ + (r2+t/2)*repmat(transpose(sin(ThetaLeftLower)),1,length(sin(Ri)));
ibx_2 = [ibx_1(end,:); ibx_1(end,:)+f*(cos(Ri))];
iby_2 = [iby_1(end,:); iby_1(end,:)+f*(sin(Ri))];
ibz_2 = [ibz_1(end,:); ibz_1(end,:)];
ibx_3 = UpperArcCenterX + (r1-t/2)*transpose(cos(ThetaUpper))*cos(Ri);
iby_3 = UpperArcCenterY + (r1-t/2)*transpose(cos(ThetaUpper))*sin(Ri);
ibz_3 = UpperArcCenterZ + (r1-t/2)*repmat(transpose(sin(ThetaUpper)),1,length(sin(Ri)));
ibx_4 = [ibx_3(end,:); ibx_3(end,:)-f*(cos(Ri))];
iby_4 = [iby_3(end,:); iby_3(end,:)-f*(sin(Ri))];
ibz_4 = [ibz_3(end,:); ibz_3(end,:)];
ibx_5 = RightLowerArcCenterX + (r2+t/2)*transpose(cos(ThetaRightLower))*cos(Ri);
iby_5 = RightLowerArcCenterY + (r2+t/2)*transpose(cos(ThetaRightLower))*sin(Ri);
ibz_5 = RightLowerArcCenterZ + (r2+t/2)*repmat(transpose(sin(ThetaRightLower)),1,length(sin(Ri)));
ibx = [ibx_1; ibx_2; ibx_3; ibx_4; ibx_5];
iby = [iby_1; iby_2; iby_3; iby_4; iby_5];
ibz = [ibz_1; ibz_2; ibz_3; ibz_4; ibz_5];

ib_surf = surf(ibx,iby,ibz,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.5,'EdgeColor','none');

% constraint
C = (phi+pi*3/4):pi/2/100:(phi+pi*5/4);

LeftLowerArcCenterX = x + (ri+r2)*cos(C);
LeftLowerArcCenterY = y + (ri+r2)*sin(C);
LeftLowerArcCenterZ = z;
UpperArcCenterX = x + (ro-r1)*cos(C);
UpperArcCenterY = y + (ro-r1)*sin(C);
UpperArcCenterZ = z - l/2;
RightLowerArcCenterX = x + (ri+r2)*cos(C);
RightLowerArcCenterY = y + (ri+r2)*sin(C);
RightLowerArcCenterZ = z - l;

cx_1 = LeftLowerArcCenterX + (r2+t/2)*transpose(cos(ThetaLeftLower))*cos(C);
cy_1 = LeftLowerArcCenterY + (r2+t/2)*transpose(cos(ThetaLeftLower))*sin(C);
cz_1 = LeftLowerArcCenterZ + (r2+t/2)*repmat(transpose(sin(ThetaLeftLower)),1,length(sin(C)));
cx_2 = [cx_1(end,:); cx_1(end,:)+f*(cos(C))];
cy_2 = [cy_1(end,:); cy_1(end,:)+f*(sin(C))];
cz_2 = [cz_1(end,:); cz_1(end,:)];
cx_3 = UpperArcCenterX + (r1-t/2)*transpose(cos(ThetaUpper))*cos(C);
cy_3 = UpperArcCenterY + (r1-t/2)*transpose(cos(ThetaUpper))*sin(C);
cz_3 = UpperArcCenterZ + (r1-t/2)*repmat(transpose(sin(ThetaUpper)),1,length(sin(C)));
cx_4 = [cx_3(end,:); cx_3(end,:)-f*(cos(C))];
cy_4 = [cy_3(end,:); cy_3(end,:)-f*(sin(C))];
cz_4 = [cz_3(end,:); cz_3(end,:)];
cx_5 = RightLowerArcCenterX + (r2+t/2)*transpose(cos(ThetaRightLower))*cos(C);
cy_5 = RightLowerArcCenterY + (r2+t/2)*transpose(cos(ThetaRightLower))*sin(C);
cz_5 = RightLowerArcCenterZ + (r2+t/2)*repmat(transpose(sin(ThetaRightLower)),1,length(sin(C)));
cx = [cx_1; cx_2; cx_3; cx_4; cx_5; x+t*cos(C);x+t*cos(C); cx_1(1,:)];
cy = [cy_1; cy_2; cy_3; cy_4; cy_5; y+t*sin(C);y+t*sin(C); cy_1(1,:)];
cz = [cz_1; cz_2; cz_3; cz_4; cz_5; z-l*(ones(1,length(C))); z-0*(ones(1,length(C))); cz_1(1,:)];

const_surf = surf(cx,cy,cz,'FaceColor',[0.9 0.85 0], 'FaceAlpha',0.5,'EdgeColor','none');

% connect constraint and spline
csx_1 = [transpose(cx_1(:,1)) transpose(cx_3(:,1)) transpose(cx_5(:,1))];
csx_1 = [csx_1; x+t*cos(C(1))*ones(1,length(csx_1))];
csy_1 = [transpose(cy_1(:,1)) transpose(cy_3(:,1)) transpose(cy_5(:,1))];
csy_1 = [csy_1; y+t*sin(C(1))*ones(1,length(csy_1))];
csz_1 = [transpose(cz_1(:,1)) transpose(cz_3(:,1)) transpose(cz_5(:,1))];
csz_1 = [csz_1; csz_1];

const_sp_surf_1 = surf(csx_1,csy_1,csz_1,'FaceColor',[0.9 0.85 0], 'FaceAlpha',0.5,'EdgeColor','none');

csx_2 = [transpose(cx_1(:,end)) transpose(cx_3(:,end)) transpose(cx_5(:,end))];
csx_2 = [csx_2; x+t*cos(C(end))*ones(1,length(csx_2))];
csy_2 = [transpose(cy_1(:,end)) transpose(cy_3(:,end)) transpose(cy_5(:,end))];
csy_2 = [csy_2; y+t*sin(C(end))*ones(1,length(csy_2))];
csz_2 = [transpose(cz_1(:,end)) transpose(cz_3(:,end)) transpose(cz_5(:,end))];
csz_2 = [csz_2; csz_2];

const_sp_surf_2 = surf(csx_2,csy_2, csz_2,'FaceColor',[0.9 0.85 0], 'FaceAlpha',0.5,'EdgeColor','none');

% top and end cap
R = 0:(pi*2-0)/100:pi*2;
tx = x + [(ri-t/2)*cos(R); (ri+t/2)*cos(R)];
ty = y + [(ri-t/2)*sin(R); (ri+t/2)*sin(R)];
tz = z - zeros(size(tx));

top_surf = surf(tx,ty,tz,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.5,'EdgeColor','none');

ez = z - l*ones(size(tx));
end_surf = surf(tx,ty,ez,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.5,'EdgeColor','none');


% add edges
% top cap
tex_1 = x + (ri-t/2)*cos(R);
tey_1 = y + (ri-t/2)*sin(R);
tez_1 = z - zeros(size(tex_1));
top_edge_1 = plot3(tex_1,tey_1,tez_1, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
tex_2 = x + (ri+t/2)*cos(R);
tey_2 = y + (ri+t/2)*sin(R);
tez_2 = z - zeros(size(tex_2));
top_edge_2 = plot3(tex_2,tey_2,tez_2, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
eex_1 = x + (ri-t/2)*cos(R);
eey_1 = y + (ri-t/2)*sin(R);
eez_1 = z - l*ones(size(eex_1));
end_edge_1 = plot3(eex_1,eey_1,eez_1, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
eex_2 = x + (ri+t/2)*cos(R);
eey_2 = y + (ri+t/2)*sin(R);
eez_2 = z - l*ones(size(eex_2));
end_edge_2 = plot3(eex_2,eey_2,eez_2, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
% constraint
cex_1 = [x+t*cos(C(1)); cx_1(:,1); cx_2(:,1); cx_3(:,1); cx_4(:,1); cx_5(:,1); x+t*cos(C(1)); x+t*cos(C(1))];
cey_1 = [y+t*sin(C(1)); cy_1(:,1); cy_2(:,1); cy_3(:,1); cy_4(:,1); cy_5(:,1); y+t*sin(C(1)); y+t*sin(C(1))];
cez_1 = [z; cz_1(:,1); cz_2(:,1); cz_3(:,1); cz_4(:,1); cz_5(:,1); z-l; z];
const_edge_1 = plot3(cex_1,cey_1,cez_1, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
cex_2 = [x+t*cos(C(end)); cx_1(:,end); cx_2(:,end); cx_3(:,end); cx_4(:,end); cx_5(:,end); x+t*cos(C(end)); x+t*cos(C(end))];
cey_2 = [y+t*sin(C(end)); cy_1(:,end); cy_2(:,end); cy_3(:,end); cy_4(:,end); cy_5(:,end); y+t*sin(C(end)); y+t*sin(C(end))];
cez_2 = [z; cz_1(:,end); cz_2(:,end); cz_3(:,end); cz_4(:,end); cz_5(:,end); z-l; z];
const_edge_2 = plot3(cex_2,cey_2,cez_2, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
cex_3 = x + t*cos(C);
cey_3 = y + t*sin(C);
cez_3 = z - zeros(size(cex_3));
const_edge_3 = plot3((cex_3),(cey_3),(cez_3), '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);
cex_4 = x + t*cos(C);
cey_4 = y + t*sin(C);
cez_4 = z - l*ones(size(cex_4));
const_edge_4 = plot3(cex_4,cey_4,cez_4, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1);

% bellow edge
Roe = 0:(pi*2)/4:pi*2;
LeftLowerArcCenterX = x + (ri+r2)*cos(Roe);
LeftLowerArcCenterY = y + (ri+r2)*sin(Roe);
LeftLowerArcCenterZ = z;
UpperArcCenterX = x + (ro-r1)*cos(Roe);
UpperArcCenterY = y + (ro-r1)*sin(Roe);
UpperArcCenterZ = z - l/2;
RightLowerArcCenterX = x + (ri+r2)*cos(Roe);
RightLowerArcCenterY = y + (ri+r2)*sin(Roe);
RightLowerArcCenterZ = z - l;

ThetaLeftLower = pi: (pi*3/2-pi)/100: pi*3/2;
ThetaUpper = pi/2: (-pi/2-pi/2)/100: -pi/2;
ThetaRightLower = pi/2:(pi-pi/2)/100: pi;


obx_1 = LeftLowerArcCenterX + (r2-t/2)*transpose(cos(ThetaLeftLower))*cos(Roe);
oby_1 = LeftLowerArcCenterY + (r2-t/2)*transpose(cos(ThetaLeftLower))*sin(Roe);
obz_1 = LeftLowerArcCenterZ + (r2-t/2)*repmat(transpose(sin(ThetaLeftLower)),1,length(sin(Roe)));
obx_2 = [obx_1(end,:); obx_1(end,:)+f*(cos(Roe))];
oby_2 = [oby_1(end,:); oby_1(end,:)+f*(sin(Roe))];
obz_2 = [obz_1(end,:); obz_1(end,:)];
obx_3 = UpperArcCenterX + (r1+t/2)*transpose(cos(ThetaUpper))*cos(Roe);
oby_3 = UpperArcCenterY + (r1+t/2)*transpose(cos(ThetaUpper))*sin(Roe);
obz_3 = UpperArcCenterZ + (r1+t/2)*repmat(transpose(sin(ThetaUpper)),1,length(sin(Roe)));
obx_4 = [obx_3(end,:); obx_3(end,:)-f*(cos(Roe))];
oby_4 = [oby_3(end,:); oby_3(end,:)-f*(sin(Roe))];
obz_4 = [obz_3(end,:); obz_3(end,:)];
obx_5 = RightLowerArcCenterX + (r2-t/2)*transpose(cos(ThetaRightLower))*cos(Roe);
oby_5 = RightLowerArcCenterY + (r2-t/2)*transpose(cos(ThetaRightLower))*sin(Roe);
obz_5 = RightLowerArcCenterZ + (r2-t/2)*repmat(transpose(sin(ThetaRightLower)),1,length(sin(Roe)));
obx = [obx_1; obx_2; obx_3; obx_4; obx_5];
oby = [oby_1; oby_2; oby_3; oby_4; oby_5];
obz = [obz_1; obz_2; obz_3; obz_4; obz_5];
ob_edge = transpose(plot3(obx, oby, obz, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1));



Rie = 0:(pi*2)/4:pi*2;

LeftLowerArcCenterX = x + (ri+r2)*cos(Rie);
LeftLowerArcCenterY = y + (ri+r2)*sin(Rie);
LeftLowerArcCenterZ = z;
UpperArcCenterX = x + (ro-r1)*cos(Rie);
UpperArcCenterY = y + (ro-r1)*sin(Rie);
UpperArcCenterZ = z - l/2;
RightLowerArcCenterX = x + (ri+r2)*cos(Rie);
RightLowerArcCenterY = y + (ri+r2)*sin(Rie);
RightLowerArcCenterZ = z - l;

ThetaLeftLower = pi: (pi*3/2-pi)/100: pi*3/2;
ThetaUpper = pi/2: (-pi/2-pi/2)/100: -pi/2;
ThetaRightLower = pi/2:(pi-pi/2)/100: pi;

ibx_1 = LeftLowerArcCenterX + (r2+t/2)*transpose(cos(ThetaLeftLower))*cos(Rie);
iby_1 = LeftLowerArcCenterY + (r2+t/2)*transpose(cos(ThetaLeftLower))*sin(Rie);
ibz_1 = LeftLowerArcCenterZ + (r2+t/2)*repmat(transpose(sin(ThetaLeftLower)),1,length(sin(Rie)));
ibx_2 = [ibx_1(end,:); ibx_1(end,:)+f*(cos(Rie))];
iby_2 = [iby_1(end,:); iby_1(end,:)+f*(sin(Rie))];
ibz_2 = [ibz_1(end,:); ibz_1(end,:)];
ibx_3 = UpperArcCenterX + (r1-t/2)*transpose(cos(ThetaUpper))*cos(Rie);
iby_3 = UpperArcCenterY + (r1-t/2)*transpose(cos(ThetaUpper))*sin(Rie);
ibz_3 = UpperArcCenterZ + (r1-t/2)*repmat(transpose(sin(ThetaUpper)),1,length(sin(Rie)));
ibx_4 = [ibx_3(end,:); ibx_3(end,:)-f*(cos(Rie))];
iby_4 = [iby_3(end,:); iby_3(end,:)-f*(sin(Rie))];
ibz_4 = [ibz_3(end,:); ibz_3(end,:)];
ibx_5 = RightLowerArcCenterX + (r2+t/2)*transpose(cos(ThetaRightLower))*cos(Rie);
iby_5 = RightLowerArcCenterY + (r2+t/2)*transpose(cos(ThetaRightLower))*sin(Rie);
ibz_5 = RightLowerArcCenterZ + (r2+t/2)*repmat(transpose(sin(ThetaRightLower)),1,length(sin(Rie)));
ibx = [ibx_1; ibx_2; ibx_3; ibx_4; ibx_5];
iby = [iby_1; iby_2; iby_3; iby_4; iby_5];
ibz = [ibz_1; ibz_2; ibz_3; ibz_4; ibz_5];
ib_edge = transpose(plot3(ibx, iby, ibz, '-', 'Color',[0,0,0,0.75], 'LineWidth', 1));



SurfHandle = [SurfHandle; [ob_surf ib_surf const_surf const_sp_surf_1 const_sp_surf_2 top_surf end_surf]];
LineHandle = [LineHandle; [top_edge_1 top_edge_2 end_edge_1 end_edge_2 const_edge_1 const_edge_2 const_edge_3 const_edge_4 ob_edge ib_edge]];

hold off

end