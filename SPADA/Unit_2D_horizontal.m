function out = Unit_2D_horizontal(ri,t,rm,l)
% InnerRadius, WallThickness, AverageRadius, UnitLength in mm
SchematicUnit = figure('visible','off');
hold on
ro = rm*2 - ri;
r1 = l/4;
r2 = l/4;
f = ro - ri - r1 - r2;

Xlim = [-ro-t, ro+t];
Ylim = [-l/2, l*3/2];

xlim(Xlim);
% ylim(Ylim);


% middle bellow

LeftLowerArcCenterX = ri + r2;
LeftLowerArcCenterY = 0;
UpperArcCenterX = ro - r1;
UpperArcCenterY = -l/2;
RightLowerArcCenterX = ri + r2;
RightLowerArcCenterY = -l;

ThetaLeftLower = pi: 0.001: pi*3/2;
ThetaUpper = pi/2: -0.001: -pi/2;
ThetaRightLower = pi/2:0.001:pi;

mbx_1 = LeftLowerArcCenterX + r2*cos(ThetaLeftLower);
mby_1 = LeftLowerArcCenterY + r2*sin(ThetaLeftLower);
mbx_2 = [mbx_1(end), mbx_1(end)+f];
mby_2 = mby_1(end)*ones(1,2);
mbx_3 = UpperArcCenterX + r1*cos(ThetaUpper);
mby_3 = UpperArcCenterY + r1*sin(ThetaUpper);
mbx_4 = [mbx_3(end), mbx_3(1)-f];
mby_4 = mby_3(end)*ones(1,2);
mbx_5 = RightLowerArcCenterX + r2*cos(ThetaRightLower);
mby_5 = RightLowerArcCenterY + r2*sin(ThetaRightLower);

mbx = [mbx_1 mbx_2 mbx_3 mbx_4 mbx_5];
mby = [mby_1 mby_2 mby_3 mby_4 mby_5];
% plot(mbx, mby, 'k--', 'LineWidth', 2)
% plot(-mbx, mby, 'k--', 'LineWidth', 2)
% 
% outer bellow
obx_1 = LeftLowerArcCenterX + (r2-t/2)*cos(ThetaLeftLower);
oby_1 = LeftLowerArcCenterY + (r2-t/2)*sin(ThetaLeftLower);
obx_2 = [obx_1(end), obx_1(end)+f];
oby_2 = oby_1(end)*ones(1,2);
obx_3 = UpperArcCenterX + (r1+t/2)*cos(ThetaUpper);
oby_3 = UpperArcCenterY + (r1+t/2)*sin(ThetaUpper);
obx_4 = [obx_3(end), obx_3(1)-f];
oby_4 = oby_3(end)*ones(1,2);
obx_5 = RightLowerArcCenterX + (r2-t/2)*cos(ThetaRightLower);
oby_5 = RightLowerArcCenterY + (r2-t/2)*sin(ThetaRightLower);

obx = [obx_1 obx_2 obx_3 obx_4 obx_5];
oby = [oby_1 oby_2 oby_3 oby_4 oby_5];
% plot(obx, oby, 'k-', 'LineWidth', 2)
% plot(-obx, oby, 'k-', 'LineWidth', 2)
% 
% inner bellow
ibx_1 = LeftLowerArcCenterX + (r2+t/2)*cos(ThetaLeftLower);
iby_1 = LeftLowerArcCenterY + (r2+t/2)*sin(ThetaLeftLower);
ibx_2 = [ibx_1(end), ibx_1(end)+f];
iby_2 = iby_1(end)*ones(1,2);
ibx_3 = UpperArcCenterX + (r1-t/2)*cos(ThetaUpper);
iby_3 = UpperArcCenterY + (r1-t/2)*sin(ThetaUpper);
ibx_4 = [ibx_3(end), ibx_3(1)-f];
iby_4 = iby_3(end)*ones(1,2);
ibx_5 = RightLowerArcCenterX + (r2+t/2)*cos(ThetaRightLower);
iby_5 = RightLowerArcCenterY + (r2+t/2)*sin(ThetaRightLower);

ibx = [ibx_1 ibx_2 ibx_3 ibx_4 ibx_5];
iby = [iby_1 iby_2 iby_3 iby_4 iby_5];
% % plot(ibx, iby, 'k-', 'LineWidth', 2)
% % plot(-ibx, iby, 'k-', 'LineWidth', 2)
% 
bx = [obx fliplr(ibx)];
by = [oby fliplr(iby)];
fill(bx,by,[0.5 0.5 0.5],'LineStyle','none', 'FaceAlpha', 0.5)
bx = [-obx fliplr(-ibx)];
by = [oby fliplr(iby)];
fill(bx,by,[0.5 0.5 0.5],'LineStyle','none', 'FaceAlpha', 0.5)
% 
% constraint
cx = [-t -t fliplr(-ibx)];
cy = [0 -l fliplr(iby)];
% plot([-t -t],[0 -l], 'k-', 'LineWidth', 2);
fill(cx,cy,[0.9 0.85 0],'LineStyle','none', 'FaceAlpha', 0.5)
% 
% 
% rest part
fill([-t -t fliplr(obx)],[0 -l fliplr(oby)],[0.5 0.5 0.5],'LineStyle','none', 'FaceAlpha', 0.3)
% 
% central axis
plot([0 0], [0 -l], 'k--', 'LineWidth', 2);
% 
% line plot
plot(mbx, mby, 'k--', 'LineWidth', 2)
plot(-mbx, mby, 'k--', 'LineWidth', 2)
plot(obx, oby, 'k-', 'LineWidth', 2)
plot(-obx, oby, 'k-', 'LineWidth', 2)
plot(ibx, iby, 'k-', 'LineWidth', 2)
plot(-ibx, iby, 'k-', 'LineWidth', 2)
plot([-t -t],[0 -l], 'k-', 'LineWidth', 2);
% 
% 
% dimension

% 
% 

% ri
quiver(0,0,ri,0, 0, 'linewidth',2,'color','r','MaxHeadSize',0.5)
text(t/2,-r2/4, '$r_{in}$', 'FontSize', 14, 'color' , 'r', 'Interpreter', 'latex')
% ro
quiver(0,-l/2,ro,0, 0, 'linewidth',2,'color','r','MaxHeadSize',0.25)
text(t/2,-(l/2+r1/3), {'$r_{ou}=2R - r_{in}$'}, 'FontSize', 14, 'color' , 'r', 'Interpreter', 'latex')
%l
quiver(0,-l/2,0,l/2, 0, 'linewidth',2,'color','r','MaxHeadSize',0.5)
quiver(0,-l/2,0,-l/2, 0, 'linewidth',2,'color','r','MaxHeadSize',0.5)
text(-t/2,-l/2, '$l$', 'FontSize', 14, 'color' , 'r', 'Interpreter', 'latex')
%t
quiver(LeftLowerArcCenterX,LeftLowerArcCenterY,-(r2-t/2)*sqrt(2)/2,-(r2-t/2)*sqrt(2)/2, 0, 'linewidth',2,'color','r','MaxHeadSize',2)
quiver(LeftLowerArcCenterX-(2*r2)*sqrt(2)/2,LeftLowerArcCenterY-(2*r2)*sqrt(2)/2,(r2-t/2)*sqrt(2)/2,(r2-t/2)*sqrt(2)/2, 0, 'linewidth',2,'color','r','MaxHeadSize',2)
text(LeftLowerArcCenterX+r2/4,LeftLowerArcCenterY-r2/4, '$t$', 'FontSize', 14, 'color' , 'r', 'Interpreter', 'latex')

% 
axis equal
maxl = max(l, ro+t);
xlim([0-maxl, 0+maxl]);
ylim([-l,0])

box on
grid on

% set(gca, 'visible', 'off')
% set(gca,'XColor', 'none','YColor','none')
% set(gca, 'PaperUnits', 'centimeters');
% set(gca, 'PaperPosition', [0 0 l*1.2 (ro+t)*2]);
% saveas(SchematicUnit, 'SchematicUnit.jpg')
exportgraphics(SchematicUnit, 'SchematicUnit.jpg')
% set(SchematicUnit, 'visible', 'on'); 
hold off

% close(SchematicUnit)
end