function createStackfigure(X1, ymatrix1)
%CREATEFIGURE(X1,YMATRIX1)
%  X1:  area x
%  YMATRIX1:  area matrix data

%  Auto-generated by MATLAB on 15-Dec-2009 15:04:21

% Create figure
figure1 = figure('XVisual',...
    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'InvertHardcopy','off',...
    'Color',[1 1 1]);

% Create axes
axes('Parent',figure1,'LineWidth',3,'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
% xlim([0 10]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim([0 1.6e-09]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim([-1 1]);
box('on');
hold('all');

% Create multiple lines using matrix input to area
area1 = area(X1,ymatrix1);
set(area1(1),'FaceColor',[1 0 0]);
set(area1(2),'FaceColor',[0 1 0]);
axis([0 10 0 0.5e-8]);


% Create xlabel
xlabel('Time (hr)','FontSize',16);

% Create ylabel
ylabel('Coating mass (ug)','FontSize',16);
%title('Mass on emitted particle');

