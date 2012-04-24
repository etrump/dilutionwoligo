figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

modelAtm.DilutionTime = 0;

global Tagg DpT modelAtm
dimerbin = 8;
MaxIndex = length(find(Tagg<modelAtm.DilutionTime));
if MaxIndex < 1
    MaxIndex = 1;
end
V0 = DpT(MaxIndex,1)^3
%VFR = DpT(:,1).^3/V0
Tagg_hr = Tagg/3600;
DpT_nm = DpT;


CpT_integrate = zeros(length(Tagg),modelAtm.NumBins);

for k = 1:length(Tagg)
for j=1:modelAtm.Pop % loop thru particle populations
    for i=1:modelAtm.NumBins  % loop thru organic constituents
        CpT_integrate(k,i) = CpT_integrate(k,i)+CpT(k,(j-1)*modelAtm.NumBins+i);
    end
    CpTot_integrate(k) = sum(CpT_integrate(k,:));
end
end

CpTot_integrate0 = CpTot_integrate(1);
VFR = CpTot_integrate/CpTot_integrate0;


plot1 = plot(Tagg_hr, DpT_nm(:,1),'LineWidth',3)

set(plot1(1),'DisplayName','Particle','Color',[0 0 1]);
%set(plot1(2),'DisplayName','Nanoparticle','Color',[0 0 1]);


% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'D_{p,ave} (nm)'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)

% Create legend
paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

legend1 = legend(axes1,'hide');
%%axis([0 5  50 90]); %BASE
%axis([0 15  0 200]);

%axis([2 5 0 1]);
%set(legend1,'Position',[0.5034 0.1988 0.37 0.1601]);
set(legend1,'Location','NorthWest')

% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/Diam.pdf'],'pdf'); 



figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure2,'LineWidth',3,'FontSize',16);
box('on');
hold('all');


plot2 = plot(Tagg_hr, VFR,'LineWidth',3)
set(plot2(1),'DisplayName','Particle','Color',[0 0 1]);

axis([0 10  0 1.01]);
% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'VFR'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/VFR.pdf'],'pdf'); 

%VFRdata120326;
VFRdata120409;

figure3 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure3,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

plot3 = plot(Tagg/3600,(CpT_integrate(:,1))./(CvT(:,dimerbin)+CpT_integrate(:,1)+CpT_integrate(:,dimerbin)),'LineWidth',3);
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'olig_{part}/(olig_{part}+mon_{vap}+mon_{part})'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/partition_olig.pdf'],'pdf'); 


figure4 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure4,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

CpTot(:,1) = 2*CpT_integrate(:,1)+CpT_integrate(:,2)+CpT_integrate(:,3)+CpT_integrate(:,4)+CpT_integrate(:,5)+CpT_integrate(:,6)+CpT_integrate(:,7)+CpT_integrate(:,8)+CpT_integrate(:,9)+CpT_integrate(:,10);
plot4 = plot(Tagg/3600,CpTot(:,1),'LineWidth',3)
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'C_{P tot}'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/Cp_tot.pdf'],'pdf'); 



figure5 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure5,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

plot5 = plot(Tagg/3600,CpT_integrate(:,1),'LineWidth',3)
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'C_{P dimer}'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/C_dimer.pdf'],'pdf');


figure6 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure6,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

ystack = [CpT_integrate(:,1) CpT_integrate(:,2) CpT_integrate(:,3) CpT_integrate(:,4) CpT_integrate(:,5) CpT_integrate(:,6) CpT_integrate(:,7) ...
    CpT_integrate(:,8),CpT_integrate(:,9),CpT_integrate(:,10)];
%ystack = [CpT(:,1) CpT(:,8),CpT(:,9),CpT(:,10)];
area1 = area(Tagg/3600,ystack,'LineWidth',3);

set(area1(1),'DisplayName','lC*=-15');
set(area1(2),'DisplayName','lC*=-2');
set(area1(3),'DisplayName','1C*-1');
set(area1(4),'DisplayName','1C*=0');
set(area1(5),'DisplayName','1C*=1');
set(area1(6),'DisplayName','lC*=2');
set(area1(7),'DisplayName','lC*=3');
set(area1(8),'DisplayName','lC*=4');
set(area1(9),'DisplayName','lC*=5');
set(area1(10),'DisplayName','lC*=6');
legend1 = legend(axes1,'show');

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/C_p composition.pdf'],'pdf');


totalMassatDilution = 0;
for i=1:modelAtm.NumBins
    totalMassatDilution = totalMassatDilution + CpT(MaxIndex,i);
end

CompatDilution = CpT(MaxIndex,:)/totalMassatDilution
