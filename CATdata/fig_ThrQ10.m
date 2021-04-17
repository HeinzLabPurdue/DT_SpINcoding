% FILE: fig_ThrQ10.m
%     : STIMind: 1-8, or 0 for all units
%
% For: RLFpaper
% From: ISH2003paper
% From: IHCON2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('Fig4_workspace.mat')

CURdir=pwd;

% global Recruit_dir RLFpaper_dir

STIMind=0;
% global SEVimpExps MILDimpExps





%%%%%%%%%%% Thresholds
ymin=-10; ymax=110;
h1=subplot(231);     % Normal
semilogx(normt(1,:),normt(2,:),'k')
hold on

for SRind=NUM_SRgroups:-1:1
	semilogx(BF_N{SRind},Thr_N{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
end
hleg1=legend('Miller et al. (1997)',sprintf('high (SR>%.1f)', ...
	SRbounds(2)),sprintf('med (%.1f<SR<=%.1f)', ...
	SRbounds(1),SRbounds(2)),sprintf('low (SR<=%.1f)',SRbounds(1)));
axis([xmin xmax ymin ymax])
title('Normal','FontSize',LabelFontSize)
ylabel('Threshold (dB SPL)','FontSize',LabelFontSize)
text(.08,.15,'NBTC','FontSize',TextFontSize,'VerticalAlignment','top','HorizontalAlignment','left');

h2=subplot(232);     % Mild
for SRind=NUM_SRgroups:-1:1
   semilogx(BF_Imild{SRind},Thr_Imild{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
   hold on
end
semilogx(normt(1,:),normt(2,:),'k')
axis([xmin xmax ymin ymax])
title('Mild Loss','FontSize',LabelFontSize)

h3=subplot(233);     % Severe
for SRind=NUM_SRgroups:-1:1
	semilogx(BF_Isev{SRind},Thr_Isev{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
semilogx(normt(1,:),normt(2,:),'k')
axis([xmin xmax ymin ymax])
title('Moderate/Severe Loss','FontSize',LabelFontSize)

%%%%%%%%%%% Q10
ymin=.3; ymax=20;
h4=subplot(234);     % Normal
for SRind=1:NUM_SRgroups
	loglog(BF_N{SRind},Q10_N{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
axis([xmin xmax ymin ymax])
text(.7,.125,sprintf('N=%d',sum([length(BF_N{1}) length(BF_N{2}) length(BF_N{3})])),'Units','norm','FontSize',TextFontSize)
ylabel('Q_{10}','FontSize',LabelFontSize,'Interpreter','tex');

h5=subplot(235);     % Mild
for SRind=1:NUM_SRgroups
	loglog(BF_Imild{SRind},Q10_Imild{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
axis([xmin xmax ymin ymax])
text(.7,.125,sprintf('N=%d',sum([length(BF_Imild{1}) length(BF_Imild{2}) length(BF_Imild{3})])),'Units','norm','FontSize',TextFontSize)
xlabel('Best Frequency (kHz)','FontSize',LabelFontSize)

h6=subplot(236);     % Severe
for SRind=1:NUM_SRgroups
	loglog(BF_Isev{SRind},Q10_Isev{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
axis([xmin xmax ymin ymax])
text(.7,.125,sprintf('N=%d',sum([length(BF_Isev{1}) length(BF_Isev{2}) length(BF_Isev{3})])),'Units','norm','FontSize',TextFontSize)



Xwidth=.31;
Xcorner=.062;
Xshift=0;
Ywidth=.42;
Ycorner=0.523;
Yshift=0.015;

set(h1,'Position',[Xcorner Ycorner Xwidth Ywidth],'XTickLabel',[],'YTick',[0:20:100], ...\
   'FontSize',AnnotFontSize,'TickLength',[TICKlength 0.025])
set(hleg1,'Position',[Xcorner+0.016    0.82    0.3500    0.0987],'FontSize',TextFontSize)
set(h2,'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth],'XTickLabel',[], ...
   'FontSize',AnnotFontSize,'YTick',[0:20:100],'YTickLabel',[],'TickLength',[TICKlength 0.025])
set(h3,'Position',[Xcorner+2*Xwidth+Xshift Ycorner Xwidth Ywidth],'XTickLabel',[], ...
   'FontSize',AnnotFontSize,'YTick',[0:20:100],'YTickLabel',[],'TickLength',[TICKlength 0.025])

set(h4,'Position',[Xcorner Ycorner-(Ywidth+Yshift) Xwidth Ywidth],'XTickLabel',[0.1 1 10], ...\
   'YTickLabel',[1 10],'FontSize',AnnotFontSize,'TickLength',[TICKlength 0.025])
set(h5,'Position',[Xcorner+Xwidth+Xshift Ycorner-(Ywidth+Yshift) Xwidth Ywidth],'XTickLabel',[0.1 1 10], ...
   'FontSize',AnnotFontSize,'YTickLabel',[],'TickLength',[TICKlength 0.025])
set(h6,'Position',[Xcorner+2*Xwidth+Xshift Ycorner-(Ywidth+Yshift) Xwidth Ywidth],'XTickLabel',[0.1 1 10], ...
   'FontSize',AnnotFontSize,'YTickLabel',[],'TickLength',[TICKlength 0.025])

drawnow
hold off

cd(CURdir)

save 'Fig4_workspace'

% eval(['cd ''' fullfile(RLFpaper_dir,'Mfiles','EPSfigs') ''''])
% orient landscape
print -dtiff ThrQ10

% print -depsc fig_ThrQ10

