clear;
clc;

saveFig= 1;
dirStruct.png= [pwd filesep 'final_figs' filesep];

nFibers= 30;
hearing_loss_dB= 30;
low_int_level= 25;
high_int_level= 60;


nh_thresh_min= 0;
nh_thresh_range= 25-.5;
nh_thresh_max= nh_thresh_min+nh_thresh_range;

hi_comp_thresh_min= nh_thresh_min+hearing_loss_dB;
hi_comp_thresh_range= 15;
hi_comp_thresh_max= hi_comp_thresh_min+hi_comp_thresh_range;

hi_expand_thresh_min= nh_thresh_min+hearing_loss_dB;
hi_expand_thresh_range= 50;
hi_expand_thresh_max= hi_expand_thresh_min+hi_expand_thresh_range;

gain_comp_group= hi_comp_thresh_min - nh_thresh_min;
gain_expand_group= hi_expand_thresh_min - nh_thresh_min;


nh_active_max_low= min(nh_thresh_max, low_int_level);
hi_comp_active_max_low= min(hi_comp_thresh_max, low_int_level+gain_comp_group);
hi_expand_active_max_low= min(hi_expand_thresh_max, low_int_level+gain_expand_group);

nh_active_max_high= min(nh_thresh_max, high_int_level);
hi_comp_active_max_high= min(hi_comp_thresh_max, high_int_level+gain_comp_group);
hi_expand_active_max_high= min(hi_expand_thresh_max, high_int_level+gain_expand_group);

%%
plt.lw= 1.5;
plt.fAlpha= .15;
plt.stim_level_Y= -.02;
plt.mulX= -0.01;
plt.mulY= 1.08;
plt.txtX= -0.2;
plt.txtY= 0.25;
plt.ttlY= 1.1;


figSize_in= [1 1 7 4];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'inches', figSize_in, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
set(gcf,figure_prop_name,figure_prop_val);


% Threshold distribution 

sp_ax(1)= subplot(2, 2, 1);
hold on;
plot([nh_thresh_min-.01 nh_thresh_min nh_thresh_max nh_thresh_max+.01], nFibers*[0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'b', 'linew', plt.lw);
plot([hi_comp_thresh_min-.01 hi_comp_thresh_min hi_comp_thresh_max hi_comp_thresh_max+.01], nFibers*[0 1 1 0]/(hi_comp_thresh_max-hi_comp_thresh_min), 'r', 'linew', plt.lw);

lHan(1)= fill([nh_thresh_min nh_active_max_low nh_active_max_low nh_thresh_min  nh_thresh_min], nFibers*[0 0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'b', 'facealpha', plt.fAlpha, 'LineStyle','none');
lHan(2)= fill([hi_comp_thresh_min hi_comp_active_max_low hi_comp_active_max_low hi_comp_thresh_min  hi_comp_thresh_min], nFibers*[0 0 1 1 0]/(hi_comp_thresh_max-hi_comp_thresh_min), 'r', 'facealpha', plt.fAlpha, 'LineStyle','none');

plot(low_int_level, plt.stim_level_Y, '^b');
plot(low_int_level+gain_comp_group, plt.stim_level_Y, '^r');
text(plt.txtX, plt.txtY, 'Compressed', 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Rotation', 90);
ttlHan= title('Low-intensity stimuli', 'Units', 'normalized');
ttlHan.Position(2)= plt.ttlY;
legend(lHan, 'NH', 'HI', 'box', 'off');
set(gca, 'XTickLabel', '');
text(plt.mulX, plt.mulY, 'x10^3', 'Units', 'normalized');

sp_ax(3)= subplot(2, 2, 3);
hold on;
plot([nh_thresh_min-.01 nh_thresh_min nh_thresh_max nh_thresh_max+.01], nFibers*[0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'linew', plt.lw);
plot([hi_expand_thresh_min-.01 hi_expand_thresh_min hi_expand_thresh_max hi_expand_thresh_max+.01], nFibers*[0 1 1 0]/(hi_expand_thresh_max-hi_expand_thresh_min), 'linew', plt.lw);

fill([nh_thresh_min nh_active_max_low nh_active_max_low nh_thresh_min  nh_thresh_min], nFibers*[0 0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'b', 'facealpha', plt.fAlpha, 'LineStyle','none');
fill([hi_expand_thresh_min hi_expand_active_max_low hi_expand_active_max_low hi_expand_thresh_min  hi_expand_thresh_min], ...
    nFibers*[0 0 1 1 0]/(hi_expand_thresh_max-hi_expand_thresh_min), 'r', 'facealpha', plt.fAlpha, 'LineStyle','none');

plot(low_int_level, plt.stim_level_Y, '^b');
plot(low_int_level+gain_expand_group, plt.stim_level_Y, '^r');
text(plt.txtX, plt.txtY, 'Expanded', 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Rotation', 90);
text(plt.mulX, plt.mulY, 'x10^3', 'Units', 'normalized');

xlabHan= xlabel('Level (dB SPL)', 'units', 'normalized');
xlabHan.Position(1)= 1.15;
ylabHan= ylabel('#AN Fibers', 'Units', 'normalized');
ylabHan.Position(1:2)= [-.1 1.2];

% Coding of low-intensity stimuli


sp_ax(2)= subplot(2, 2, 2);
hold on;
plot([nh_thresh_min-.01 nh_thresh_min nh_thresh_max nh_thresh_max+.01], nFibers*[0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'linew', plt.lw);
plot([hi_comp_thresh_min-.01 hi_comp_thresh_min hi_comp_thresh_max hi_comp_thresh_max+.01], nFibers*[0 1 1 0]/(hi_comp_thresh_max-hi_comp_thresh_min), 'linew', plt.lw);

fill([nh_thresh_min nh_active_max_high nh_active_max_high nh_thresh_min  nh_thresh_min], nFibers*[0 0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'b', 'facealpha', plt.fAlpha, 'LineStyle','none');
fill([hi_comp_thresh_min hi_comp_active_max_high hi_comp_active_max_high hi_comp_thresh_min  hi_comp_thresh_min], ...
    nFibers*[0 0 1 1 0]/(hi_comp_thresh_max-hi_comp_thresh_min), 'r', 'facealpha', plt.fAlpha, 'LineStyle','none');

plot(high_int_level, plt.stim_level_Y, '^b');
plot(high_int_level+gain_comp_group, plt.stim_level_Y, '^r');
ttlHan= title('High-intensity stimuli', 'Units', 'normalized');
ttlHan.Position(2)= plt.ttlY;
set(gca, 'XTickLabel', '');
text(plt.mulX, plt.mulY, 'x10^3', 'Units', 'normalized');

sp_ax(4)= subplot(2, 2, 4);
hold on;
plot([nh_thresh_min-.01 nh_thresh_min nh_thresh_max nh_thresh_max+.01], nFibers*[0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'linew', plt.lw);
plot([hi_expand_thresh_min-.01 hi_expand_thresh_min hi_expand_thresh_max hi_expand_thresh_max+.01], nFibers*[0 1 1 0]/(hi_expand_thresh_max-hi_expand_thresh_min), 'linew', plt.lw);

fill([nh_thresh_min nh_active_max_high nh_active_max_high nh_thresh_min  nh_thresh_min], nFibers*[0 0 1 1 0]/(nh_thresh_max-nh_thresh_min), 'b', 'facealpha', plt.fAlpha, 'LineStyle','none');
fill([hi_expand_thresh_min hi_expand_active_max_high hi_expand_active_max_high hi_expand_thresh_min  hi_expand_thresh_min], ...
    nFibers*[0 0 1 1 0]/(hi_expand_thresh_max-hi_expand_thresh_min), 'r', 'facealpha', plt.fAlpha, 'LineStyle','none');

plot(high_int_level, plt.stim_level_Y, '^b');
plot(high_int_level+gain_expand_group, plt.stim_level_Y, '^r');
text(plt.mulX, plt.mulY, 'x10^3', 'Units', 'normalized');

linkaxes(sp_ax);
set(findall(gcf,'-property','FontSize'),'FontSize', 12);
set(findall(gcf,'-property','TickLength'),'TickLength', [.02 .03]);
set(findall(gcf,'-property','XTick'),'XTick', 0:20:100);

%% define new axes for AB
Xcorner_AB= .09;
Xwidth_AB= .36;
Xshift_horz= .12;

Ycorner_AB= .11;
Ywidth_AB= .35;
Yshift_AB= .09;

% B
set(sp_ax(3),'Position',[Xcorner_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% A
set(sp_ax(1),'Position',[Xcorner_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow
Xcorner_CD= Xcorner_AB+Xwidth_AB+Xshift_horz;

% D
set(sp_ax(4),'Position',[Xcorner_CD Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% C
set(sp_ax(2),'Position',[Xcorner_CD Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow

%%
if saveFig
    print([dirStruct.png 'FigS2_vowel_conso_hypo'], '-dpng',  '-r600');
end
