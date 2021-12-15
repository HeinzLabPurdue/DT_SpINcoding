clear;
clc;

saveFig= 0;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];
use_trifilt= 1;
npoints_trifilt= 11;

[sig, fsOrg]= audioread(fullfile(pwd, 'stimuli', 'SNR_0', 'FLN_Stim_S_P.wav'));
fs= 24e3;

% rep:         d      e    n         g    a    m    l    e        m    a    n    d         s    m    i    l     e     d      e         s
tVals_ortho= [.012, .072, .12, nan, .2, .27, .35, .42, .46, nan, .52, .59, .64, .69, nan, .79, .9, .97, 1.02, 1.06, 1.12, 1.18, nan, 1.256];
orthographic_rep= 'den gamle mand smilede s';

sig= helper.gen_resample(sig, fsOrg, fs);
t= (1:length(sig))/fs;


figSize_cm= [55 5 8.5 7];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

%% F0, Formants, and level
danish.voiced_boundaries= helper.find_voicing_boundaries(sig, fs, 0, .13);
tFormants= 0.0:.005:max(t);

temp_f0= load(['mat_data' filesep 'danish_pitch.mat']);
temp_f0= temp_f0.pitch_data;
danish.voiced_inds= any(tFormants>danish.voiced_boundaries(:,1) & tFormants<danish.voiced_boundaries(:,2), 1);
danish.trajectory.f0= zeros(size(tFormants));
danish.trajectory.f0(danish.voiced_inds)= interp1([temp_f0.time], [temp_f0.est], tFormants(danish.voiced_inds), 'pchip');
danish.voiced_inds= danish.voiced_inds & danish.trajectory.f0>95 & danish.trajectory.f0<150; % make talker - looking at estiamtes, This removes the edge estimates.


temp_formants= load(['mat_data' filesep 'danish_formant.mat']);
temp_formants= temp_formants.formant_data;
danish.trajectory.f1= interp1([temp_formants.time], [temp_formants.f1], tFormants, 'pchip');
danish.trajectory.f2= interp1([temp_formants.time], [temp_formants.f2], tFormants, 'pchip');
danish.trajectory.f3= interp1([temp_formants.time], [temp_formants.f3], tFormants, 'pchip');

danish.trajectory.f0(~danish.voiced_inds)= nan;
danish.trajectory.f1(~danish.voiced_inds)= nan;
danish.trajectory.f2(~danish.voiced_inds)= nan;
danish.trajectory.f3(~danish.voiced_inds)= nan;

if use_trifilt
    danish.trajectory.f0= helper.trifilt(danish.trajectory.f0, npoints_trifilt);
    danish.trajectory.f1= helper.trifilt(danish.trajectory.f1, npoints_trifilt);
    danish.trajectory.f2= helper.trifilt(danish.trajectory.f2, npoints_trifilt);
    danish.trajectory.f3= helper.trifilt(danish.trajectory.f3, npoints_trifilt);
end



[stim_level, t_level]= helper.gen_get_spl_vals(sig, fs, 40e-3, 0.5);
%%

plt.lw= 0.75;
plt.lw2= 1.0;
plt.lw3= 1.5;
plt.tWindow= 50e-3;
plt.fracOVlap= .95;
plt.PSD_Range= 50;
plt.tick_len= [.03 .02];

%%

sp_ax(1)= subplot(211);
plot(t, sig, 'linew', plt.lw);
ylabel('Amplitude (a.u.)');
ylim([-.29 .24]);
set(gca, 'XTickLabel', '');

txtHan= nan(length(tVals_ortho), 1);
for letVar=1:length(tVals_ortho)
    txtHan(letVar)= text(tVals_ortho(letVar), .2, orthographic_rep(letVar));
end

yyaxis right;
plot(t_level, stim_level, 'k', 'linew', plt.lw2);
set(gca, 'YColor', 'k');
ylabel('Level (dB SPL)');
box off;
ttlHan(1)= text(.0, 1.12, 'A', 'Units', 'normalized');

sp_ax(2)= subplot(212);
hold on;
helper.plot_spectrogram(sig, fs, plt.tWindow, plt.fracOVlap);
ylim([0 8]);
ttlHan(2)= text(.0, 1.06, 'B', 'Units', 'normalized');

ColBar_Spec= max(get(colorbar, 'Limits'));
caxis([ColBar_Spec-plt.PSD_Range ColBar_Spec]);
colorbar off;
%set(gca, 'TickDir', 'Out');
box off;

line(tFormants, danish.trajectory.f1/1e3, 'linew', plt.lw2, 'color', 'r');
line(tFormants, danish.trajectory.f2/1e3, 'linew', plt.lw2, 'color', 'r');
line(tFormants, danish.trajectory.f3/1e3, 'linew', plt.lw2, 'color', 'r');

yyaxis right;
line(tFormants, danish.trajectory.f0, 'linew', plt.lw2, 'color', 'k');
linkaxes(sp_ax, 'x')
xlim([-.005 1.305]);
ylabel('F0 (Hz)');
set(gca, 'YColor', 'k');

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(txtHan,'FontSize', 10);
set(ttlHan,'FontSize', 11, 'FontWeight', 'bold');
set(findall(gcf,'-property','TickLength'),'TickLength', [.02 .01]);

%%
Xcorner= .135;
Xwidth= .725;
Ycorner= .145;
Ywidth1= .26;
Ywidth2= .48;
Yshift= .065;

set(sp_ax(2),'Position',[Xcorner Ycorner Xwidth Ywidth2])
drawnow

set(sp_ax(1),'Position',[Xcorner Ycorner+Yshift+Ywidth2 Xwidth Ywidth1])
drawnow

%%
if saveFig
    fName= [dirStruct.png 'FigS1_stim_spectrogram_ortho'];
    print(fName, '-dpng',  '-r600');
    saveas(gcf, [dirStruct.eps 'FigS1_stim_spectrogram_ortho'], 'epsc');
end