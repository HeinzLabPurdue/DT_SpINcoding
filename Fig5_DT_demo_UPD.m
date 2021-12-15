clear;
clc;

saveFig= 1;
figSize_cm= [10 5 11.6 7];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];
dirStruct.loading_dir= [pwd filesep 'ANF_Data' filesep];

all_chinID= [379 360];
allChinSpikeData = [];
if ~exist('all_chinID', 'var')
    allfiles=dir([dirStruct.loading_dir '*.mat']);
else
    for fileVar=1:length(all_chinID)
        curChinID= all_chinID(fileVar);
        allfiles(fileVar)=dir([dirStruct.loading_dir '*' num2str(curChinID) '*']); %#ok<SAGROW>
    end
end
for fileVar=1:length(allfiles)
    temp = load([dirStruct.loading_dir allfiles(fileVar).name]);
    allChinSpikeData= [allChinSpikeData; temp.spike_data']; %#ok<AGROW>
end
cur_snr_allChinSpikeData= allChinSpikeData(strcmp({allChinSpikeData.noise}, 'SSN'));
[~, uniqInds]= unique([ [cur_snr_allChinSpikeData.chinID]', [cur_snr_allChinSpikeData.track]', [cur_snr_allChinSpikeData.unit]', [cur_snr_allChinSpikeData.SPL]'], 'rows');
uniqSpikeData= cur_snr_allChinSpikeData(uniqInds);

[sig, fsOrg]= audioread(['mat_data' filesep 'FLN_Stim_S_P.wav']);

fs= 20e3;
sig= helper.gen_resample(sig, fsOrg, fs);
sig= helper.gen_rescale(sig, 65);
tSig_ms= (1:length(sig))/fs*1e3;

anl.tStart_ms= 235;
anl.tEnd_ms= 325;
anl.tZoom_ms= 10;

anl.tBinEdges_s= (anl.tStart_ms-anl.tZoom_ms)/1e3:1/fs:(anl.tEnd_ms+anl.tZoom_ms)/1e3;
anl.tBinCenters_ms= (( anl.tBinEdges_s(1:end-1) + anl.tBinEdges_s(2:end) )/2)*1e3;
anl.validWindowInds= (anl.tBinCenters_ms>anl.tStart_ms) & (anl.tBinCenters_ms<anl.tEnd_ms);

% Other potential example: [235-325 ms] [HI: 360 4 5] [NH: 379 2 7]
% NH example
nh.chin_track_unit= [379, 2, 7];
nh.ind= ismember([[uniqSpikeData.chinID]', [uniqSpikeData.track]', [uniqSpikeData.unit]'], nh.chin_track_unit, 'rows');
nh.data= uniqSpikeData(nh.ind);
nh.uR_pos_plot= histcounts( cell2mat(nh.data.SpikeTrains{1,1}), anl.tBinEdges_s) / length(nh.data.SpikeTrains{1,1}) * fs;
nh.uR_neg_plot= histcounts( cell2mat(nh.data.SpikeTrains{1,2}), anl.tBinEdges_s) / length(nh.data.SpikeTrains{1,2}) * fs;
nh.uR_pos_window= nh.uR_pos_plot(anl.validWindowInds);
nh.uR_neg_window= nh.uR_neg_plot(anl.validWindowInds);
nh.dT_window= (nh.uR_pos_window-nh.uR_neg_window)/2;
nh.tipThresh= nh.data.thresh_dB;
nh.tailThresh= nh.data.thresh_dB + nh.data.TTR_dB;
nh.tailFreq= nh.data.TC.freqkHz(dsearchn(nh.data.TC.TCfit(:), nh.tailThresh));

% HI example
hi.chin_track_unit= [360, 4, 6];
hi.ind= ismember([[uniqSpikeData.chinID]', [uniqSpikeData.track]', [uniqSpikeData.unit]'], hi.chin_track_unit, 'rows');
hi.data= uniqSpikeData(hi.ind);
hi.uR_pos_plot= histcounts( cell2mat(hi.data.SpikeTrains{1,1}), anl.tBinEdges_s) / length(hi.data.SpikeTrains{1,1}) * fs;
hi.uR_neg_plot= histcounts( cell2mat(hi.data.SpikeTrains{1,2}), anl.tBinEdges_s) / length(hi.data.SpikeTrains{1,2}) * fs;
hi.uR_pos_window= hi.uR_pos_plot(anl.validWindowInds);
hi.uR_neg_window= hi.uR_neg_plot(anl.validWindowInds);
hi.dT_window= (hi.uR_pos_window-hi.uR_neg_window)/2;
hi.tipThresh= hi.data.thresh_dB;
hi.tailThresh= hi.data.thresh_dB + hi.data.TTR_dB;
hi.tailFreq= hi.data.TC.freqkHz(dsearchn(hi.data.TC.TCfit(:), hi.tailThresh));

% Signal
tSig_plot= tSig_ms( (tSig_ms > (anl.tStart_ms-anl.tZoom_ms)) & (tSig_ms < (anl.tEnd_ms+anl.tZoom_ms)) );
tSig_window= tSig_ms( (tSig_ms > anl.tStart_ms) & (tSig_ms < anl.tEnd_ms) );

sig_plot= sig( (tSig_ms > (anl.tStart_ms-anl.tZoom_ms)) & (tSig_ms < (anl.tEnd_ms+anl.tZoom_ms)) );
sig_window_nh= sig( (tSig_ms > anl.tStart_ms) & (tSig_ms < anl.tEnd_ms) );

sig= helper.gen_rescale(sig, 80);
sig_window_hi= sig( (tSig_ms > anl.tStart_ms) & (tSig_ms < anl.tEnd_ms) );


%% Plot

plt.AmpStim= max([max(nh.uR_pos_plot), max(nh.uR_neg_plot), max(hi.uR_pos_plot), max(hi.uR_neg_plot)]);
plt.ShiftSig= 0.9*plt.AmpStim;
plt.ScaleSig= 0.33*plt.AmpStim/max(abs(sig_plot));
plt.Shift_uR= 1.33*plt.AmpStim;
plt.lwCF= 1.0;
plt.lw1= 1.25;
plt.lwTC= 1.5;
plt.lw2= 2;
plt.mrkSize0= 3;
plt.mrkSize= 4;
plt.mrkSizeCF= 6;
plt.mrkSize2= 8;
plt.freqTick= [.1 .5 1 4];
plt.tick_len= [.03 .03];
plt.fSize= 9;
plt.ttl_fSize= 11;
plt.TCylim= [-11 90];
plt.PSDxlim= [.07 4.1];
plt.TCthreshTick= [0 40 80];

figure(1);
clf;
sp_ax(1)= subplot(211);
hold on;

plot(tSig_plot, -plt.ShiftSig + plt.ScaleSig*sig_plot, 'Color', helper.get_color('wg'));
lineHan(3)= plot(tSig_window, -plt.ShiftSig + plt.ScaleSig*sig_window_nh, 'Color', helper.get_color('k'));

plot(anl.tBinCenters_ms, plt.Shift_uR + nh.uR_pos_plot, 'Color', helper.get_color('wg'));
plot(anl.tBinCenters_ms, plt.Shift_uR - nh.uR_neg_plot, 'Color', helper.get_color('wg'));
lineHan(1)= plot(anl.tBinCenters_ms(anl.validWindowInds), plt.Shift_uR + nh.uR_pos_window, 'Color', helper.get_color('b'));
plot(anl.tBinCenters_ms(anl.validWindowInds), plt.Shift_uR - nh.uR_neg_window, 'Color', helper.get_color('lb'));

plot(anl.tBinCenters_ms, hi.uR_pos_plot, 'Color', helper.get_color('wg'));
plot(anl.tBinCenters_ms, -hi.uR_neg_plot, 'Color', helper.get_color('wg'));
lineHan(2)= plot(anl.tBinCenters_ms(anl.validWindowInds), hi.uR_pos_window, 'Color', helper.get_color('r'));
plot(anl.tBinCenters_ms(anl.validWindowInds), -hi.uR_neg_window, 'Color', helper.get_color('lr'));

yTickScale= .1*fs;
ytickValsDefault= [0 yTickScale];
ytickVals= [ytickValsDefault plt.Shift_uR+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault/1e3, 1, 2)), 'UniformOutput', false);
time_tick_vals= [anl.tStart_ms anl.tEnd_ms];

xlim([anl.tStart_ms-anl.tZoom_ms anl.tEnd_ms+2.5*anl.tZoom_ms]);
xLabHan= xlabel('Time (ms)', 'Units', 'normalized');
xLabHan.Position(2)= -.1;

ylabel('Rate (spikes/s)')
txtHan(1)= text(.05, .99, 'A', 'Units', 'normalized', 'FontWeight', 'bold');
set(gca, 'ytick', ytickVals, 'YTickLabel', ytickLabs, 'XTick', time_tick_vals);
ylim([-.34 0.6]*fs);
[legHan, icons]= legend(lineHan, 'NH', 'HI', 'Stim.', 'location', 'east', 'box', 'off');
count= 4;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
legHan.Position(1:2)= [.83 .65];
text(-.06, 1.01, 'x 10^3', 'Units', 'normalized');


sp_ax(2)= subplot(234);
yyaxis left;
hold on;
[~,~,lHan_sigPSD]= helper.plot_dpss_psd(sig_window_nh, fs, 'xunit', 'khz', 'yrange', 50);
set(lHan_sigPSD, 'Color', 1.2*helper.get_color('gray'), 'LineWidth', plt.lw1);
set(gca, 'YColor', 'k');
text(.1, -48, 'F_0', 'HorizontalAlignment', 'center');
text(0.7, -46, 'F_1', 'HorizontalAlignment', 'center');
text(1.3, -53, 'F_2', 'HorizontalAlignment', 'center');
ylim([-80 -25]);

yyaxis right;
plot(nh.data.TC.freqkHz, nh.data.TC.TCfit, ':', 'Color', 'b', 'LineWidth', plt.lwTC);
plot(hi.data.TC.freqkHz, hi.data.TC.TCfit, ':', 'Color', 'r', 'LineWidth', plt.lwTC);
ylim([-1 85])
xlim(plt.PSDxlim);
xlab_han= xlabel('Frequency (kHz)', 'Units', 'normalized');
xlab_han.Position(1:2)= [1.15, -0.14];
set(gca, 'XTick', plt.freqTick, 'TickLength', plt.tick_len, 'YColor', 'k', 'YTick', plt.TCthreshTick);
txtHan(2)= text(.05, .99, 'B', 'Units', 'normalized', 'FontWeight', 'bold');
plot([nh.tailFreq nh.tailFreq], [nh.tailThresh nh.tipThresh], '-', 'LineWidth', plt.lw1, 'Color', helper.get_color('b'));
plot(nh.tailFreq, nh.tailThresh-1, '^', 'LineWidth', plt.lw1, 'Color', 'b', 'MarkerSize', plt.mrkSize0);
plot(nh.tailFreq, nh.tipThresh+1, 'v', 'LineWidth', plt.lw1, 'Color', 'b', 'MarkerSize', plt.mrkSize0);
plot(nh.data.CF_Hz/1e3, nh.tipThresh-1, 'p', 'LineWidth', plt.lwCF, 'Color', 'b', 'MarkerSize', plt.mrkSizeCF);

plot([hi.tailFreq hi.tailFreq], [hi.tailThresh hi.tipThresh], '-', 'LineWidth', plt.lw1, 'Color', helper.get_color('r'));
plot(hi.tailFreq, hi.tailThresh-1, '^', 'LineWidth', plt.lw1, 'Color', 'r', 'MarkerSize', plt.mrkSize0);
plot(hi.tailFreq, hi.tipThresh+1, 'v', 'LineWidth', plt.lw1, 'Color', 'r', 'MarkerSize', plt.mrkSize0);
plot(hi.data.CF_Hz/1e3, hi.tipThresh-1, 'p', 'LineWidth', plt.lwCF, 'Color', 'r', 'MarkerSize', plt.mrkSizeCF);
ylabel('TC thresh. (dB SPL)');

sp_ax(3)= subplot(235);
% yyaxis left;
hold on;
[Pxx_nh_dB,freq_nh,lHan_dT_PSD_nh]= helper.plot_dpss_psd(nh.dT_window/rms(nh.dT_window), fs, 'xunit', 'khz');
set(lHan_dT_PSD_nh, 'Color', helper.get_color('b'), 'LineWidth', plt.lw1);
[~, max_ind_nh]= max(Pxx_nh_dB);
arrowHan(1)= text(freq_nh(max_ind_nh), Pxx_nh_dB(max_ind_nh)+3.5, '\downarrow', 'Color', helper.get_color('b'), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
arrowTxtHan(1)= text(freq_nh(max_ind_nh), Pxx_nh_dB(max_ind_nh)+5.5, 'F_2', 'Color', 'k', 'HorizontalAlignment', 'center');

[Pxx_hi_dB,freq_hi,lHan_dT_PSD_hi]= helper.plot_dpss_psd(hi.dT_window/rms(hi.dT_window), fs, 'xunit', 'khz', 'yrange', 30);
set(lHan_dT_PSD_hi, 'Color', helper.get_color('lr'), 'LineWidth', plt.lw1);
[~, max_ind_hi]= max(Pxx_hi_dB);
arrowHan(2)= text(freq_hi(max_ind_hi), Pxx_hi_dB(max_ind_hi)+3.5, '\downarrow', 'Color', helper.get_color('r'), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
arrowTxtHan(2)= text(freq_hi(max_ind_hi), Pxx_hi_dB(max_ind_hi)+5.5, 'F_1', 'Color', helper.get_color('k'), 'HorizontalAlignment', 'center');
set(gca, 'YColor', 'k', 'XTick', plt.freqTick);
ylabel('');
ylim([-50 -22]);
xlim(plt.PSDxlim);
xlabel('');

txtHan(3)= text(.05, .99, 'C', 'Units', 'normalized', 'FontWeight', 'bold');

set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);
set(txtHan, 'fontsize', plt.ttl_fSize);
set(arrowHan, 'fontsize', 14);


%% Load AN data
anl.NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
anl.HIchins= [358 360 361 362 367 370];
allChinData= helper.load_danish_chin_data(anl);

cf_kHz= [allChinData.CF_Hz]/1e3;
TTR_dB= [allChinData.TTR_dB];
nhInds= ismember([allChinData.chinID], anl.NHchins);
hiInds= ismember([allChinData.chinID], anl.HIchins);
CFcenters= [.5 1 2 4 8];


sp_ax(4)= subplot(236);
hold on;
lHan(1)= plot(cf_kHz(nhInds), TTR_dB(nhInds), 'x', 'color', helper.get_color('b'), 'MarkerSize', plt.mrkSize);
lHan(2)= plot(cf_kHz(hiInds), TTR_dB(hiInds), '+', 'color', helper.get_color('r'), 'MarkerSize', plt.mrkSize);
set(gca, 'xscale', 'log', 'xtick', CFcenters, 'TickLength', plt.tick_len, 'YTick', plt.TCthreshTick);
xlab_hanD= xlabel('CF (kHz)', 'Units', 'normalized');
xlab_hanD.Position(2)= xlab_han.Position(2);
ylabel('TTR (dB)');
txtHan(4)= text(.05, .99, 'D', 'Units', 'normalized', 'FontWeight', 'bold');
[legHan, icons]= legend(lHan, 'NH', 'HI', 'Location', 'northwest', 'box', 'off');
icons(1).Position(1)= icons(1).Position(1)-0.2;
icons(2).Position(1)= icons(2).Position(1)-0.2;

%% define new axes for AB
Xcorner_X= .09;
Xwidth_X= .225;
Xshift_X= .09;
Ycorner_X= .116;
Ywidth_X= .355;
Yshift_X= .13;

set(sp_ax(1),'Position',[Xcorner_X Ycorner_X+Ywidth_X+Yshift_X 3*Xwidth_X+2.4*Xshift_X 1.03*Ywidth_X])
drawnow

set(sp_ax(2),'Position',[Xcorner_X Ycorner_X Xwidth_X Ywidth_X])
drawnow

set(sp_ax(3),'Position',[Xcorner_X+Xwidth_X+1.5*Xshift_X Ycorner_X Xwidth_X Ywidth_X])
drawnow

set(sp_ax(4),'Position',[Xcorner_X+2*Xwidth_X+2.4*Xshift_X Ycorner_X Xwidth_X Ywidth_X])
drawnow


if saveFig
    print([dirStruct.png 'Fig3_dt_demo_TTR'], '-dpng',  '-r600');
    saveas(gcf, [dirStruct.eps 'Fig3_dt_demo_TTR'], 'epsc');
end
