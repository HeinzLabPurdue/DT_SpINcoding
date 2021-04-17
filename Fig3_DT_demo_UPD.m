clear;
clc;

saveFig= 1;
figSize_cm= [10 5 17 12];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

dirStruct.png= [pwd filesep 'final_figs' filesep];
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
sig= gen_resample(sig, fsOrg, fs);
sig= gen_rescale(sig, 65);
tSig= (1:length(sig))/fs;

anl.tStart= 235e-3;
anl.tEnd= 325e-3;
anl.tZoom= 10e-3;

anl.tBinEdges= (anl.tStart-anl.tZoom):1/fs:(anl.tEnd+anl.tZoom);
anl.tBinCenters= ( anl.tBinEdges(1:end-1) + anl.tBinEdges(2:end) )/2;
anl.validWindowInds= (anl.tBinCenters>anl.tStart) & (anl.tBinCenters<anl.tEnd);

% Other potential example: [235-325 ms] [HI: 360 4 5] [NH: 379 2 7]
% NH example
nh.chin_track_unit= [379, 2, 7];
nh.ind= ismember([[uniqSpikeData.chinID]', [uniqSpikeData.track]', [uniqSpikeData.unit]'], nh.chin_track_unit, 'rows');
nh.data= uniqSpikeData(nh.ind);
nh.uR_pos_plot= histcounts( cell2mat(nh.data.SpikeTrains{1,1}), anl.tBinEdges) / length(nh.data.SpikeTrains{1,1}) * fs;
nh.uR_neg_plot= histcounts( cell2mat(nh.data.SpikeTrains{1,2}), anl.tBinEdges) / length(nh.data.SpikeTrains{1,2}) * fs;
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
hi.uR_pos_plot= histcounts( cell2mat(hi.data.SpikeTrains{1,1}), anl.tBinEdges) / length(hi.data.SpikeTrains{1,1}) * fs;
hi.uR_neg_plot= histcounts( cell2mat(hi.data.SpikeTrains{1,2}), anl.tBinEdges) / length(hi.data.SpikeTrains{1,2}) * fs;
hi.uR_pos_window= hi.uR_pos_plot(anl.validWindowInds);
hi.uR_neg_window= hi.uR_neg_plot(anl.validWindowInds);
hi.dT_window= (hi.uR_pos_window-hi.uR_neg_window)/2;
hi.tipThresh= hi.data.thresh_dB;
hi.tailThresh= hi.data.thresh_dB + hi.data.TTR_dB;
hi.tailFreq= hi.data.TC.freqkHz(dsearchn(hi.data.TC.TCfit(:), hi.tailThresh));

% Signal
tSig_plot= tSig( (tSig > (anl.tStart-anl.tZoom)) & (tSig < (anl.tEnd+anl.tZoom)) );
tSig_window= tSig( (tSig > anl.tStart) & (tSig < anl.tEnd) );

sig_plot= sig( (tSig > (anl.tStart-anl.tZoom)) & (tSig < (anl.tEnd+anl.tZoom)) );
sig_window_nh= sig( (tSig > anl.tStart) & (tSig < anl.tEnd) );

sig= gen_rescale(sig, 80);
sig_window_hi= sig( (tSig > anl.tStart) & (tSig < anl.tEnd) );


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
plt.TCylim= [-11 95];
plt.PSDxlim= [.07 4.1];

figure(1);
clf;
sp_ax(1)= subplot(211);
hold on;

plot(tSig_plot, -plt.ShiftSig + plt.ScaleSig*sig_plot, 'Color', get_color('wg'));
lineHan(3)= plot(tSig_window, -plt.ShiftSig + plt.ScaleSig*sig_window_nh, 'Color', get_color('k'));

plot(anl.tBinCenters, plt.Shift_uR + nh.uR_pos_plot, 'Color', get_color('wg'));
plot(anl.tBinCenters, plt.Shift_uR - nh.uR_neg_plot, 'Color', get_color('wg'));
lineHan(1)= plot(anl.tBinCenters(anl.validWindowInds), plt.Shift_uR + nh.uR_pos_window, 'Color', get_color('b'));
plot(anl.tBinCenters(anl.validWindowInds), plt.Shift_uR - nh.uR_neg_window, 'Color', get_color('lb'));

plot(anl.tBinCenters, hi.uR_pos_plot, 'Color', get_color('wg'));
plot(anl.tBinCenters, -hi.uR_neg_plot, 'Color', get_color('wg'));
lineHan(2)= plot(anl.tBinCenters(anl.validWindowInds), hi.uR_pos_window, 'Color', get_color('r'));
plot(anl.tBinCenters(anl.validWindowInds), -hi.uR_neg_window, 'Color', get_color('lr'));

yTickScale= .1*fs;
ytickValsDefault= [-yTickScale 0 yTickScale];
ytickVals= [ytickValsDefault plt.Shift_uR+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault/1e3, 1, 2)), 'UniformOutput', false);

xlim([anl.tStart-anl.tZoom anl.tEnd+2*anl.tZoom]);
xlabel('Time (s)');
ylabel('Rate (spikes/s)')
txtHan(1)= text(.05, .99, 'A', 'Units', 'normalized', 'FontWeight', 'bold');
set(gca, 'ytick', ytickVals, 'YTickLabel', ytickLabs);
ylim([-.4 0.6]*fs)
[legHan, icons]= legend(lineHan, 'NH', 'HI', 'Stim.', 'location', 'east', 'box', 'off');
count= 4;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
legHan.Position(1:2)= [.87 .62];
text(-.05, 1.01, 'x 10^3', 'Units', 'normalized');


sp_ax(2)= subplot(234);
yyaxis left;
hold on;
[~,~,lHan_sigPSD]= plot_dpss_psd(sig_window_nh, fs, 'xunit', 'khz', 'yrange', 50);
set(lHan_sigPSD, 'Color', 1.2*get_color('gray'), 'LineWidth', plt.lw1);
set(gca, 'YColor', 'k');
text(.1, -50, 'F_0', 'HorizontalAlignment', 'center');
text(0.7, -47, 'F_1', 'HorizontalAlignment', 'center');
text(1.3, -54, 'F_2', 'HorizontalAlignment', 'center');
ylim([-80 -30]);

yyaxis right;
plot(nh.data.TC.freqkHz, nh.data.TC.TCfit, ':', 'Color', 'b', 'LineWidth', plt.lwTC);
plot(hi.data.TC.freqkHz, hi.data.TC.TCfit, ':', 'Color', 'r', 'LineWidth', plt.lwTC);
ylim(plt.TCylim);
xlim(plt.PSDxlim);
xlab_han= xlabel('Frequency (kHz)', 'Units', 'normalized');
xlab_han.Position(1)= 1.12;
set(gca, 'XTick', plt.freqTick, 'TickLength', plt.tick_len, 'YColor', 'k');
txtHan(2)= text(.05, .99, 'B', 'Units', 'normalized', 'FontWeight', 'bold');
plot([nh.tailFreq nh.tailFreq], [nh.tailThresh nh.tipThresh], '-', 'LineWidth', plt.lw1, 'Color', get_color('b'));
plot(nh.tailFreq, nh.tailThresh-1, '^', 'LineWidth', plt.lw1, 'Color', 'b', 'MarkerSize', plt.mrkSize0);
plot(nh.tailFreq, nh.tipThresh+1, 'v', 'LineWidth', plt.lw1, 'Color', 'b', 'MarkerSize', plt.mrkSize0);
plot(nh.data.CF_Hz/1e3, nh.tipThresh-1, 'p', 'LineWidth', plt.lwCF, 'Color', 'b', 'MarkerSize', plt.mrkSizeCF);

plot([hi.tailFreq hi.tailFreq], [hi.tailThresh hi.tipThresh], '-', 'LineWidth', plt.lw1, 'Color', get_color('r'));
plot(hi.tailFreq, hi.tailThresh-1, '^', 'LineWidth', plt.lw1, 'Color', 'r', 'MarkerSize', plt.mrkSize0);
plot(hi.tailFreq, hi.tipThresh+1, 'v', 'LineWidth', plt.lw1, 'Color', 'r', 'MarkerSize', plt.mrkSize0);
plot(hi.data.CF_Hz/1e3, hi.tipThresh-1, 'p', 'LineWidth', plt.lwCF, 'Color', 'r', 'MarkerSize', plt.mrkSizeCF);
ylabel('TC Thresh. (dB SPL)');

sp_ax(3)= subplot(235);
% yyaxis left;
hold on;
[Pxx_nh_dB,freq_nh,lHan_dT_PSD_nh]= plot_dpss_psd(nh.dT_window/rms(nh.dT_window), fs, 'xunit', 'khz');
set(lHan_dT_PSD_nh, 'Color', get_color('b'), 'LineWidth', plt.lw1);
[~, max_ind_nh]= max(Pxx_nh_dB);
arrowHan(1)= text(freq_nh(max_ind_nh), Pxx_nh_dB(max_ind_nh)+3, '\downarrow', 'Color', get_color('b'), 'HorizontalAlignment', 'center');
arrowTxtHan(1)= text(freq_nh(max_ind_nh), Pxx_nh_dB(max_ind_nh)+6, 'F_2', 'Color', 'k', 'HorizontalAlignment', 'center');

[Pxx_hi_dB,freq_hi,lHan_dT_PSD_hi]= plot_dpss_psd(hi.dT_window/rms(hi.dT_window), fs, 'xunit', 'khz', 'yrange', 30);
set(lHan_dT_PSD_hi, 'Color', get_color('lr'), 'LineWidth', plt.lw1);
[~, max_ind_hi]= max(Pxx_hi_dB);
arrowHan(2)= text(freq_hi(max_ind_hi), Pxx_hi_dB(max_ind_hi)+3, '\downarrow', 'Color', get_color('r'), 'HorizontalAlignment', 'center');
arrowTxtHan(2)= text(freq_hi(max_ind_hi), Pxx_hi_dB(max_ind_hi)+6, 'F_1', 'Color', get_color('k'), 'HorizontalAlignment', 'center');
set(gca, 'YColor', 'k', 'XTick', plt.freqTick);
ylabel('');
ylim([-50 -22]);
xlim(plt.PSDxlim);
xlabel('');

txtHan(3)= text(.05, .99, 'C', 'Units', 'normalized', 'FontWeight', 'bold');

set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);
set(txtHan, 'fontsize', plt.ttl_fSize);
set(arrowHan, 'fontsize', 20);


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
lHan(1)= plot(cf_kHz(nhInds), TTR_dB(nhInds), 'x', 'color', get_color('b'), 'MarkerSize', plt.mrkSize);
lHan(2)= plot(cf_kHz(hiInds), TTR_dB(hiInds), '+', 'color', get_color('r'), 'MarkerSize', plt.mrkSize);
set(gca, 'xscale', 'log', 'xtick', CFcenters, 'TickLength', plt.tick_len);
xlabel('CF (kHz)');
ylabel('TTR (dB)');
txtHan(4)= text(.05, .99, 'D', 'Units', 'normalized', 'FontWeight', 'bold');
[legHan, icons]= legend(lHan, 'NH', 'HI', 'Location', 'northwest', 'box', 'off');
icons(1).Position(1)= icons(1).Position(1)-0.2;
icons(2).Position(1)= icons(2).Position(1)-0.2;

%% define new axes for AB
Xcorner_X= .06;
Xwidth_X= .245;
Xshift_X= .08;
Ycorner_X= .09;
Ywidth_X= .39;
Yshift_X= .1;

set(sp_ax(1),'Position',[Xcorner_X Ycorner_X+Ywidth_X+Yshift_X 3*Xwidth_X+2.4*Xshift_X Ywidth_X])
drawnow

set(sp_ax(2),'Position',[Xcorner_X Ycorner_X Xwidth_X Ywidth_X])
drawnow

set(sp_ax(3),'Position',[Xcorner_X+Xwidth_X+1.5*Xshift_X Ycorner_X Xwidth_X Ywidth_X])
drawnow

set(sp_ax(4),'Position',[Xcorner_X+2*Xwidth_X+2.4*Xshift_X Ycorner_X Xwidth_X Ywidth_X])
drawnow


if saveFig
    print([dirStruct.png 'Fig3_dt_demo_TTR'], '-dpng',  '-r600');
end
