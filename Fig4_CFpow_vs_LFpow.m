clear;
clc;

saveFig= 1;
saveStats= 1;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];

figSize_cm= [15 5 17 15];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

%% All other NH chins
dirStruct.Root_savingDir='';

anl.NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
anl.HIchins= [358 360 361 362 367 370];
all_chinDanishData= helper.load_danish_chin_data(anl);

%%
fs= 20e3;
anl.stimDuration= 1300e-3;
anl.fs= fs; % 1/(bin resolution) for uRate/histcount
anl.nw= 3;
anl.tStart= 0; % 0.23;
anl.tEnd= 1.3; %0.66; %1.3;
anl.binRes= 1/anl.fs;
anl.tEdge_hist= anl.tStart:anl.binRes:anl.tEnd;
anl.tBinCenter= (anl.tEdge_hist(1:end-1) + anl.tEdge_hist(2:end))/2;

anl.LF_co= 400;
anl.CFlowlim_kHz= .7;
anl.CFuplim_kHz= 5.1;

%% Load stim
[sig, fsSig]= audioread(['mat_data' filesep 'FLN_Stim_S_P.wav']);
sig= gen_resample(sig, fsSig, fs);
danish.voiced_boundaries= helper.find_voicing_boundaries(sig, fs, 0, .13);

temp_f0= load(['mat_data' filesep 'danish_pitch.mat']);
temp_f0= temp_f0.pitch_data;
danish.voiced_inds= any(anl.tBinCenter>danish.voiced_boundaries(:,1) & anl.tBinCenter<danish.voiced_boundaries(:,2), 1);

%% get voiced driven rate
VoicedDrivenRate= nan(length(all_chinDanishData), 1);
PSD_data_dB_cf= nan(length(all_chinDanishData), 1);
PSD_data_dB_lf= nan(length(all_chinDanishData), 1);

parfor unitVar= 1:length(all_chinDanishData)
    cur_unit_data= all_chinDanishData(unitVar);
    cur_CF_Hz= cur_unit_data.CF_Hz;
    
    if ~isempty(cur_unit_data.pos) & ~isnan(cur_CF_Hz) %#ok<AND2>
        %% SN
        temp_uRate_pos= histcounts( cell2mat(cur_unit_data.pos)-cur_unit_data.delay, anl.tEdge_hist) / length(cur_unit_data.pos); %#ok<PFBNS>
        temp_uRate_neg= histcounts( cell2mat(cur_unit_data.neg)-cur_unit_data.delay, anl.tEdge_hist) / length(cur_unit_data.neg);
        temp_uRate_sum= (temp_uRate_pos + temp_uRate_neg) / 2;
        temp_uRate_diff= (temp_uRate_pos - temp_uRate_neg)/2;
        temp_uRate_diff= temp_uRate_diff .* danish.voiced_inds; %#ok<PFBNS>
        
        [dT_PSD, freq_PSD]= plot_dpss_psd(temp_uRate_diff, fs, 'nw', anl.nw);
        
        VoicedDrivenRate(unitVar)= sum(temp_uRate_sum .* danish.voiced_inds);
        PSD_data_dB_cf(unitVar)= helper.psd_cf_energy(cur_CF_Hz, dT_PSD, freq_PSD, fs);
        PSD_data_dB_lf(unitVar)= helper.psd_f0_energy(anl.LF_co, dT_PSD, freq_PSD, fs);
        
    end
end


%% Plot
plt.lw1= 1.25;
plt.lw2= 2;
plt.lw0= 0.75;
plt.mrkSize= 4;
plt.freqTick= [.5 1 2 3 4 5];
plt.octRange4Avging= .33;
plt.plotVar= true;
plt.MINpts= 2;
plt.avgType= 'mean'; % 'weighted-mean' | 'mean' | median
plt.fSize= 9;
plt.ttl_fSize= 11;
plt.yLim_DynRange= 15;

all_cf_kHz= [all_chinDanishData.CF_Hz]/1e3;
all_Q10_local= [all_chinDanishData.Q10local];
all_Q10_global= [all_chinDanishData.Q10global];
all_thresh_dB= [all_chinDanishData.thresh_dB];
all_TTR_dB= [all_chinDanishData.TTR_dB];
all_SR_persec= [all_chinDanishData.SR];

validRateInds= VoicedDrivenRate>0 & all_cf_kHz(:)>anl.CFlowlim_kHz & all_cf_kHz(:)<anl.CFuplim_kHz;
nhInds= ismember([all_chinDanishData.chinID], anl.NHchins)' & validRateInds;
hiInds= ismember([all_chinDanishData.chinID], anl.HIchins)' & validRateInds;

all_sr= sum(danish.voiced_boundaries(:,2) - danish.voiced_boundaries(:,1)) * [all_chinDanishData.SR];

figure(1);
clf;

sp_ax(1)= subplot(4, 1, 1);
hold on;

plot(all_cf_kHz(nhInds), VoicedDrivenRate(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[meanRateNH, meanRateFreqNH, lHan] = helper.octAVG(all_cf_kHz(nhInds), VoicedDrivenRate(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw2, 'LineStyle', '-', 'Color', 'b');

plot(all_cf_kHz(hiInds), VoicedDrivenRate(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[meanRateHI, meanRateFreqHI, lHan]= helper.octAVG( all_cf_kHz(hiInds), VoicedDrivenRate(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw2, 'LineStyle', '-', 'Color', 'r');

plot(all_cf_kHz(nhInds), all_sr(nhInds), '.', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize, 'linew', plt.lw0);
[meanRateNH_SR, meanRateFreqNH_SR, lHan]= helper.octAVG( all_cf_kHz(nhInds), all_sr(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw1, 'linestyle', '--', 'Color', 'b');

plot(all_cf_kHz(hiInds), all_sr(hiInds), '.', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize, 'linew', plt.lw0);
[meanRateHI_SR, meanRateFreqHI_SR, lHan]= helper.octAVG( all_cf_kHz(hiInds), all_sr(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw1, 'linestyle', '--', 'Color', 'r');

set(gca, 'XScale', 'log', 'XTick', plt.freqTick);
ylim([0 151]);
ylabel('Rate_{voiced} (spikes/s)');


cfRange_min_kHz= 0;
cfRange_max_kHz= 5;
valid_cf= (all_cf_kHz(:)>cfRange_min_kHz) & (all_cf_kHz(:)<cfRange_max_kHz);
[~, pVoicedFrac]= ttest2(VoicedDrivenRate(nhInds&valid_cf), VoicedDrivenRate(hiInds&valid_cf));
fprintf('Driven rate p=%.3f: %.1f < CF < %.1f kHz\n', pVoicedFrac, cfRange_min_kHz, cfRange_max_kHz);

lHan(1)= plot(nan, nan, '-', 'Color', 'b', 'MarkerSize', plt.mrkSize, 'linew', plt.lw2);
lHan(2)= plot(nan, nan, '-', 'Color', 'r', 'MarkerSize', plt.mrkSize, 'linew', plt.lw2);
lHan(3)= plot(nan, nan, '--', 'Color', 'b', 'MarkerSize', plt.mrkSize, 'linew', plt.lw1);
lHan(4)= plot(nan, nan, '--', 'Color', 'r', 'MarkerSize', plt.mrkSize, 'linew', plt.lw1);

legend_str= {'NH DR', 'HI DR', 'NH SR', 'HI SR'};
[legHan, icons]= leg_helper.columnlegend(lHan, length(lHan), legend_str, 'Location', 'NorthEast', 'boxoff');
legHan.Position(2)= .89;

count= 5;
icons(count).XData= mean(icons(count).XData) + [-.04 +.04]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [-.04 +.04]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [-.04 +.04]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [-.04 +.04]; % count= count+1;


sp_ax(2)= subplot(4, 1, 2);
hold on;
plot(all_cf_kHz(nhInds), PSD_data_dB_cf(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan(1)]= helper.octAVG( all_cf_kHz(nhInds), PSD_data_dB_cf(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan(1), 'linew', plt.lw2, 'linestyle', '-', 'Color', 'b');
plot(all_cf_kHz(hiInds), PSD_data_dB_cf(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan(2)]= helper.octAVG( all_cf_kHz(hiInds), PSD_data_dB_cf(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan(2), 'linew', plt.lw2, 'linestyle', '-', 'Color', 'r');
set(gca, 'XScale', 'log', 'XTick', plt.freqTick);
ylabel('Power_{CF} (dB)');
ylim([max(ylim)-plt.yLim_DynRange*.67 max(ylim)])
[legHan, icons]= legend(lHan, 'NH', 'HI', 'Location', 'southwest', 'box', 'off');
count= 3;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; count= count+2;
icons(count).XData= mean(icons(count).XData) + [0.0 +.25]; 
clear lHan;

cfRange_min_kHz= 0;
cfRange_max_kHz= 3;
valid_cf= (all_cf_kHz(:)>cfRange_min_kHz) & (all_cf_kHz(:)<cfRange_max_kHz);
[~, pVoicedFrac]= ttest2(PSD_data_dB_cf(nhInds&valid_cf), PSD_data_dB_cf(hiInds&valid_cf));
fprintf('Pow_nearCF p=%.4f: %.1f < CF < %.1f kHz\n', pVoicedFrac, cfRange_min_kHz, cfRange_max_kHz);

cfRange_min_kHz= 3;
cfRange_max_kHz= 10;
valid_cf= (all_cf_kHz(:)>cfRange_min_kHz) & (all_cf_kHz(:)<cfRange_max_kHz);
[~, pVoicedFrac]= ttest2(PSD_data_dB_cf(nhInds&valid_cf), PSD_data_dB_cf(hiInds&valid_cf));
fprintf('Pow_nearCF p=%.4f: %.1f < CF < %.1f kHz\n', pVoicedFrac, cfRange_min_kHz, cfRange_max_kHz);


sp_ax(3)= subplot(4, 1, 3);
hold on;
plot(all_cf_kHz(nhInds), PSD_data_dB_lf(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan]= helper.octAVG( all_cf_kHz(nhInds), PSD_data_dB_lf(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw2, 'linestyle', '-', 'Color', 'b');
plot(all_cf_kHz(hiInds), PSD_data_dB_lf(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan]= helper.octAVG( all_cf_kHz(hiInds), PSD_data_dB_lf(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw2, 'linestyle', '-', 'Color', 'r');
set(gca, 'XScale', 'log', 'XTick', plt.freqTick);
ylabel('Power_{LF} (dB)');
ylim([max(ylim)-plt.yLim_DynRange*.67 max(ylim)])

cfRange_min_kHz= .6;
cfRange_max_kHz= 5;
valid_cf= (all_cf_kHz(:)>cfRange_min_kHz) & (all_cf_kHz(:)<cfRange_max_kHz);
[~, pVoicedFrac]= ttest2(PSD_data_dB_lf(nhInds&valid_cf), PSD_data_dB_lf(hiInds&valid_cf));
fprintf('PfracLF p=%.4f: %.1f < CF < %.1f kHz\n', pVoicedFrac, cfRange_min_kHz, cfRange_max_kHz);

sp_ax(4)= subplot(4, 1, 4);
hold on;
plot(all_cf_kHz(nhInds), PSD_data_dB_cf(nhInds) - PSD_data_dB_lf(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan]= helper.octAVG( all_cf_kHz(nhInds), PSD_data_dB_cf(nhInds) - PSD_data_dB_lf(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw2, 'linestyle', '-', 'Color', 'b');
plot(all_cf_kHz(hiInds), PSD_data_dB_cf(hiInds) - PSD_data_dB_lf(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan]= helper.octAVG( all_cf_kHz(hiInds), PSD_data_dB_cf(hiInds) - PSD_data_dB_lf(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'linew', plt.lw2, 'linestyle', '-', 'Color', 'r');
set(gca, 'XScale', 'log', 'XTick', plt.freqTick);
ylabel('Power_{CF / LF} (dB)');
ylim([max(ylim)-plt.yLim_DynRange max(ylim)])
linkaxes(sp_ax, 'x');
xlim(sp_ax(1), [anl.CFlowlim_kHz, anl.CFuplim_kHz]);
xlabel('Characteristic Frequency (kHz)');

set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);
txtHan= add_subplot_letter(4, 1, plt.ttl_fSize, .02, .95);

%
Xcorner_X= .075;
Xwidth_X= .89;
Ycorner_X= .064;
Ywidth_X= .192;
Yshift_X= .05;

% D
set(sp_ax(4),'Position',[Xcorner_X, Ycorner_X, Xwidth_X, Ywidth_X])
drawnow
% C
set(sp_ax(3),'Position',[Xcorner_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner_X, Ycorner_X+2*Ywidth_X+2*Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% A
set(sp_ax(1),'Position',[Xcorner_X, Ycorner_X+3*Ywidth_X+3*Yshift_X, Xwidth_X, Ywidth_X])
drawnow

if saveFig
    print([dirStruct.png 'Fig4_dt_CFvsLF'], '-dpng',  '-r600');
end


%% 
thresh_cutoff= -inf;
freq_kHz_cutoffs= [.7 5.1];

valid_lowCFinds= (all_cf_kHz(:) > min(freq_kHz_cutoffs)) & (all_cf_kHz(:) < max(freq_kHz_cutoffs));
valid_thresh_inds= all_thresh_dB(:)>thresh_cutoff;

nh_valid_all= nhInds&valid_lowCFinds&valid_thresh_inds;
hi_valid_all= hiInds&valid_lowCFinds&valid_thresh_inds;

chinStatus= nan(length(nhInds), 1);
chinStatus(ismember([all_chinDanishData.chinID], anl.NHchins)')= 1;
chinStatus(ismember([all_chinDanishData.chinID], anl.HIchins)')= 0;
all_chinIDs= [all_chinDanishData.chinID];

% with NaN version 
% notNanInds= ~isnan( all_Q10_local(:) + all_TTR_dB(:) + PSD_data_dB_cf(:) + chinStatus(:) + all_SR_persec(:) + all_thresh_dB(:) );
notNanInds= ones(size(all_Q10_local(:)));
valid_CF_rate_thresh_inds= valid_lowCFinds & valid_thresh_inds & validRateInds;
valid_tab_inds= notNanInds & valid_CF_rate_thresh_inds;

var_tab_chinID= all_chinIDs(valid_tab_inds);
var_tab_CF_kHz= all_cf_kHz(valid_tab_inds);
var_tab_Q10_local= all_Q10_local(valid_tab_inds);
var_tab_Q10_global= all_Q10_global(valid_tab_inds);
var_tab_TTR_dB= all_TTR_dB(valid_tab_inds);
var_tab_CF2LF_dB= PSD_data_dB_cf(valid_tab_inds) - PSD_data_dB_lf(valid_tab_inds);
var_tab_CFpow_dB= PSD_data_dB_cf(valid_tab_inds);
var_tab_LFpow_dB= PSD_data_dB_lf(valid_tab_inds);
var_tab_hearingStatus= chinStatus(valid_tab_inds);
var_tab_SR= all_SR_persec(valid_tab_inds);
var_tab_thresh_dB= all_thresh_dB(valid_tab_inds);
var_tab_drivenRate= VoicedDrivenRate(valid_tab_inds);

% without NaN version 
tbl_DT_CF2LF_withNan = table(var_tab_chinID(:), var_tab_CF_kHz(:), var_tab_thresh_dB(:), var_tab_Q10_local(:), var_tab_Q10_global(:), var_tab_TTR_dB(:), ...
    var_tab_SR(:), var_tab_hearingStatus(:), var_tab_CF2LF_dB(:), var_tab_CFpow_dB(:), var_tab_LFpow_dB(:), var_tab_drivenRate(:),  ...
    'VariableNames', {'ChinID', 'CF_kHz', 'Threshold_dBSPL', 'Q10local',       'Q10global',       'TTR_dB',       ...
    'SpontRate',      'HearingStatus',           'CF2LFpow',           'CFpow',                'LFpow',             'DrivenRate'});

if saveStats
    writetable(tbl_DT_CF2LF_withNan , [dirStruct.stats 'Fig4_tbl_DT_CF2LF_withNan.txt'], 'Delimiter', 'tab');
end