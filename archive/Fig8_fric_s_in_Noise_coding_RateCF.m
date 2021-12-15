% Look at modulation filtered meanrates.
clear;
clc;

saveFig= 0;
saveStats= 0;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];
dirStruct.loading_dir= ['ANF_Data' filesep];

use_LinSub1_QuadSub0= 1; % Reviewer wanted to see the effect of sqrt [ corr(SN,S).^2 - corr(SN,N).^2 ]

chinIDs= [321 322 325 338 341 343 346 347 354 355 358 360 361 362 367 370 373 379];
HIchins= [358 360 361 362 367 370];

figSize_cm= [15 5 17.6 12];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

%% Important params
fs= 5e3;
anl.N_lp=2;
anl.modFreq= 32;
anl.noise2use= 'SSN';
anl.tStart= 740e-3;
anl.tEnd= 840e-3;
anl.plot_tStart= 150e-3;
anl.plot_tEnd= 1.0;
anl.tBinEdges= (anl.plot_tStart):(1/fs):(anl.plot_tEnd);
anl.tBinCenters= ( anl.tBinEdges(1:end-1) + anl.tBinEdges(2:end) )/2;

%% stats params
stat_params.cf_min= .5;
stat_params.cf_max= 8;
stat_params.fric_range_kHz= [3 8];
stat_params.valid_range_vals= nan(3,2);
stat_params.invalid_range_vals= nan(3,2);

%%
meanrate_binwidth=1/fs; % good for mod freq upto 1/2/binWidth
all_ChinSpikeData=[];

for chinVar=1:length(chinIDs)
    curChinID=chinIDs(chinVar);
    curFileMeta=dir(sprintf('%s*%d*',dirStruct.loading_dir,curChinID));
    if isempty(curFileMeta)
        warning('Skipping chin %d. No directory found. Run data_SNRenv_analysis on this chin.', curChinID);
    else
        if length(curFileMeta)>1
            error('should not come here');
        else
            curFile=[curFileMeta.folder filesep curFileMeta.name];
        end
        
        load(curFile);
        if isfield(spike_data, 'thresh') % earlier spike_data created using mr_sEPSM do not have thresh
            spike_data=rmfield(spike_data, 'thresh');
        end
        all_ChinSpikeData=[all_ChinSpikeData, spike_data];  %#ok<*AGROW>
    end
end

all_ChinSpikeData= all_ChinSpikeData(strcmp({all_ChinSpikeData.noise}, anl.noise2use));
all_chin_track_unit_snr_mat= [[all_ChinSpikeData.chinID]', [all_ChinSpikeData.track]', [all_ChinSpikeData.unit]', [all_ChinSpikeData.SNR]'];

unique_chin_track_unit_snr_mat=unique(all_chin_track_unit_snr_mat, 'rows');
fricative_corr_Data=repmat(struct(...
    'CF_Hz', nan, 'SR', nan, 'SNR', nan, 'chinID', nan, 'track', nan', 'unit', nan, ...
    's_sn', nan, 's_n', nan, 'sn_n', nan), ...
    size(unique_chin_track_unit_snr_mat, 1), 1);


%% Main loop
parfor iterVar=1:size(unique_chin_track_unit_snr_mat,1) % 18 is the first 0 dB SNR ind for unique_chin_snr_track_unit_mat
% for iterVar= 1:size(unique_chin_track_unit_snr_mat,1) % 18 is the first 0 dB SNR ind for unique_chin_snr_track_unit_mat
    
    cur_ind=find(ismember(all_chin_track_unit_snr_mat, unique_chin_track_unit_snr_mat(iterVar, :), 'rows'));
    curMetaData=all_ChinSpikeData(cur_ind);
    curSpikeData=curMetaData.SpikeTrains;
    
    if ismember(curMetaData.chinID, HIchins)
        fricative_corr_Data(iterVar).hearingStatus= 'HI';
    else
        fricative_corr_Data(iterVar).hearingStatus= 'NH';
    end
    
    %% check new or old NEL system
    if curMetaData.chinID>368
        cur_delay= 2.52e-3; % compensate for inv-calib FIR filter delay
    else
        cur_delay= 0;
    end
    
    [S_stim_plus, fs_stim]=audioread(fullfile('stimuli', 'SNR_0', 'FLN_Stim_S_P.wav'));
    t_stim=(1:length(S_stim_plus))/fs_stim;
    dur_stim=length(S_stim_plus)/fs_stim;
    
    %% Load spikes
    S_rate_plus=sort(cell2mat(curSpikeData{1,1})) - cur_delay;
    S_spike_minus=sort(cell2mat(curSpikeData{1,2})) - cur_delay;
    
    N_spike_plus=sort(cell2mat(curSpikeData{2,1})) - cur_delay;
    N_spike_minus=sort(cell2mat(curSpikeData{2,2})) - cur_delay;
    
    SN_spike_plus=sort(cell2mat(curSpikeData{3,1})) - cur_delay;
    SN_spike_minus=sort(cell2mat(curSpikeData{3,2})) - cur_delay;
    
    %% Create rate vector using histogram
    rate_binWidth=min(meanrate_binwidth);
    hist_edges=0:rate_binWidth:dur_stim;
    t_hist=mean([hist_edges(1:end-1); hist_edges(2:end)],1);
    
    valid_inds= (t_hist>anl.tStart) & (t_hist<anl.tEnd);
    
    S_rate_plus= histcounts(S_rate_plus,hist_edges)/numel(curSpikeData{1,1});
    S_rate_minus= histcounts(S_spike_minus,hist_edges)/numel(curSpikeData{1,2});
    S_rate_env= (S_rate_plus+S_rate_minus)/2;
    
    N_rate_plus=histcounts(N_spike_plus,hist_edges)/numel(curSpikeData{2,1});
    N_rate_minus=histcounts(N_spike_minus,hist_edges)/numel(curSpikeData{2,2});
    N_rate_env=(N_rate_plus+N_rate_minus)/2;
    
    SN_rate_plus=histcounts(SN_spike_plus,hist_edges)/numel(curSpikeData{3,1});
    SN_rate_minus=histcounts(SN_spike_minus,hist_edges)/numel(curSpikeData{3,2});
    SN_rate_env=(SN_rate_plus+SN_rate_minus)/2;
    
    
    [corrStruct, rmsStruct]= helper.uRate_lpf_corr_S_N_SN(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
        SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, anl.modFreq, anl.N_lp, anl.tStart, anl.tEnd);
    
    fricative_corr_Data(iterVar).CF_Hz= curMetaData.CF_Hz;
    fricative_corr_Data(iterVar).thresh_dB= curMetaData.thresh_dB;
    fricative_corr_Data(iterVar).TTR_dB= curMetaData.TTR_dB;
    fricative_corr_Data(iterVar).Q10local= curMetaData.Q10local;
    fricative_corr_Data(iterVar).Q10global= curMetaData.Q10global;
    fricative_corr_Data(iterVar).SR= curMetaData.SR;
    fricative_corr_Data(iterVar).SNR= curMetaData.SNR;
    fricative_corr_Data(iterVar).chinID= curMetaData.chinID;
    fricative_corr_Data(iterVar).track= curMetaData.track;
    fricative_corr_Data(iterVar).unit= curMetaData.unit;
    fricative_corr_Data(iterVar).noise= curMetaData.noise;
    
    fricative_corr_Data(iterVar).DR_SN= nansum(SN_rate_env(valid_inds)) / (anl.tEnd-anl.tStart);
    fricative_corr_Data(iterVar).DR_S= nansum(S_rate_env(valid_inds)) / (anl.tEnd-anl.tStart);
    fricative_corr_Data(iterVar).DR_N= nansum(N_rate_env(valid_inds)) / (anl.tEnd-anl.tStart);

    
    fricative_corr_Data(iterVar).s_sn= (corrStruct.s_sn.pos + corrStruct.s_sn.neg)/2;
    fricative_corr_Data(iterVar).s_n= (corrStruct.s_n.pos + corrStruct.s_n.neg)/2;
    fricative_corr_Data(iterVar).sn_n= (corrStruct.sn_n.pos + corrStruct.sn_n.neg)/2;
    
    fricative_corr_Data(iterVar).rms_S= sqrt(rmsStruct.s.pos*rmsStruct.s.neg);
    fricative_corr_Data(iterVar).rms_SN= sqrt(rmsStruct.sn.pos*rmsStruct.sn.neg);
    fricative_corr_Data(iterVar).rms_N= sqrt(rmsStruct.n.pos*rmsStruct.n.neg);
end

%%
nhInds= strcmp({fricative_corr_Data.hearingStatus}, 'NH');
hiInds= strcmp({fricative_corr_Data.hearingStatus}, 'HI') ;
SR= [fricative_corr_Data.SR];
CF_Hz= [fricative_corr_Data.CF_Hz];
CF_kHz= CF_Hz/1e3;
Q10local= [fricative_corr_Data.Q10local];
Q10global= [fricative_corr_Data.Q10global];
thresh_dB= [fricative_corr_Data.thresh_dB];
TTR_dB= [fricative_corr_Data.TTR_dB];
corr_s_sn= [fricative_corr_Data.s_sn];
corr_s_n= [fricative_corr_Data.s_n];
rms_s= [fricative_corr_Data.rms_S];
rms_sn= [fricative_corr_Data.rms_SN];

Rate_SN_Hz= [fricative_corr_Data.DR_SN];
Rate_N_Hz= [fricative_corr_Data.DR_N];
Rate_SNrelN= max(1, Rate_SN_Hz - Rate_N_Hz);
% Rate_SNrelN= max(-inf, Rate_SN_Hz - Rate_N_Hz);
% warning('remove this monstrosity')

if use_LinSub1_QuadSub0
    corr_snr_2plot= (corr_s_sn-corr_s_n);
else 
    corr_snr_2plot= sqrt(corr_s_sn.^2 - corr_s_n.^2);
end

snrs2use= sort(unique([fricative_corr_Data.SNR]), 'descend');

SRboundary= 18;
hsrInds= SR>SRboundary;
lmsrInds= SR<SRboundary;


%%

demo.NH_chin_track_unit_snr= [373, 1, 6, 0]; % One more NH example = [322, 3, 2] | [355, 2, 1, 0]
demo.HI_chin_track_unit_snr= [361, 1, 2, 0]; % HI example = [361, 1, 3]
demo.nhInd= find(ismember(all_chin_track_unit_snr_mat, demo.NH_chin_track_unit_snr, 'rows'));
demo.hiInd= find(ismember(all_chin_track_unit_snr_mat, demo.HI_chin_track_unit_snr, 'rows'));

demo.nh_data= all_ChinSpikeData(demo.nhInd);
demo.hi_data= all_ChinSpikeData(demo.hiInd);
demo.delay= 0;

demo.nh_uRate_pos_S= histcounts( cell2mat(demo.nh_data.SpikeTrains{1,1}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{1,1}) * fs;
demo.nh_uRate_pos_N= histcounts( cell2mat(demo.nh_data.SpikeTrains{2,1}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{2,1}) * fs;
demo.nh_uRate_pos_SN= histcounts( cell2mat(demo.nh_data.SpikeTrains{3,1}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{3,1}) * fs;

demo.hi_uRate_pos_S= histcounts( cell2mat(demo.hi_data.SpikeTrains{1,1}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{1,1}) * fs;
demo.hi_uRate_pos_N= histcounts( cell2mat(demo.hi_data.SpikeTrains{2,1}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{2,1}) * fs;
demo.hi_uRate_pos_SN= histcounts( cell2mat(demo.hi_data.SpikeTrains{3,1}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{3,1}) * fs;

demo.nh_uRate_neg_S= histcounts( cell2mat(demo.nh_data.SpikeTrains{1,2}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{1,2}) * fs;
demo.nh_uRate_neg_N= histcounts( cell2mat(demo.nh_data.SpikeTrains{2,2}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{2,2}) * fs;
demo.nh_uRate_neg_SN= histcounts( cell2mat(demo.nh_data.SpikeTrains{3,2}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{3,2}) * fs;

demo.hi_uRate_neg_S= histcounts( cell2mat(demo.hi_data.SpikeTrains{1,2}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{1,2}) * fs;
demo.hi_uRate_neg_N= histcounts( cell2mat(demo.hi_data.SpikeTrains{2,2}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{2,2}) * fs;
demo.hi_uRate_neg_SN= histcounts( cell2mat(demo.hi_data.SpikeTrains{3,2}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{3,2}) * fs;

curFilt= designfilt('lowpassiir','FilterOrder',anl.N_lp, 'PassbandFrequency',anl.modFreq,'PassbandRipple',0.2, 'SampleRate',fs);
demo.nhAshift= 1.3*max([demo.nh_uRate_pos_S, demo.nh_uRate_pos_N, demo.nh_uRate_pos_SN]);
demo.hiAshift= max([demo.hi_uRate_pos_S, demo.hi_uRate_pos_N, demo.hi_uRate_pos_SN]);

demo.nh_LPenv_S= 4*filtfilt(curFilt, ( demo.nh_uRate_pos_S + demo.nh_uRate_neg_S) / 2);
demo.nh_LPenv_SN= 4*filtfilt(curFilt, ( demo.nh_uRate_pos_SN + demo.nh_uRate_neg_SN) / 2);
demo.nh_LPenv_N= 4*filtfilt(curFilt, ( demo.nh_uRate_pos_N + demo.nh_uRate_neg_N) / 2);

demo.hi_LPenv_S= 4*filtfilt(curFilt, ( demo.hi_uRate_pos_S + demo.hi_uRate_neg_S) / 2);
demo.hi_LPenv_SN= 4*filtfilt(curFilt, ( demo.hi_uRate_pos_SN + demo.hi_uRate_neg_SN) / 2);
demo.hi_LPenv_N= 4*filtfilt(curFilt, ( demo.hi_uRate_pos_N + demo.hi_uRate_neg_N) / 2);


%% Get stims for SNR plot
snr2use = 0;
noise2use = 'SSN';

[stimS, fsSig]=audioread(fullfile('stimuli', 'SNR_0', 'FLN_Stim_S_P.wav'));
[stimN, ~]=audioread(fullfile('stimuli', 'SNR_0', sprintf('%s_Stim0dB_N_P.wav', noise2use)));

stimS= helper.gen_rescale(stimS, 65);
stimN= helper.gen_rescale(stimN, 65 - snr2use);
tStim= (1:length(stimS))/fsSig;

tValidInds= tStim>anl.tStart & tStim<anl.tEnd;
tPlotInds= tStim>anl.plot_tStart & tStim<anl.plot_tEnd;

%% Plot
plt.plotVar= true;
plt.octRange4Avging= 0.75;
plt.MINpts= 1;
plt.lw0= .75;
plt.lw1= 1;
plt.lw2= 1.5;
plt.leg_fSize= 8;
plt.fSize= 9;
plt.ttl_fSize= 11;
plt.xtick_vals_cf_kHz= [.5 1 2 4 8];
plt.xtick_labs_cf_kHz= cellfun(@(x) num2str(x), num2cell(plt.xtick_vals_cf_kHz), 'UniformOutput', false);
plt.nw= 5;
plt.colGray= [150 150 150]/255;
plt.txtHan_X= .05;
plt.txtHan_Y= 1.00;
plt.mrkSize= 4;
plt.mrkSizeLeg= 6;
plt.ytick_psd= -80:10:-60;
plt.txt_shift_time= .05;
plt.txt_shift_Y_nh= .15*fs;
plt.txt_shift_Y_hi= .4*fs;
plt.tick_len= [.03 .02];
plt.ytick_tc= [0 40 80];
plt.time_xlim= [.15 1];

figure(1);
clf;
sp_ax= nan(length(snrs2use)*3, 1);
sp_bx= nan(length(snrs2use)*4, 1);

%  Spectra
sp_ax(1)= subplot(length(snrs2use), 3, 1);
yyaxis left;
hold on;
plot(demo.nh_data.TC.freqkHz, demo.nh_data.TC.TCfit, '-', 'Color', 'b', 'linew', plt.lw2);
plot(demo.hi_data.TC.freqkHz, demo.hi_data.TC.TCfit, '-', 'Color', 'r', 'linew', plt.lw2);
set(gca, 'YColor', helper.get_color('k'), 'XTick', plt.xtick_vals_cf_kHz, 'yTick', plt.ytick_tc);
ylim([20 140]);
xlim([.495 10.5]);
ylab_tc_han= ylabel('TC thresh. (dB SPL)', 'Units', 'normalized');


yyaxis right;
[~,~, lHan(1)]= helper.plot_dpss_psd(stimS(tValidInds), fsSig, 'nw', plt.nw, 'xunit', 'khz');
hold on
set(lHan(1), 'color', helper.get_color('g'), 'linew', plt.lw1);
[~,~, lHan(2)]= helper.plot_dpss_psd(stimN(tValidInds), fsSig, 'nw', plt.nw, 'xunit', 'khz');
ylabHan_psd= ylabel('PSD (dB/Hz)', 'Units', 'normalized');
ylabHan_psd.Position(1:2)= [1.16 .55];
xlabHan_psd= xlabel('Frequency (kHz)', 'Units', 'normalized');
xlabHan_psd.Position(2)= -.14;

set(gca, 'XTick', plt.xtick_vals_cf_kHz, 'YColor', 'k', 'YTick', plt.ytick_psd);
set(lHan(2), 'color', helper.get_color('prp'), 'linew', plt.lw1);
[legHanA, iconsA]= legend('TC-NH', 'TC-HI', 'S (/s/)', 'N', 'box', 'off', 'Location', 'southwest');
iconsA(5).XData= mean(iconsA(5).XData) + [0 +.25];
iconsA(7).XData= mean(iconsA(7).XData) + [0 +.25];
iconsA(9).XData= mean(iconsA(9).XData) + [0 +.25];
iconsA(11).XData= mean(iconsA(11).XData) + [0 +.25];
legHanA.Position(1:2)= [.16 .87];

pink_freq= [.5 1 2 4 8 12];
pink_amp= polyval([-3 -67], log2(pink_freq ./ min(pink_freq))); 
plot(pink_freq, pink_amp, '--', 'color', helper.get_color('pink'), 'linew', 2);

xlim([.5 11]);
ylim([-88 -58]);
ttlHan(1)= text(plt.txtHan_X, plt.txtHan_Y, 'A', 'Units', 'normalized');


% Demo NH
sp_ax(4)= subplot(length(snrs2use), 3, 4);
hold on;
plot(anl.tBinCenters, demo.nhAshift + demo.nh_uRate_pos_S, 'color', plt.colGray);
plot(anl.tBinCenters, demo.nhAshift + demo.nh_LPenv_S, 'color', helper.get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, demo.nh_uRate_pos_SN, 'color', plt.colGray);
plot(anl.tBinCenters, demo.nh_LPenv_SN, 'color', helper.get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, -demo.nhAshift + demo.nh_uRate_pos_N, 'color', plt.colGray);
plot(anl.tBinCenters, -demo.nhAshift + demo.nh_LPenv_N, 'color', helper.get_color('k'), 'linew', plt.lw2);
plot(tStim(tPlotInds), -1.5*demo.nhAshift + stimS(tPlotInds)/max(stimS(tPlotInds))*demo.nhAshift/5, 'color', helper.get_color('g'));
ttlHan(4)= text(plt.txtHan_X, plt.txtHan_Y, 'B. NH (0 dB SNR)', 'Units', 'normalized');
xlim(plt.time_xlim)
ylim([-.63 .65] * fs);
plot([anl.tStart anl.tStart], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
plot([anl.tEnd anl.tEnd], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
text(-.15, 1.05, 'x10^3', 'Units', 'normalized');

text(anl.tEnd+plt.txt_shift_time, demo.nhAshift+plt.txt_shift_Y_nh, 'S', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, 0+plt.txt_shift_Y_nh, 'SN', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, -demo.nhAshift+plt.txt_shift_Y_nh, 'N', 'FontWeight', 'bold');

yTickScale= 0.2 * fs;
ytickValsDefault= [0 yTickScale];
ytickVals= [ytickValsDefault demo.nhAshift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault/1e3, 1, 2)), 'UniformOutput', false);
set(sp_ax(4), 'YTick', ytickVals, 'YTickLabel', ytickLabs);

% Demo HI
sp_ax(7)= subplot(length(snrs2use), 3, 7);
hold on;

plot(anl.tBinCenters, demo.hiAshift + demo.hi_uRate_pos_S, 'color', plt.colGray);
plot(anl.tBinCenters, demo.hiAshift + demo.hi_LPenv_S, 'color', helper.get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, demo.hi_uRate_pos_SN, 'color', plt.colGray);
plot(anl.tBinCenters, demo.hi_LPenv_SN, 'color', helper.get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, -demo.hiAshift + demo.hi_uRate_pos_N, 'color', plt.colGray);
plot(anl.tBinCenters, -demo.hiAshift + demo.hi_LPenv_N, 'color', helper.get_color('k'), 'linew', plt.lw2);

plot(tStim(tPlotInds), -1.6*demo.hiAshift + stimS(tPlotInds)/max(stimS(tPlotInds))*demo.hiAshift/5, 'color', helper.get_color('g'), 'linew', plt.lw0);
xlim(plt.time_xlim)
ylim([-1.67 2] * fs);
plot([anl.tStart anl.tStart], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
plot([anl.tEnd anl.tEnd], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
ttlHan(7)= text(plt.txtHan_X, plt.txtHan_Y, 'C. HI (0 dB SNR)', 'Units', 'normalized');

text(anl.tEnd+plt.txt_shift_time, demo.hiAshift+plt.txt_shift_Y_hi, 'S', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, 0+plt.txt_shift_Y_hi, 'SN', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, -demo.hiAshift+plt.txt_shift_Y_hi, 'N', 'FontWeight', 'bold');
text(-.15, 1.05, 'x10^3', 'Units', 'normalized');

yTickScale= 0.4 * fs;
ytickValsDefault= [0 yTickScale];
ytickVals= [ytickValsDefault demo.hiAshift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault/1e3, 1, 2)), 'UniformOutput', false);
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs);

minCorrVal= 1e-3;
corrYtick_vals= [.001 .01 .1];

corr_snr_raw= corr_snr_2plot;
corr_snr_2plot(corr_snr_2plot<minCorrVal)= minCorrVal;

SPletters= 'DEFGHI';
for snrVar= 1:length(snrs2use)
    curSNR= snrs2use(snrVar);
    curSNRinds= [fricative_corr_Data.SNR]==curSNR;
    validCFinds= (CF_kHz>3) & (CF_kHz<8);
    
    
    %%
    figure(1)
    sp_ax(3*snrVar-1)= subplot(length(snrs2use), 3, 3*snrVar-1);
    hold on;
    plot(CF_kHz(hsrInds & nhInds & curSNRinds), corr_snr_2plot(hsrInds & nhInds & curSNRinds), '*', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & nhInds & curSNRinds), corr_snr_2plot(hsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');
    
    plot(CF_kHz(lmsrInds & nhInds & curSNRinds), corr_snr_2plot(lmsrInds & nhInds & curSNRinds), 's', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & nhInds & curSNRinds), corr_snr_2plot(lmsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');
    
    plot(CF_kHz(hsrInds & hiInds & curSNRinds), corr_snr_2plot(hsrInds & hiInds & curSNRinds), '*', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & hiInds & curSNRinds), corr_snr_2plot(hsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');
    
    plot(CF_kHz(lmsrInds & hiInds & curSNRinds), corr_snr_2plot(lmsrInds & hiInds & curSNRinds), 's', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & hiInds & curSNRinds), corr_snr_2plot(lmsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');

    
    set(gca, 'XScale', 'log', 'XTick', plt.xtick_vals_cf_kHz, 'YScale', 'log', 'YTick', corrYtick_vals);
    ttlHan(3*snrVar-1)= text(plt.txtHan_X, plt.txtHan_Y, SPletters(snrVar), 'Units', 'normalized');
    
    ttlHan(9+snrVar)= text(-.37, .25, sprintf('%d dB SNR', curSNR), 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Rotation', 90); 
    
    cfs_in_fric= (CF_kHz>stat_params.fric_range_kHz(1)) & (CF_kHz<stat_params.fric_range_kHz(2));
    cfs_out_fric= ( (CF_kHz>stat_params.cf_min) & (CF_kHz<stat_params.fric_range_kHz(1)) ) |  ( (CF_kHz<stat_params.cf_max) & (CF_kHz>stat_params.fric_range_kHz(2)) );
    
    [~, stat_params.valid_range_vals(snrVar, 1)]= ttest2(corr_snr_2plot(hsrInds & nhInds & curSNRinds & cfs_in_fric), corr_snr_2plot(hsrInds & hiInds & curSNRinds & cfs_in_fric));
    [~, stat_params.valid_range_vals(snrVar, 2)]= ttest2(corr_snr_2plot(lmsrInds & nhInds & curSNRinds & cfs_in_fric), corr_snr_2plot(lmsrInds & hiInds & curSNRinds & cfs_in_fric));
    
    [~, stat_params.invalid_range_vals(snrVar, 1)]= ttest2(corr_snr_2plot(hsrInds & nhInds & curSNRinds & cfs_out_fric), corr_snr_2plot(hsrInds & hiInds & curSNRinds & cfs_out_fric));
    [~, stat_params.invalid_range_vals(snrVar, 2)]= ttest2(corr_snr_2plot(lmsrInds & nhInds & curSNRinds & cfs_out_fric), corr_snr_2plot(lmsrInds & hiInds & curSNRinds & cfs_out_fric));
    
    
    
    %% Driven rate
    
    sp_ax(3*snrVar)= subplot(length(snrs2use), 3, 3*snrVar);
    hold on;
    plot(CF_kHz(hsrInds & nhInds & curSNRinds), Rate_SNrelN(hsrInds & nhInds & curSNRinds), '*', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & nhInds & curSNRinds), Rate_SNrelN(hsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');
    
    plot(CF_kHz(lmsrInds & nhInds & curSNRinds), Rate_SNrelN(lmsrInds & nhInds & curSNRinds), 's', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & nhInds & curSNRinds), Rate_SNrelN(lmsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');
    
    plot(CF_kHz(hsrInds & hiInds & curSNRinds), Rate_SNrelN(hsrInds & hiInds & curSNRinds), '*', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & hiInds & curSNRinds), Rate_SNrelN(hsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');
    
    plot(CF_kHz(lmsrInds & hiInds & curSNRinds), Rate_SNrelN(lmsrInds & hiInds & curSNRinds), 's', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & hiInds & curSNRinds), Rate_SNrelN(lmsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');
    
    ttlHan(3*snrVar)= text(plt.txtHan_X, plt.txtHan_Y, SPletters(3+snrVar), 'Units', 'normalized');
    set(gca, 'XScale', 'log', 'XTick', plt.xtick_vals_cf_kHz, 'YScale', 'log', 'YTick', [1 10 40]);
    ylim([0 50])
    xlim([.495 10.5]);
    
end

%%
figure(1);

set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);
set(ttlHan,'FontSize', plt.ttl_fSize);
set(legHanA,'FontSize', plt.leg_fSize);

axes(sp_ax(8))
xlab_han= xlabel('Characteristic frequency (kHz)', 'Units', 'normalized');
xlab_han.Position(1)= 1.15;

lHan(1)= plot(nan, nan, '*--', 'Color', 'b', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
lHan(2)= plot(nan, nan, 's-', 'Color', 'b', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
lHan(3)= plot(nan, nan, '*--', 'Color', 'r', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
lHan(4)= plot(nan, nan, 's-', 'Color', 'r', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
[legHan, icons]= legend(lHan, 'NH-HSR', 'NH-LMSR', 'HI-HSR', 'HI-LMSR', 'box', 'off', 'Location', 'northwest');
legHan.FontSize= plt.leg_fSize;
legHan.Position(1:2)= [.445 .20];
icons(5).XData= mean(icons(5).XData) + [-0.22 +.22];
icons(7).XData= mean(icons(7).XData) + [-0.22 +.22];
icons(9).XData= mean(icons(9).XData) + [-0.22 +.22];
icons(11).XData= mean(icons(11).XData) + [-0.22 +.22];

axes(sp_ax(2));
ttlHan(13)= title('Fricative-coding fidelity', 'Units', 'normalized');
ttlHan(13).Position(1:2)= [1.15 1.05];

axes(sp_ax(5));
ylab_corr= ylabel('corr(S,SN) - corr(N,SN)', 'Units', 'normalized');
ylab_corr.Position(1)= -.2;

axes(sp_ax(6));
ylabel('rate(SN)- rate(N), in spikes/s');


axes(sp_ax(7));
ylab_rate_han= ylabel('Rate (spikes/s)', 'Units', 'normalized');
ylab_rate_han.Position(2)= 1.1;
xlabel('Time (sec)');
linkaxes(sp_ax([3 5]), 'x')
xlim(sp_ax(4), [anl.plot_tStart anl.plot_tEnd]);


linkaxes(sp_ax([2 5 8]));
xlim(sp_ax(2), [.5 1.01*max(CF_kHz)]);

linkaxes(sp_ax([3 6 9]));
xlim(sp_ax(3), [.5 1.01*max(CF_kHz)]);

ylab_tc_han.Position(1)= -.19;
ylab_rate_han.Position(1)= -.19;

set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len);

%% define new axes for AB
Xcorner= .067;
Xwidth= .228;
Xshift_AD= .16;
Xshift_DG= .08;
Ycorner= .088;
Ywidth= .235;
Yshift= .08;

% G
set(sp_ax(7),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow

% H
set(sp_ax(8),'Position',[Xcorner+Xwidth+Xshift_AD Ycorner Xwidth Ywidth])
drawnow

% I
set(sp_ax(9),'Position',[Xcorner+2*Xwidth+Xshift_AD+Xshift_DG Ycorner Xwidth Ywidth])
drawnow

% D
set(sp_ax(4),'Position',[Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% E
set(sp_ax(5),'Position',[Xcorner+Xwidth+Xshift_AD Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% F
set(sp_ax(6),'Position',[Xcorner+2*Xwidth+Xshift_AD+Xshift_DG Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% A
set(sp_ax(1),'Position',[Xcorner Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift_AD Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

% C
set(sp_ax(3),'Position',[Xcorner+2*Xwidth+Xshift_AD+Xshift_DG Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

if saveFig
    if use_LinSub1_QuadSub0
        print([dirStruct.png 'Fig8_fric_SpIN_DR'], '-dpng',  '-r600');
        saveas(gcf, [dirStruct.eps 'Fig8_fric_SpIN_DR'], 'epsc');
    else
        print([dirStruct.png 'Fig8_fric_SpIN_QuadSub_DR'], '-dpng',  '-r600');
    end
end

%% Create table
all_snr= [fricative_corr_Data.SNR];
all_chinIDs= [fricative_corr_Data.chinID];
xTrack = num2cell([fricative_corr_Data.track]');
xUnit = num2cell([fricative_corr_Data.unit]');
xChinID = num2cell([fricative_corr_Data.chinID]');
UnitID= cellfun(@(x,y,z) sprintf('%d_%d_%d', x, y, z), xChinID, xTrack, xUnit, 'UniformOutput', false);

nonNanInds= ones(size(validCFinds(:)));
valid_fric_inds= nonNanInds(:) & validCFinds(:);

var_tab_chinID= all_chinIDs(valid_fric_inds);
var_tab_unitID= UnitID(valid_fric_inds);
var_tab_CF_kHz= CF_kHz(valid_fric_inds);
var_tab_thresh_dB= thresh_dB(valid_fric_inds);
var_tab_Q10_local= Q10local(valid_fric_inds);
var_tab_Q10_global= Q10global(valid_fric_inds);
var_tab_TTR_dB= TTR_dB(valid_fric_inds);
var_tab_SR= SR(valid_fric_inds);
var_tab_hearingStatus= nhInds(valid_fric_inds);
var_tab_resp_corr= corr_snr_2plot(valid_fric_inds);
var_tab_resp_raw= corr_snr_raw(valid_fric_inds);
var_tab_stim_snr= all_snr(valid_fric_inds);
var_tab_DrivenRate_corr= Rate_SN_Hz(valid_fric_inds);
var_tab_DR_SN= Rate_SN_Hz(valid_fric_inds);
var_tab_DR_N= Rate_N_Hz(valid_fric_inds);
var_tab_DR_SNrelN= Rate_SNrelN(valid_fric_inds);

table_fric_in_noise_withNan= table(var_tab_unitID(:), var_tab_chinID(:), var_tab_CF_kHz(:), var_tab_thresh_dB(:), var_tab_Q10_local(:), var_tab_Q10_global(:), var_tab_TTR_dB(:), ...
    var_tab_SR(:), var_tab_hearingStatus(:), var_tab_stim_snr(:), var_tab_resp_corr(:), var_tab_resp_raw(:),  var_tab_DR_SN(:),  var_tab_DR_N(:), var_tab_DR_SNrelN(:), ...
    'VariableNames', {                  'UnitID',       'ChinID',         'CF_kHz', 'Threshold_dBSPL',              'Q10local',           'Q10global',         'TTR_dB',       ...
    'SpontRate',      'HearingStatus',           'StimSNR',           'RespCorr',          'RespCorrRaw',          'DR_SN',          'DR_N',             'DR_SNrelN'});

if saveStats
    writetable(table_fric_in_noise_withNan, [dirStruct.stats filesep 'Fig8_table_fric_in_noise_DR.txt'], 'Delimiter', 'tab');
end