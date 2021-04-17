% Look at modulation filtered meanrates.
clear;
clc;

saveFig= 1;
saveStats= 1;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];
dirStruct.loading_dir= ['ANF_Data' filesep];

chinIDs= [321 322 325 338 341 343 346 347 354 355 358 360 361 362 367 370 373 379];
HIchins= [358 360 361 362 367 370];

figSize_cm= [15 5 17 14];
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
    
    [S_stim_plus, fs_stim]=audioread(curMetaData.StimsFNames{1,1}{1});
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
corr_snr_corrected= (corr_s_sn-corr_s_n);

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

demo.nh_uRate_pos_S= histcounts( cell2mat(demo.nh_data.SpikeTrains{1,1}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{1,1});
demo.nh_uRate_pos_N= histcounts( cell2mat(demo.nh_data.SpikeTrains{2,1}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{2,1});
demo.nh_uRate_pos_SN= histcounts( cell2mat(demo.nh_data.SpikeTrains{3,1}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{3,1});

demo.hi_uRate_pos_S= histcounts( cell2mat(demo.hi_data.SpikeTrains{1,1}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{1,1});
demo.hi_uRate_pos_N= histcounts( cell2mat(demo.hi_data.SpikeTrains{2,1}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{2,1});
demo.hi_uRate_pos_SN= histcounts( cell2mat(demo.hi_data.SpikeTrains{3,1}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{3,1});

demo.nh_uRate_neg_S= histcounts( cell2mat(demo.nh_data.SpikeTrains{1,2}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{1,2});
demo.nh_uRate_neg_N= histcounts( cell2mat(demo.nh_data.SpikeTrains{2,2}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{2,2});
demo.nh_uRate_neg_SN= histcounts( cell2mat(demo.nh_data.SpikeTrains{3,2}) - demo.delay, anl.tBinEdges) / length(demo.nh_data.SpikeTrains{3,2});

demo.hi_uRate_neg_S= histcounts( cell2mat(demo.hi_data.SpikeTrains{1,2}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{1,2});
demo.hi_uRate_neg_N= histcounts( cell2mat(demo.hi_data.SpikeTrains{2,2}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{2,2});
demo.hi_uRate_neg_SN= histcounts( cell2mat(demo.hi_data.SpikeTrains{3,2}) - demo.delay, anl.tBinEdges) / length(demo.hi_data.SpikeTrains{3,2});

curFilt= designfilt('lowpassiir','FilterOrder',anl.N_lp, 'PassbandFrequency',anl.modFreq,'PassbandRipple',0.2, 'SampleRate',fs);
demo.nhAshift= max([demo.nh_uRate_pos_S, demo.nh_uRate_pos_N, demo.nh_uRate_pos_SN]);
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

[stimS, fsSig]=audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_09_09-Q321_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
[stimN, ~]=audioread(sprintf('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_09_09-Q321_AN_NH/Signals/MH/SNRenv/SNR_0/%s_Stim0dB_N_P.wav', noise2use));

stimS= gen_rescale(stimS, 65);
stimN= gen_rescale(stimN, 65 - snr2use);
tStim= (1:length(stimS))/fsSig;

tValidInds= tStim>anl.tStart & tStim<anl.tEnd;
tPlotInds= tStim>anl.plot_tStart & tStim<anl.plot_tEnd;

%% Plot
plt.plotVar= true;
plt.octRange4Avging= 0.5;
plt.MINpts= 1;
plt.lw0= .75;
plt.lw1= 1;
plt.lw2= 1.5;
plt.leg_fSize= 7;
plt.fSize= 10;
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
plt.txt_shift_Y_nh= .15;
plt.txt_shift_Y_hi= .4;
plt.tick_len= [.02 .02];
plt.ytick_tc= [0 50 100 150];

figure(1);
clf;
sp_ax= nan(length(snrs2use), 1);
sp_bx= nan(length(snrs2use)*4, 1);

%  Spectra
sp_ax(1)= subplot(length(snrs2use), 2, 1);
yyaxis left;
hold on;
plot(demo.nh_data.TC.freqkHz, demo.nh_data.TC.TCfit, '-', 'Color', 'b', 'linew', plt.lw2);
plot(demo.hi_data.TC.freqkHz, demo.hi_data.TC.TCfit, '-', 'Color', 'r', 'linew', plt.lw2);
set(gca, 'YColor', get_color('k'), 'XTick', plt.xtick_vals_cf_kHz, 'yTick', plt.ytick_tc);
ylim([20 140]);
xlim([.495 10.5]);
ylab_tc_han= ylabel('TC Thresh. (dB SPL)', 'Units', 'normalized');


yyaxis right;
[~,~, lHan(1)]= plot_dpss_psd(stimS(tValidInds), fsSig, 'nw', plt.nw, 'xunit', 'khz');
hold on
set(lHan(1), 'color', get_color('g'), 'linew', plt.lw1);
[~,~, lHan(2)]= plot_dpss_psd(stimN(tValidInds), fsSig, 'nw', plt.nw, 'xunit', 'khz');
set(gca, 'XTick', plt.xtick_vals_cf_kHz, 'YColor', 'k', 'YTick', plt.ytick_psd);
set(lHan(2), 'color', get_color('prp'), 'linew', plt.lw1);
[legHanA, iconsA]= legend('TC-NH', 'TC-HI', 'S (/s/)', 'N', 'box', 'off', 'Location', 'southwest');
iconsA(5).XData= mean(iconsA(5).XData) + [0 +.25];
iconsA(7).XData= mean(iconsA(7).XData) + [0 +.25];
iconsA(9).XData= mean(iconsA(9).XData) + [0 +.25];
iconsA(11).XData= mean(iconsA(11).XData) + [0 +.25];
legHanA.Position(1:2)= [.04 .75];

pink_freq= [.5 1 2 4 8 12];
pink_amp= polyval([-3 -67], log2(pink_freq ./ min(pink_freq))); 
plot(pink_freq, pink_amp, '--', 'color', get_color('pink'), 'linew', 2);

xlim([.5 11]);
ylim([-88 -58]);
ttlHan(1)= text(plt.txtHan_X, plt.txtHan_Y, 'A', 'Units', 'normalized');


% Demo NH
sp_ax(3)= subplot(length(snrs2use), 2, 3);
hold on;
plot(anl.tBinCenters, demo.nhAshift + demo.nh_uRate_pos_S, 'color', plt.colGray);
plot(anl.tBinCenters, demo.nhAshift + demo.nh_LPenv_S, 'color', get_color('k'), 'linew', plt.lw2);
text(anl.tEnd+plt.txt_shift_time, demo.nhAshift+plt.txt_shift_Y_nh, 'S', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, 0+plt.txt_shift_Y_nh, 'SN', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, -demo.nhAshift+plt.txt_shift_Y_nh, 'N', 'FontWeight', 'bold');

plot(anl.tBinCenters, demo.nh_uRate_pos_SN, 'color', plt.colGray);
plot(anl.tBinCenters, demo.nh_LPenv_SN, 'color', get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, -demo.nhAshift + demo.nh_uRate_pos_N, 'color', plt.colGray);
plot(anl.tBinCenters, -demo.nhAshift + demo.nh_LPenv_N, 'color', get_color('k'), 'linew', plt.lw2);
plot(tStim(tPlotInds), -1.5*demo.nhAshift + stimS(tPlotInds)/max(stimS(tPlotInds))*demo.nhAshift/5, 'color', get_color('g'));
ttlHan(3)= text(plt.txtHan_X, plt.txtHan_Y, 'B. NH (0 dB SNR)', 'Units', 'normalized');
ylim([-.52 .5]);
plot([anl.tStart anl.tStart], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
plot([anl.tEnd anl.tEnd], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);

yTickScale= 0.15;
ytickValsDefault= [0 yTickScale];
ytickVals= [ytickValsDefault demo.nhAshift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault, 1, 2)), 'UniformOutput', false);
set(sp_ax(3), 'YTick', ytickVals, 'YTickLabel', ytickLabs);

% Demo HI
sp_ax(5)= subplot(length(snrs2use), 2, 5);
hold on;

plot(anl.tBinCenters, demo.hiAshift + demo.hi_uRate_pos_S, 'color', plt.colGray);
plot(anl.tBinCenters, demo.hiAshift + demo.hi_LPenv_S, 'color', get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, demo.hi_uRate_pos_SN, 'color', plt.colGray);
plot(anl.tBinCenters, demo.hi_LPenv_SN, 'color', get_color('k'), 'linew', plt.lw2);

plot(anl.tBinCenters, -demo.hiAshift + demo.hi_uRate_pos_N, 'color', plt.colGray);
plot(anl.tBinCenters, -demo.hiAshift + demo.hi_LPenv_N, 'color', get_color('k'), 'linew', plt.lw2);

text(anl.tEnd+plt.txt_shift_time, demo.hiAshift+plt.txt_shift_Y_hi, 'S', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, 0+plt.txt_shift_Y_hi, 'SN', 'FontWeight', 'bold');
text(anl.tEnd+plt.txt_shift_time, -demo.hiAshift+plt.txt_shift_Y_hi, 'N', 'FontWeight', 'bold');

plot(tStim(tPlotInds), -1.6*demo.hiAshift + stimS(tPlotInds)/max(stimS(tPlotInds))*demo.hiAshift/5, 'color', get_color('g'), 'linew', plt.lw0);
ylim([-1.6 2]);
plot([anl.tStart anl.tStart], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
plot([anl.tEnd anl.tEnd], [min(ylim) max(ylim)], 'm-.', 'linew', plt.lw0);
ttlHan(5)= text(plt.txtHan_X, plt.txtHan_Y, 'C. HI (0 dB SNR)', 'Units', 'normalized');

yTickScale= 0.3;
ytickValsDefault= [0 yTickScale];
ytickVals= [ytickValsDefault demo.hiAshift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault, 1, 2)), 'UniformOutput', false);
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs);

minCorrVal= 1e-3;
corrYtick_vals= [.001 .01 .1];

corr_snr_raw= corr_snr_corrected;
corr_snr_corrected(corr_snr_corrected<minCorrVal)= minCorrVal;

SPletters= 'DEF';
for snrVar= 1:length(snrs2use)
    curSNR= snrs2use(snrVar);
    curSNRinds= [fricative_corr_Data.SNR]==curSNR;
    validCFinds= (CF_kHz>3) & (CF_kHz<8);
    
    
    %%
    figure(1)
    sp_ax(2*snrVar)= subplot(length(snrs2use), 2, 2*snrVar);
    hold on;
    plot(CF_kHz(hsrInds & nhInds & curSNRinds), corr_snr_corrected(hsrInds & nhInds & curSNRinds), '*', 'Color', get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & nhInds & curSNRinds), corr_snr_corrected(hsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');
    
    plot(CF_kHz(lmsrInds & nhInds & curSNRinds), corr_snr_corrected(lmsrInds & nhInds & curSNRinds), 's', 'Color', get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & nhInds & curSNRinds), corr_snr_corrected(lmsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');
    
    plot(CF_kHz(hsrInds & hiInds & curSNRinds), corr_snr_corrected(hsrInds & hiInds & curSNRinds), '*', 'Color', get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & hiInds & curSNRinds), corr_snr_corrected(hsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');
    
    plot(CF_kHz(lmsrInds & hiInds & curSNRinds), corr_snr_corrected(lmsrInds & hiInds & curSNRinds), 's', 'Color', get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & hiInds & curSNRinds), corr_snr_corrected(lmsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');
    
    set(gca, 'XScale', 'log', 'XTick', plt.xtick_vals_cf_kHz, 'YScale', 'log', 'YTick', corrYtick_vals);
    ttlHan(2*snrVar)= text(plt.txtHan_X, plt.txtHan_Y, sprintf('%s. SNR = %d dB', SPletters(snrVar), curSNR), 'Units', 'normalized');
    
    cfs_in_fric= (CF_kHz>stat_params.fric_range_kHz(1)) & (CF_kHz<stat_params.fric_range_kHz(2));
    cfs_out_fric= ( (CF_kHz>stat_params.cf_min) & (CF_kHz<stat_params.fric_range_kHz(1)) ) |  ( (CF_kHz<stat_params.cf_max) & (CF_kHz>stat_params.fric_range_kHz(2)) );
    
    [~, stat_params.valid_range_vals(snrVar, 1)]= ttest2(corr_snr_corrected(hsrInds & nhInds & curSNRinds & cfs_in_fric), corr_snr_corrected(hsrInds & hiInds & curSNRinds & cfs_in_fric));
    [~, stat_params.valid_range_vals(snrVar, 2)]= ttest2(corr_snr_corrected(lmsrInds & nhInds & curSNRinds & cfs_in_fric), corr_snr_corrected(lmsrInds & hiInds & curSNRinds & cfs_in_fric));
    
    [~, stat_params.invalid_range_vals(snrVar, 1)]= ttest2(corr_snr_corrected(hsrInds & nhInds & curSNRinds & cfs_out_fric), corr_snr_corrected(hsrInds & hiInds & curSNRinds & cfs_out_fric));
    [~, stat_params.invalid_range_vals(snrVar, 2)]= ttest2(corr_snr_corrected(lmsrInds & nhInds & curSNRinds & cfs_out_fric), corr_snr_corrected(lmsrInds & hiInds & curSNRinds & cfs_out_fric));
    
end

%%
figure(1);
xlabel('Characteristic Frequency (kHz)');

set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);
set(ttlHan,'FontSize', 11);

lHan(1)= plot(nan, nan, '*--', 'Color', 'b', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
lHan(2)= plot(nan, nan, 's-', 'Color', 'b', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
lHan(3)= plot(nan, nan, '*--', 'Color', 'r', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
lHan(4)= plot(nan, nan, 's-', 'Color', 'r', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
[legHan, icons]= legend(lHan, 'NH-HSR', 'NH-LMSR', 'HI-HSR', 'HI-LMSR', 'box', 'off', 'Location', 'northwest');
legHan.FontSize= plt.leg_fSize;
legHan.Position(1:2)= [.61 .20];
icons(5).XData= mean(icons(5).XData) + [-0.22 +.22];
icons(7).XData= mean(icons(7).XData) + [-0.22 +.22];
icons(9).XData= mean(icons(9).XData) + [-0.22 +.22];
icons(11).XData= mean(icons(11).XData) + [-0.22 +.22];


axes(sp_ax(5));
ylab_rate_han= ylabel('Rate (spikes/bin)', 'Units', 'normalized');
ylab_rate_han.Position(2)= 1.1;

axes(sp_ax(4));
ylabel('Fricative-coding fidelity');


axes(sp_ax(5));
xlabel('Time (sec)');
linkaxes(sp_ax([3 5]), 'x')
xlim(sp_ax(3), [anl.plot_tStart anl.plot_tEnd]);

linkaxes(sp_ax([3 5]), 'x');
linkaxes(sp_ax([2 4 6]));
xlim(sp_ax(end), [1 1.01*max(CF_kHz)]);

ylab_tc_han.Position(1)= -.12;
ylab_rate_han.Position(1)= -.12;

set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len);

% define new axes for AB
Xcorner= .075;
Xwidth= .4;
Xshift= .12;
Ycorner= .078;
Ywidth= .265;
Yshift= .055;

% D
set(sp_ax(5),'Position', [Xcorner Ycorner Xwidth Ywidth])
drawnow

% C
set(sp_ax(3),'Position', [Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% A
set(sp_ax(1),'Position', [Xcorner Ycorner+2.15*Ywidth+2*Yshift Xwidth .85*Ywidth])
drawnow

% F
set(sp_ax(6),'Position', [Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow

% E
set(sp_ax(4),'Position', [Xcorner+Xwidth+Xshift Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% B
set(sp_ax(2),'Position', [Xcorner+Xwidth+Xshift Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

if saveFig
    print([dirStruct.png 'Fig8_fric_SpIN'], '-dpng',  '-r600');
end

%% Create table
all_snr= [fricative_corr_Data.SNR];
all_chinIDs= [fricative_corr_Data.chinID];
xTrack = num2cell([fricative_corr_Data.track]');
xUnit = num2cell([fricative_corr_Data.unit]');
xChinID = num2cell([fricative_corr_Data.chinID]');
UnitID= cellfun(@(x,y,z) sprintf('%d_%d_%d', x, y, z), xChinID, xTrack, xUnit, 'UniformOutput', false);

nonNatInds= ones(size(validCFinds(:)));
valid_fric_inds= nonNatInds(:) & validCFinds(:);

var_tab_chinID= all_chinIDs(valid_fric_inds);
var_tab_unitID= UnitID(valid_fric_inds);
var_tab_CF_kHz= CF_kHz(valid_fric_inds);
var_tab_thresh_dB= thresh_dB(valid_fric_inds);
var_tab_Q10_local= Q10local(valid_fric_inds);
var_tab_Q10_global= Q10global(valid_fric_inds);
var_tab_TTR_dB= TTR_dB(valid_fric_inds);
var_tab_SR= SR(valid_fric_inds);
var_tab_hearingStatus= nhInds(valid_fric_inds);
var_tab_resp_corr= corr_snr_corrected(valid_fric_inds);
var_tab_resp_raw= corr_snr_raw(valid_fric_inds);
var_tab_stim_snr= all_snr(valid_fric_inds);

table_fric_in_noise_withNan= table(var_tab_unitID(:), var_tab_chinID(:), var_tab_CF_kHz(:), var_tab_thresh_dB(:), var_tab_Q10_local(:), var_tab_Q10_global(:), var_tab_TTR_dB(:), ...
    var_tab_SR(:), var_tab_hearingStatus(:), var_tab_stim_snr(:), var_tab_resp_corr(:), var_tab_resp_raw(:),  ...
    'VariableNames', {                  'UnitID',       'ChinID',         'CF_kHz', 'Threshold_dBSPL',              'Q10local',           'Q10global',         'TTR_dB',       ...
    'SpontRate',      'HearingStatus',           'StimSNR',           'RespCorr',          'RespCorrRaw'});

if saveStats
    writetable(table_fric_in_noise_withNan, [dirStruct.stats filesep 'Fig8_table_fric_in_noise_withNan.txt'], 'Delimiter', 'tab');
end