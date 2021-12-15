% Look at modulation filtered meanrates.
clear;
clc;

saveFig= 0;
saveStats= 0;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];
dirStruct.loading_dir= ['ANF_Data' filesep];

chinIDs= [321 322 325 338 341 343 346 347 354 355 358 360 361 362 367 370 373 379];
HIchins= [358 360 361 362 367 370];

figSize_cm= [-10 5 8.5 8.5];
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
    
    
    %     [corrStruct, rmsStruct]= helper.uRate_lpf_corr_S_N_SN(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
    %         SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, anl.modFreq, anl.N_lp, anl.tStart, anl.tEnd);
    
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
    
    
    %     fricative_corr_Data(iterVar).s_sn= (corrStruct.s_sn.pos + corrStruct.s_sn.neg)/2;
    %     fricative_corr_Data(iterVar).s_n= (corrStruct.s_n.pos + corrStruct.s_n.neg)/2;
    %     fricative_corr_Data(iterVar).sn_n= (corrStruct.sn_n.pos + corrStruct.sn_n.neg)/2;
    %
    %     fricative_corr_Data(iterVar).rms_S= sqrt(rmsStruct.s.pos*rmsStruct.s.neg);
    %     fricative_corr_Data(iterVar).rms_SN= sqrt(rmsStruct.sn.pos*rmsStruct.sn.neg);
    %     fricative_corr_Data(iterVar).rms_N= sqrt(rmsStruct.n.pos*rmsStruct.n.neg);
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
% rms_s= [fricative_corr_Data.rms_S];
% rms_sn= [fricative_corr_Data.rms_SN];

Rate_SN_Hz= [fricative_corr_Data.DR_SN];
Rate_N_Hz= [fricative_corr_Data.DR_N];
Rate_SNrelN= max(1, Rate_SN_Hz - Rate_N_Hz);
% Rate_SNrelN= max(-inf, Rate_SN_Hz - Rate_N_Hz);
% warning('remove this monstrosity')

snrs2use= sort(unique([fricative_corr_Data.SNR]), 'descend');

SRboundary= 18;
hsrInds= SR>SRboundary;
lmsrInds= SR<SRboundary;


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
plt.lw0= .5;
plt.lw1= 0.85;
plt.lw2= 1;
plt.leg_fSize= 7;
plt.fSize= 9;
plt.ttl_fSize= 10;
plt.xtick_vals_cf_kHz= [.5 1 2 4 8];
plt.rate_ytick= [0 100 200];
plt.mrkSize= 3;
plt.mrkSizeLeg= 4;
plt.tick_len= [.035 .015];
plt.txtX= -0.44;
plt.txtY= 0.0;

figure(1);
clf;
sp_ax= nan(length(snrs2use)*2, 1);
ttlHan= nan(5, 1); % 3 SNRs + SN + N 

for snrVar= 1:length(snrs2use)
    curSNR= snrs2use(snrVar);
    curSNRinds= [fricative_corr_Data.SNR]==curSNR;
    validCFinds= (CF_kHz>3) & (CF_kHz<8);
    
    %% Driven rate: Noisy speech
    sp_num= 2*snrVar-1;
    sp_ax(sp_num)= subplot(length(snrs2use), 2, sp_num);
    hold on;

    plot(CF_kHz(hsrInds & nhInds & curSNRinds), Rate_SN_Hz(hsrInds & nhInds & curSNRinds), '*', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    plot(CF_kHz(lmsrInds & nhInds & curSNRinds), Rate_SN_Hz(lmsrInds & nhInds & curSNRinds), 's', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    plot(CF_kHz(hsrInds & hiInds & curSNRinds), Rate_SN_Hz(hsrInds & hiInds & curSNRinds), '*', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    plot(CF_kHz(lmsrInds & hiInds & curSNRinds), Rate_SN_Hz(lmsrInds & hiInds & curSNRinds), 's', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);

    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & nhInds & curSNRinds), Rate_SN_Hz(hsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');
    
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & nhInds & curSNRinds), Rate_SN_Hz(lmsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');
    
    [~, ~, lHan_hi_hsr] = helper.octAVG(CF_kHz(hsrInds & hiInds & curSNRinds), Rate_SN_Hz(hsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_hi_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');
    
    [~, ~, lHan_hi_lmsr] = helper.octAVG(CF_kHz(lmsrInds & hiInds & curSNRinds), Rate_SN_Hz(lmsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_hi_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');
        
    set(gca, 'XScale', 'log', 'XTick', plt.xtick_vals_cf_kHz, 'YTick', plt.rate_ytick, 'TickDir', 'out','FontSize', plt.fSize);
    if snrVar<length(snrs2use)
        set(gca, 'XTickLabel', '');
    end
    
    %% Driven rate: Noise alone
    sp_num= 2*snrVar;
    sp_ax(sp_num)= subplot(length(snrs2use), 2, sp_num);
    hold on;
    plot(CF_kHz(hsrInds & nhInds & curSNRinds), Rate_N_Hz(hsrInds & nhInds & curSNRinds), '*', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    plot(CF_kHz(lmsrInds & nhInds & curSNRinds), Rate_N_Hz(lmsrInds & nhInds & curSNRinds), 's', 'Color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    plot(CF_kHz(hsrInds & hiInds & curSNRinds), Rate_N_Hz(hsrInds & hiInds & curSNRinds), '*', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);
    plot(CF_kHz(lmsrInds & hiInds & curSNRinds), Rate_N_Hz(lmsrInds & hiInds & curSNRinds), 's', 'Color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw0);

    [~, ~, lHan_nh_hsr] = helper.octAVG(CF_kHz(hsrInds & nhInds & curSNRinds), Rate_N_Hz(hsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');
    
    [~, ~, lHan_nh_lmsr] = helper.octAVG(CF_kHz(lmsrInds & nhInds & curSNRinds), Rate_N_Hz(lmsrInds & nhInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_nh_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');
    
    [~, ~, lHan_hi_hsr] = helper.octAVG(CF_kHz(hsrInds & hiInds & curSNRinds), Rate_N_Hz(hsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_hi_hsr, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');
    
    [~, ~, lHan_hi_lmsr] = helper.octAVG(CF_kHz(lmsrInds & hiInds & curSNRinds), Rate_N_Hz(lmsrInds & hiInds & curSNRinds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan_hi_lmsr, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');
    
    set(gca, 'XScale', 'log', 'XTick', plt.xtick_vals_cf_kHz, 'YTick', plt.rate_ytick, 'YTickLabel', '', 'TickDir', 'out','FontSize', plt.fSize);
    if snrVar<length(snrs2use)
        set(gca, 'XTickLabel', '');
    end
    
    
end

%%
linkaxes(sp_ax(3:6));
ylim(sp_ax(3), [-1 201])

linkaxes(sp_ax(1:2));
ylim(sp_ax(1), [-1 251])

axes(sp_ax(5))
xlab_han= xlabel('Characteristic frequency (kHz)', 'Units', 'normalized');
xlab_han.Position(1:2)= [1.15, -.265];

axes(sp_ax(1));
doLeg= 1;
if doLeg
    lHan(1)= plot(nan, nan, '*--', 'Color', 'b', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
    lHan(2)= plot(nan, nan, 's-', 'Color', 'b', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
    lHan(3)= plot(nan, nan, '*--', 'Color', 'r', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
    lHan(4)= plot(nan, nan, 's-', 'Color', 'r', 'LineWidth', plt.lw0, 'markersize', plt.mrkSizeLeg);
    
    legend_str= {'HSR', 'LMSR', 'HSR', 'LMSR'};
    
    [legHan, icons]= legend(lHan, legend_str, 'box', 'off', 'Location', 'south','Orientation','horizontal');
    legHan.Position(1:2)= [.1 .91];
    icons(1).FontSize= plt.leg_fSize;
    icons(2).FontSize= plt.leg_fSize;
    icons(3).FontSize= plt.leg_fSize;
    icons(4).FontSize= plt.leg_fSize;
    icons(1).Position(1)= icons(1).Position(1) - .03;
    icons(2).Position(1)= icons(2).Position(1) - .03;
    icons(3).Position(1)= icons(3).Position(1) - .03;
    icons(4).Position(1)= icons(4).Position(1) - .03;
    
    count= 5;
    icons(count).XData= mean(icons(count).XData) + [-.03, .04]; count= count+2;
    icons(count).XData= mean(icons(count).XData) + [-.03, .03]; count= count+2;
    icons(count).XData= mean(icons(count).XData) + [-.03, .04]; count= count+2;
    icons(count).XData= mean(icons(count).XData) + [-.03, .03]; % count= count+1;
end

%%

axes(sp_ax(1));
ttlHan(1)= text(plt.txtX, plt.txtY, sprintf('%d dB SNR', snrs2use(1)), 'Units', 'normalized', 'VerticalAlignment', 'middle', 'Rotation', 90, 'FontWeight', 'bold');
ttlHan(4)= text(.5, 1.25, 'Noisy speech', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontWeight', 'bold');

axes(sp_ax(2));
ttlHan(5)= text(.5, 1.25, 'Noise alone', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontWeight', 'bold');

axes(sp_ax(3));
ylab_rate_han= ylabel('Rate (spikes/s)', 'Units', 'normalized');
ttlHan(2)= text(plt.txtX, plt.txtY, sprintf('%d dB SNR', snrs2use(2)), 'Units', 'normalized', 'VerticalAlignment', 'middle', 'Rotation', 90, 'FontWeight', 'bold');

axes(sp_ax(5));
ttlHan(3)= text(plt.txtX, plt.txtY, sprintf('%d dB SNR', snrs2use(3)), 'Units', 'normalized', 'VerticalAlignment', 'middle', 'Rotation', 90, 'FontWeight', 'bold');

set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len);
set(ttlHan,'FontSize', plt.ttl_fSize);

%% define new axes for AB
Xcorner= .18;
Xwidth= .375;
Xshift= .05;
Ycorner= .12;
Ywidth= .235;
Yshift= .05;

set(sp_ax(5),'Position',[Xcorner, Ycorner, Xwidth, Ywidth])
drawnow

set(sp_ax(6),'Position',[Xcorner+Xwidth+Xshift, Ycorner, Xwidth, Ywidth])
drawnow

set(sp_ax(3),'Position',[Xcorner, Ycorner+Ywidth+1*Yshift, Xwidth, Ywidth])
drawnow

set(sp_ax(4),'Position',[Xcorner+Xwidth+Xshift, Ycorner+Ywidth+1*Yshift, Xwidth, Ywidth])
drawnow

set(sp_ax(1),'Position',[Xcorner, Ycorner+2*Ywidth+2*Yshift, Xwidth, Ywidth])
drawnow

set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift, Ycorner+2*Ywidth+2*Yshift, Xwidth, Ywidth])
drawnow


if saveFig
    print([dirStruct.png 'FigS3_fric_DR_SN_N'], '-dpng',  '-r600');
    saveas(gcf, [dirStruct.eps 'FigS3_fric_DR_SN_N'], 'epsc');
end

