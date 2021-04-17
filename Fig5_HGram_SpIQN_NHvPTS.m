clear;
clc;

saveFig= 1;
saveStats= 1;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];

figSize_cm= [15 5 17 12];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

fs= 20e3; % Sampling frequency for stim, uRate
allSNRs= [inf 0 -5];
nHarms2IncludeNearFormants_F1= 3;
nHarms2IncludeNearFormants_F2_F3= 3;
use_FracSig0_FracHGram1_SynInd2_FracLPpow3= 1;
octRange4Avging= .33;
avgType= 'mean'; % 'weighted-mean' | 'mean' | median

%% anl params
anl.stimPeriod= 1800e-3; % stim rep rate (stim on + stim off)
anl.stimDuration= 1300e-3;
anl.fs= fs; % 1/(bin resolution) for uRate/histcount
anl.nw= 1.5;
anl.tStart= 0; % 0.23;
anl.tEnd= 1.3; %0.66; %1.3;
anl.binRes= 1/anl.fs;
anl.tEdge_hist= anl.tStart:anl.binRes:anl.tEnd;
anl.tBinCenter= (anl.tEdge_hist(1:end-1) + anl.tEdge_hist(2:end))/2;
anl.nHarmonics= 35;
anl.f0lpfBW= 20;
anl.LPFpow_cutoff= 3.5e3;


%% Load stim
[sig, fsSig]= audioread(['mat_data' filesep 'FLN_Stim_S_P.wav']);
sig= gen_resample(sig, fsSig, fs);
danish.voiced_boundaries= find_voicing_boundaries(sig, fs, 0, .13);

temp_f0= load(['mat_data' filesep 'danish_pitch.mat']);
temp_f0= temp_f0.pitch_data;
danish.voiced_inds= any(anl.tBinCenter>danish.voiced_boundaries(:,1) & anl.tBinCenter<danish.voiced_boundaries(:,2), 1);
danish.trajectory.f0= zeros(size(anl.tBinCenter));
danish.trajectory.f0(danish.voiced_inds)= interp1([temp_f0.time], [temp_f0.est], anl.tBinCenter(danish.voiced_inds), 'pchip');
danish.voiced_inds= danish.voiced_inds & danish.trajectory.f0>95 & danish.trajectory.f0<150; % Remove unrealistic values
danish.trajectory.f0(~danish.voiced_inds)= 0;

temp_formants= load(['mat_data' filesep 'danish_formant.mat']);
temp_formants= temp_formants.formant_data;
danish.trajectory.f1= interp1([temp_formants.time], [temp_formants.f1], anl.tBinCenter, 'pchip');
danish.trajectory.f1(~danish.voiced_inds)= nan;
danish.trajectory.f2= interp1([temp_formants.time], [temp_formants.f2], anl.tBinCenter, 'pchip');
danish.trajectory.f2(~danish.voiced_inds)= nan;
danish.trajectory.f3= interp1([temp_formants.time], [temp_formants.f3], anl.tBinCenter, 'pchip');
danish.trajectory.f3(~danish.voiced_inds)= nan;


%% compute trajectory powers
harmMat= (1:anl.nHarmonics)' + zeros(size(danish.trajectory.f0));
f1_harmInds= danish.trajectory.f1./danish.trajectory.f0;
f2_harmInds= danish.trajectory.f2./danish.trajectory.f0;
f3_harmInds= danish.trajectory.f3./danish.trajectory.f0;

if nHarms2IncludeNearFormants_F1==3
    f1_harmInds= [floor(f1_harmInds-.5); round(f1_harmInds); ceil(f1_harmInds+.5)]';
elseif nHarms2IncludeNearFormants_F1==5
    f1_harmInds= [floor(f1_harmInds-1.5); floor(f1_harmInds-.5); round(f1_harmInds); ceil(f1_harmInds+.5); ceil(f1_harmInds+1.5)]';
else
    error('only 3 or 5 (for the time being)');
end

if nHarms2IncludeNearFormants_F2_F3==3
    f2_harmInds= [floor(f2_harmInds-.5); round(f2_harmInds); ceil(f2_harmInds+.5)]';
    f3_harmInds= [floor(f3_harmInds-.5); round(f3_harmInds); ceil(f3_harmInds+.5)]';
elseif nHarms2IncludeNearFormants_F2_F3==5
    f2_harmInds= [floor(f2_harmInds-1.5); floor(f2_harmInds-.5); round(f2_harmInds); ceil(f2_harmInds+.5); ceil(f2_harmInds+1.5)]';
    f3_harmInds= [floor(f3_harmInds-1.5); floor(f3_harmInds-.5); round(f3_harmInds); ceil(f3_harmInds+.5); ceil(f3_harmInds+1.5)]';
else
    error('only 3 or 5 (for the time being)');
end


f1_harmInds(f1_harmInds<=0)= nan;
f2_harmInds(f2_harmInds<=0)= nan;
f3_harmInds(f3_harmInds<=0)= nan;

nPoints= size(f1_harmInds, 1);

f1_harm_mask= zeros(anl.nHarmonics, length(anl.tBinCenter));
f2_harm_mask= zeros(anl.nHarmonics, length(anl.tBinCenter));
f3_harm_mask= zeros(anl.nHarmonics, length(anl.tBinCenter));
nf_harm_mask= zeros(anl.nHarmonics, length(anl.tBinCenter));

temp_f1_subX= f1_harmInds(:);
temp_f1_subX(temp_f1_subX==0)= nan;
temp_f1_subY= repmat((1:nPoints)', nHarms2IncludeNearFormants_F1, 1);
temp_valid_f1inds= find(~(isnan(temp_f1_subX) | isnan(temp_f1_subY)));
temp_f1_inds= sub2ind(size(f1_harm_mask), temp_f1_subX(temp_valid_f1inds), temp_f1_subY(temp_valid_f1inds));
f1_harm_mask(temp_f1_inds)= 1;

temp_f2_subX= f2_harmInds(:);
temp_f2_subX(temp_f2_subX==0)= nan;
temp_f2_subY= repmat((1:nPoints)', nHarms2IncludeNearFormants_F2_F3, 1);
temp_valid_f2inds= find(~(isnan(temp_f2_subX) | isnan(temp_f2_subY)));
temp_f2_inds= sub2ind(size(f2_harm_mask), temp_f2_subX(temp_valid_f2inds), temp_f2_subY(temp_valid_f2inds));
f2_harm_mask(temp_f2_inds)= 1;

temp_f3_subX= f3_harmInds(:);
temp_f3_subX(temp_f3_subX==0)= nan;
temp_f3_subY= repmat((1:nPoints)', nHarms2IncludeNearFormants_F2_F3, 1);
temp_valid_f3inds= find(~(isnan(temp_f3_subX) | isnan(temp_f3_subY)));
temp_f3_inds= sub2ind(size(f3_harm_mask), temp_f3_subX(temp_valid_f3inds), temp_f3_subY(temp_valid_f3inds));
f3_harm_mask(temp_f3_inds)= 1;


%%
dirStruct.loading_dir= ['ANF_Data' filesep];

NHchins= [321 322 325 338 341 343 346 347 354 355 373 379];
HIchins= [358 360 361 362 367 370];
all_chinID= [NHchins, HIchins];

allChinSpikeData = [];
if ~exist('all_chinID', 'var')
    allfiles=dir([dirStruct.loading_dir '*.mat']);
else
    for fileVar=1:length(all_chinID)
        curChinID= all_chinID(fileVar);
        allfiles(fileVar)=dir([dirStruct.loading_dir '*' num2str(curChinID) '*']);
    end
end
for fileVar=1:length(allfiles)
    temp = load([dirStruct.loading_dir allfiles(fileVar).name]);
    allChinSpikeData= [allChinSpikeData; temp.spike_data']; %#ok<AGROW>
end


%%
if use_FracSig0_FracHGram1_SynInd2_FracLPpow3==3
    d_lpf_pow= designfilt('lowpassiir','FilterOrder', 4, ...
        'HalfPowerFrequency', (anl.LPFpow_cutoff)/(anl.fs/2), 'DesignMethod','butter');
else
    d_lpf_pow= nan;
end

%%
figure(1);
clf;

plt.plotVar= true;
plt.octRange4Avging= octRange4Avging;
plt.MINpts= 1;
plt.lw= 1.5;
plt.xtick_freq_val= [500 1e3 2e3 3e3 5e3 10e3];
plt.xtick_freq_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_freq_val/1e3), 'uniformoutput', false);
plt.mrkSize= 4;
spLetters= 'ABCDEFGHI';
spLetters= reshape(reshape(spLetters, 3, 3)', 9, 1)';
plt.tick_len= [0.025 0.025];
count= 0;
txtHan= nan(length(allSNRs)*3, 1);
plt.markerHSR= 'x';
plt.markerLMSR= 'x';
%%

stat_params.cf_min= .5;
stat_params.cf_max= 5;
stat_params.f1_range_kHz= [.5 1];
stat_params.f2_range_kHz= [1 2.5];
stat_params.f3_range_kHz= [2 3];
stat_params.valid_range_vals= nan(3,3);
stat_params.invalid_range_vals= nan(3,3);

%%
Fig5_Fpow_table= [];
for snrVar= 1:length(allSNRs)
    snr2use= allSNRs(snrVar);
    
    noise2use= 'SSN';
    if isinf(snr2use)
        cur_snr_allChinSpikeData= allChinSpikeData(strcmp({allChinSpikeData.noise}, noise2use));
        snr_str= 'Quiet';
    else
        cur_snr_allChinSpikeData= allChinSpikeData([allChinSpikeData.SNR]==snr2use & strcmp({allChinSpikeData.noise}, noise2use));
        snr_str= sprintf('%d dB SNR', snr2use);
    end
    
    [~, uniqInds]= unique([ [cur_snr_allChinSpikeData.chinID]', [cur_snr_allChinSpikeData.track]', [cur_snr_allChinSpikeData.unit]', [cur_snr_allChinSpikeData.SPL]'], 'rows');
    uniqSpikeData= cur_snr_allChinSpikeData(uniqInds);
    
    %% Construct PSTH
    all_SN_d_t= nan(length(uniqSpikeData), length(anl.tBinCenter));
    all_N_d_t= nan(length(uniqSpikeData), length(anl.tBinCenter));
    
    all_SN_d_t_filt= nan(length(uniqSpikeData), length(anl.tBinCenter));
    
    d_lp_f0= designfilt('lowpassiir','FilterOrder', 3, ...
        'HalfPowerFrequency', (anl.f0lpfBW/2)/(anl.fs/2), 'DesignMethod','butter');
    
    all_SN_HGram= cell(length(uniqSpikeData), 1);
    
    F1_strength= nan(length(uniqSpikeData), 1);
    F2_strength= nan(length(uniqSpikeData), 1);
    F3_strength= nan(length(uniqSpikeData), 1);
    
    F1_strength_NF= nan(length(uniqSpikeData), 1);
    F2_strength_NF= nan(length(uniqSpikeData), 1);
    F3_strength_NF= nan(length(uniqSpikeData), 1);
    
    SN_VoicedDrivenRate= nan(length(uniqSpikeData), 1);
    N_VoicedDrivenRate= nan(length(uniqSpikeData), 1);
    
    
    parfor unitVar= 1:length(uniqSpikeData)
        % for unitVar= 10
        cur_unit_data= uniqSpikeData(unitVar);
        cur_CF_Hz= cur_unit_data.CF_Hz;
        
        if cur_unit_data.chinID > 368
            cur_delay= 4.59e-3;
        else
            cur_delay= 0;
        end
        
        if ~isempty(cur_unit_data.SpikeTrains{1,1})
            %% SN
            if isinf(snr2use)
                SpikeTrain_SN_pos= cur_unit_data.SpikeTrains{1,1};
                SpikeTrain_SN_neg= cur_unit_data.SpikeTrains{1,2};
            else
                SpikeTrain_SN_pos= cur_unit_data.SpikeTrains{3,1};
                SpikeTrain_SN_neg= cur_unit_data.SpikeTrains{3,2};
            end
            
            temp_SN_uRate_pos= histcounts( cell2mat(SpikeTrain_SN_pos)-cur_delay, anl.tEdge_hist) / length(SpikeTrain_SN_pos);
            temp_SN_uRate_neg= histcounts( cell2mat(SpikeTrain_SN_neg)-cur_delay, anl.tEdge_hist) / length(SpikeTrain_SN_neg);
            temp_SN_uRate_sum= (temp_SN_uRate_pos + temp_SN_uRate_neg) / 2;
            temp_SN_uRate_diff= (temp_SN_uRate_pos - temp_SN_uRate_neg)/2;
            temp_SN_uRate_diff= temp_SN_uRate_diff .* danish.voiced_inds;
            
            SN_VoicedDrivenRate(unitVar)= sum(temp_SN_uRate_sum .* danish.voiced_inds);
            
            temp_SN_Image= nan(anl.nHarmonics, length(anl.tBinCenter));
            
            if ismember(use_FracSig0_FracHGram1_SynInd2_FracLPpow3, [2, 3])
                for harmVar= 1:anl.nHarmonics
                    [~, temp_SN_Image(harmVar, :)]= get_trajectory_signal(temp_SN_uRate_diff/anl.binRes, anl.fs, harmVar*danish.trajectory.f0, d_lp_f0);
                end
                
                if use_FracSig0_FracHGram1_SynInd2_FracLPpow3==2
                    F1_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f1_harm_mask==1)) .^2, 1))) / SN_VoicedDrivenRate(unitVar);
                    F2_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f2_harm_mask==1)) .^2, 1))) / SN_VoicedDrivenRate(unitVar);
                    F3_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f3_harm_mask==1)) .^2, 1))) / SN_VoicedDrivenRate(unitVar);
                elseif use_FracSig0_FracHGram1_SynInd2_FracLPpow3==3
                    cur_lpf_pow= rms(filter(d_lpf_pow, temp_SN_uRate_diff/anl.binRes));
                    F1_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f1_harm_mask==1)) .^2, 1))) / cur_lpf_pow;
                    F2_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f2_harm_mask==1)) .^2, 1))) / cur_lpf_pow;
                    F3_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f3_harm_mask==1)) .^2, 1))) / cur_lpf_pow;
                end
            elseif ismember(use_FracSig0_FracHGram1_SynInd2_FracLPpow3, [0 1])
                for harmVar= 1:anl.nHarmonics
                    [temp_SN_Image(harmVar, :), ~]= get_trajectory_signal(temp_SN_uRate_diff/anl.binRes, anl.fs, harmVar*danish.trajectory.f0, d_lp_f0);
                end
                
                if use_FracSig0_FracHGram1_SynInd2_FracLPpow3==1
                    F1_strength(unitVar)= nanmean(nanmean( (temp_SN_Image .* (f1_harm_mask==1)) .^2, 1)) / nanmean(nanmean( temp_SN_Image .^2, 1));
                    F2_strength(unitVar)= nanmean(nanmean( (temp_SN_Image .* (f2_harm_mask==1)) .^2, 1)) / nanmean(nanmean( temp_SN_Image .^2, 1));
                    F3_strength(unitVar)= nanmean(nanmean( (temp_SN_Image .* (f3_harm_mask==1)) .^2, 1)) / nanmean(nanmean( temp_SN_Image .^2, 1));
                elseif use_FracSig0_FracHGram1_SynInd2_FracLPpow3==0
                    F1_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f1_harm_mask==1)) .^2, 1)));
                    F2_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f2_harm_mask==1)) .^2, 1)));
                    F3_strength(unitVar)= sqrt(nanmean(nanmean( (temp_SN_Image .* (f3_harm_mask==1)) .^2, 1)));
                end
                
            end
        end
        %% N
        SpikeTrain_N_pos= cur_unit_data.SpikeTrains{2,1};
        SpikeTrain_N_neg= cur_unit_data.SpikeTrains{2,2};
        
        temp_N_uRate_pos= histcounts( cell2mat(SpikeTrain_N_pos)-cur_delay, anl.tEdge_hist) / length(SpikeTrain_N_pos);
        temp_N_uRate_neg= histcounts( cell2mat(SpikeTrain_N_neg)-cur_delay, anl.tEdge_hist) / length(SpikeTrain_N_neg);
        temp_N_uRate_sum= (temp_N_uRate_pos + temp_N_uRate_neg) / 2;
        temp_N_uRate_diff= (temp_N_uRate_pos - temp_N_uRate_neg)/2;

        N_VoicedDrivenRate(unitVar)= sum(temp_N_uRate_sum .* danish.voiced_inds);
        
        temp_N_Image= nan(anl.nHarmonics, length(anl.tBinCenter));
        if ismember(use_FracSig0_FracHGram1_SynInd2_FracLPpow3, [2, 3])
            for harmVar= 1:anl.nHarmonics
                [~, temp_N_Image(harmVar, :)]= get_trajectory_signal(temp_N_uRate_diff/anl.binRes, anl.fs, harmVar*danish.trajectory.f0, d_lp_f0);
            end
            
            if use_FracSig0_FracHGram1_SynInd2_FracLPpow3==2
                F1_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f1_harm_mask==1)) .^2, 1))) / N_VoicedDrivenRate(unitVar);
                F2_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f2_harm_mask==1)) .^2, 1))) / N_VoicedDrivenRate(unitVar);
                F3_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f3_harm_mask==1)) .^2, 1))) / N_VoicedDrivenRate(unitVar);
            elseif use_FracSig0_FracHGram1_SynInd2_FracLPpow3==3
                cur_lpf_pow= rms(filter(d_lpf_pow, temp_N_uRate_diff/anl.binRes));
                F1_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f1_harm_mask==1)) .^2, 1))) / cur_lpf_pow;
                F2_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f2_harm_mask==1)) .^2, 1))) / cur_lpf_pow;
                F3_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f3_harm_mask==1)) .^2, 1))) / cur_lpf_pow;
            end
        elseif ismember(use_FracSig0_FracHGram1_SynInd2_FracLPpow3, [0 1])
            for harmVar= 1:anl.nHarmonics
                [temp_N_Image(harmVar, :), ~]= get_trajectory_signal(temp_N_uRate_diff/anl.binRes, anl.fs, harmVar*danish.trajectory.f0, d_lp_f0);
            end
            
            if use_FracSig0_FracHGram1_SynInd2_FracLPpow3==1
                F1_strength_NF(unitVar)= nanmean(nanmean( (temp_N_Image .* (f1_harm_mask==1)) .^2, 1)) / nanmean(nanmean( temp_N_Image .^2, 1));
                F2_strength_NF(unitVar)= nanmean(nanmean( (temp_N_Image .* (f2_harm_mask==1)) .^2, 1)) / nanmean(nanmean( temp_N_Image .^2, 1));
                F3_strength_NF(unitVar)= nanmean(nanmean( (temp_N_Image .* (f3_harm_mask==1)) .^2, 1)) / nanmean(nanmean( temp_N_Image .^2, 1));
            elseif use_FracSig0_FracHGram1_SynInd2_FracLPpow3==0
                F1_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f1_harm_mask==1)) .^2, 1)));
                F2_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f2_harm_mask==1)) .^2, 1)));
                F3_strength_NF(unitVar)= sqrt(nanmean(nanmean( (temp_N_Image .* (f3_harm_mask==1)) .^2, 1)));
            end
        end
    end
    %% Apply equal weight to CF normalize for CF distribution
    all_cf_Hz= [uniqSpikeData.CF_Hz];
    nhInds= ismember([uniqSpikeData.chinID], NHchins);
    hiInds= ismember([uniqSpikeData.chinID], HIchins);
    
    all_SR= [uniqSpikeData.SR];
    all_chinIDs= [uniqSpikeData.chinID];
    xTrack = num2cell([uniqSpikeData.track]');
    xUnit = num2cell([uniqSpikeData.unit]');
    xChinID = num2cell([uniqSpikeData.chinID]');
    UnitID= cellfun(@(x,y,z) sprintf('%d_%d_%d', x, y, z), xChinID, xTrack, xUnit, 'UniformOutput', false);

    all_thresh_dB= [uniqSpikeData.thresh_dB];
    all_Q10_local= [uniqSpikeData.Q10local];
    all_TTR_dB= [uniqSpikeData.TTR_dB];
    all_hearing= cell(length([uniqSpikeData.chinID]), 1);
    all_hearing(nhInds)= {'NH'};
    all_hearing(hiInds)= {'HI'};
    lowSR= all_SR<18;
    highSR= all_SR>=18;
    valid_inds= (all_SR > -inf) & SN_VoicedDrivenRate(:)'>30;
    
    nhInds= valid_inds & nhInds;
    hiInds= valid_inds & hiInds; % Some of the HI neurons are normal
    
    
    %% Plot 
    % F1
    axTemp= subplot(length(allSNRs), 3, (snrVar-1)*3+1);
    
    if isinf(snr2use)
        F1_strength_rel= F1_strength;
    else
        F1_strength_rel= F1_strength-F1_strength_NF;
    end
    F1_strength_rel(F1_strength_rel<0)= 0;
    
    yyaxis left;
    hold on;
    plot(all_cf_Hz(nhInds & highSR), F1_strength_rel(nhInds & highSR), plt.markerHSR, 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
    plot(all_cf_Hz(nhInds & lowSR), F1_strength_rel(nhInds & lowSR), plt.markerLMSR, 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
    [~, ~, lHan] = helper.octAVG(all_cf_Hz(nhInds), F1_strength_rel(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar, [], stat_params.cf_min);
    set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw, 'color', 'b');
    
    plot(all_cf_Hz(hiInds & highSR), F1_strength_rel(hiInds & highSR), plt.markerHSR, 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
    plot(all_cf_Hz(hiInds & lowSR), F1_strength_rel(hiInds & lowSR), plt.markerLMSR, 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
    [~, ~, lHan] = helper.octAVG(all_cf_Hz(hiInds), F1_strength_rel(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar, [], stat_params.cf_min);
    set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw, 'color', 'r');
    set(gca, 'YColor', 'k', 'TickLength', plt.tick_len);
    ylim([0 0.6])
    
    yyaxis right;
    plot(danish.trajectory.f1, anl.tBinCenter, 'k');
    set(gca, 'XScale', 'log', 'YColor', 'k', 'xtick', plt.xtick_freq_val, 'XTickLabel', plt.xtick_freq_lab);
    
    sp_ax_X((snrVar-1)*3+1)= axTemp;
    sp_ax_leftY((snrVar-1)*3+1)= axTemp.YAxis(1);
    sp_ax_rightY((snrVar-1)*3+1)= axTemp.YAxis(2);
    
    count= count+1;
    txtHan(count)= text(.05, 1.05, sprintf('%s', spLetters(count)), 'Units', 'normalized', 'FontWeight', 'bold');
    
    % F2
    axTemp= subplot(length(allSNRs), 3, (snrVar-1)*3+2);
    hold on;
    
    F2_strength_rel= F2_strength-F2_strength_NF;
    F2_strength_rel(F2_strength_rel<0)= 0;
    
    yyaxis left;
    plot(all_cf_Hz(nhInds & highSR), F2_strength_rel(nhInds & highSR), plt.markerHSR, 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
    plot(all_cf_Hz(nhInds & lowSR), F2_strength_rel(nhInds & lowSR), plt.markerLMSR, 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
    [~, ~, lHan] = helper.octAVG(all_cf_Hz(nhInds), F2_strength_rel(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar, [], stat_params.cf_min);
    set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw, 'color', 'b');
    
    plot(all_cf_Hz(hiInds & highSR), F2_strength_rel(hiInds & highSR), plt.markerHSR, 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
    plot(all_cf_Hz(hiInds & lowSR), F2_strength_rel(hiInds & lowSR), plt.markerLMSR, 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
    [~, ~, lHan] = helper.octAVG(all_cf_Hz(hiInds), F2_strength_rel(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar, [], stat_params.cf_min);
    set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw, 'color', 'r');
    set(gca, 'YColor', 'k', 'TickLength', plt.tick_len);
    ylim([0 0.175])
    
    yyaxis right;
    plot(danish.trajectory.f2, anl.tBinCenter, 'k');
    set(gca, 'XScale', 'log', 'YColor', 'k', 'xtick', plt.xtick_freq_val, 'XTickLabel', plt.xtick_freq_lab);
    
    sp_ax_X((snrVar-1)*3+2)= axTemp;
    sp_ax_leftY((snrVar-1)*3+2)= axTemp.YAxis(1);
    sp_ax_rightY((snrVar-1)*3+2)= axTemp.YAxis(2);
    
    count= count+1;
    txtHan(count)= text(.05, 1.05, sprintf('%s', spLetters(count)), 'Units', 'normalized', 'FontWeight', 'bold');
    
    % F3
    axTemp= subplot(length(allSNRs), 3, (snrVar-1)*3+3);
    hold on;
    
    F3_strength_rel= F3_strength-F3_strength_NF;
    F3_strength_rel(F3_strength_rel<0)=0 ;
    
    yyaxis left;
    plot(all_cf_Hz(nhInds & highSR), F3_strength_rel(nhInds & highSR), plt.markerHSR, 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
    plot(all_cf_Hz(nhInds & lowSR), F3_strength_rel(nhInds & lowSR), plt.markerLMSR, 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
    [~, ~, lHan] = helper.octAVG(all_cf_Hz(nhInds), F3_strength_rel(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar, [], stat_params.cf_min);
    set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw, 'color', 'b');
    
    plot(all_cf_Hz(hiInds & highSR), F3_strength_rel(hiInds & highSR), plt.markerHSR, 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
    plot(all_cf_Hz(hiInds & lowSR), F3_strength_rel(hiInds & lowSR), plt.markerLMSR, 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
    [~, ~, lHan] = helper.octAVG(all_cf_Hz(hiInds), F3_strength_rel(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar, [], stat_params.cf_min);
    set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw, 'color', 'r');
    set(gca, 'YColor', 'k', 'TickLength', plt.tick_len);
    
    ylim([0 0.075])
    
    yyaxis right;
    plot(danish.trajectory.f3, anl.tBinCenter, 'k');
    set(gca, 'XScale', 'log', 'YColor', 'k', 'xtick', plt.xtick_freq_val, 'XTickLabel', plt.xtick_freq_lab);
    
    sp_ax_X((snrVar-1)*3+3)= axTemp; %#ok<*SAGROW>
    sp_ax_leftY((snrVar-1)*3+3)= axTemp.YAxis(1);
    sp_ax_rightY((snrVar-1)*3+3)= axTemp.YAxis(2);
    
    count= count+1;
    txtHan(count)= text(.05, 1.05, sprintf('%s', spLetters(count)), 'Units', 'normalized', 'FontWeight', 'bold');
    
    %% stats
    
    allCFs_kHz= all_cf_Hz/1e3;
    cfs_in_F1= (allCFs_kHz>stat_params.f1_range_kHz(1)) & (allCFs_kHz<stat_params.f1_range_kHz(2));
    cfs_in_F2= (allCFs_kHz>stat_params.f2_range_kHz(1)) & (allCFs_kHz<stat_params.f2_range_kHz(2));
    cfs_in_F3= (allCFs_kHz>stat_params.f3_range_kHz(1)) & (allCFs_kHz<stat_params.f3_range_kHz(2));
    
    cfs_out_F1= ( (allCFs_kHz>stat_params.cf_min) & (allCFs_kHz<stat_params.f1_range_kHz(1)) ) |  ( (allCFs_kHz<stat_params.cf_max) & (allCFs_kHz>stat_params.f1_range_kHz(2)) );
    cfs_out_F2= ( (allCFs_kHz>stat_params.cf_min) & (allCFs_kHz<stat_params.f2_range_kHz(1)) ) |  ( (allCFs_kHz<stat_params.cf_max) & (allCFs_kHz>stat_params.f2_range_kHz(2)) );
    cfs_out_F3= ( (allCFs_kHz>stat_params.cf_min) & (allCFs_kHz<stat_params.f3_range_kHz(1)) ) |  ( (allCFs_kHz<stat_params.cf_max) & (allCFs_kHz>stat_params.f3_range_kHz(2)) );
    
    [~, stat_params.valid_range_vals(snrVar, 1)]= ttest2(F1_strength_rel(nhInds & cfs_in_F1), F1_strength_rel(hiInds & cfs_in_F1));
    [~, stat_params.valid_range_vals(snrVar, 2)]= ttest2(F2_strength_rel(nhInds & cfs_in_F2), F2_strength_rel(hiInds & cfs_in_F2));
    [~, stat_params.valid_range_vals(snrVar, 3)]= ttest2(F3_strength_rel(nhInds & cfs_in_F3), F3_strength_rel(hiInds & cfs_in_F3));
    
    [~, stat_params.invalid_range_vals(snrVar, 1)]= ttest2(F1_strength_rel(nhInds & cfs_out_F1), F1_strength_rel(hiInds & cfs_out_F1));
    [~, stat_params.invalid_range_vals(snrVar, 2)]= ttest2(F2_strength_rel(nhInds & cfs_out_F2), F2_strength_rel(hiInds & cfs_out_F2));
    [~, stat_params.invalid_range_vals(snrVar, 3)]= ttest2(F3_strength_rel(nhInds & cfs_out_F3), F3_strength_rel(hiInds & cfs_out_F3));

    var_snr= repmat(snr2use, numel(allCFs_kHz), 1);
    cur_snr_tbl= table(all_chinIDs(:), UnitID(:), allCFs_kHz(:), all_SR(:), all_hearing(:), all_thresh_dB(:), all_Q10_local(:), all_TTR_dB(:), ...
    var_snr(:), F1_strength_rel(:), F2_strength_rel(:), F3_strength_rel(:), ...
    'VariableNames', {  'ChinID',      'UnitID',   'CF_kHz',    'SpontRate',  'HearingStatus',      'thresh_dB',        'Q10local',      'TTR_dB', ...
      'SNR',        'F1pow',            'F2pow',            'F3pow'});
    Fig5_Fpow_table= [Fig5_Fpow_table; cur_snr_tbl]; %#ok<AGROW>

end

if saveStats
    writetable(Fig5_Fpow_table, [dirStruct.stats 'Fig5_Fpow_table.txt'], 'Delimiter', 'tab');
end

subplot(length(allSNRs), 3, 1);
ttlHan(1)= text(-.35, .25, 'Quiet', 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Rotation', 90);
subplot(length(allSNRs), 3, 4);
ttlHan(2)= text(-.35, 0.09, '0 dB SNR', 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Rotation', 90);
subplot(length(allSNRs), 3, 7);
ttlHan(3)= text(-.35, .08, '-5 dB SNR', 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Rotation', 90);

subplot(length(allSNRs), 3, 1);
ttlHan(4)= text(0.5, 1.075, 'F_1', 'Units', 'normalized', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
subplot(length(allSNRs), 3, 2);
ttlHan(5)= text(0.5, 1.075, 'F_2', 'Units', 'normalized', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
subplot(length(allSNRs), 3, 3);
ttlHan(6)= text(0.5, 1.075, 'F_3', 'Units', 'normalized', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');


%%
if any(isinf(allSNRs)) && ismember(use_FracSig0_FracHGram1_SynInd2_FracLPpow3, [0 1 2 3])
    linkprop(sp_ax_leftY(4:3:end), 'limits');
    linkprop(sp_ax_leftY(5:3:end), 'limits');
    linkprop(sp_ax_leftY(6:3:end), 'limits');
else
    linkprop(sp_ax_leftY(1:3:end), 'limits');
    linkprop(sp_ax_leftY(2:3:end), 'limits');
    linkprop(sp_ax_leftY(3:3:end), 'limits');
end
linkprop(sp_ax_rightY, 'limits');
sp_ax_rightY(end).Limits=[0 1.3];

linkaxes(sp_ax_X([1 4 7]), 'x');
xlim(sp_ax_X(1), [150 5.5e3]);
linkaxes(sp_ax_X([2 3 5 6 8 9]), 'x');
xlim(sp_ax_X(9), [500 5.5e3]);

axes(sp_ax_X(4));
yyaxis left
ylabel('Formant power');

axes(sp_ax_X(6));
yyaxis right;
ylabel('Time (s)');

axes(sp_ax_X(end-1));
xlabel('Characteristic Frequency or Frequency (kHz)');

axes(sp_ax_X(end));
lHan= plot(nan, nan, '-b', nan, nan, '-r', nan, nan, '-k');
set(lHan(1:2), 'linew', plt.lw);

[legHan, icons]= legend(lHan, 'NH', 'HI', 'Traj.', 'Location', 'northwest', 'box', 'off');
icons(4).XData= mean(icons(4).XData) + [0 +.25];
icons(6).XData= mean(icons(6).XData) + [0 +.25];
icons(8).XData= mean(icons(8).XData) + [0 +.25];
legHan.Position(1:2)= [.7 .2];



set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(legHan,'FontSize', 8);
set(txtHan,'FontSize', 10);
set(ttlHan,'FontSize', 13);
%% define new axes for AB
Xcorner= .09;
Xwidth= .215;
Xshift= .1;
Ycorner= .08;
Ywidth= .24;
Yshift= .08;

% G
set(sp_ax_X(7),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow

% H
set(sp_ax_X(8),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow

% I
set(sp_ax_X(9),'Position',[Xcorner+2*Xwidth+2*Xshift Ycorner Xwidth Ywidth])
drawnow

% D
set(sp_ax_X(4),'Position',[Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% E
set(sp_ax_X(5),'Position',[Xcorner+Xwidth+Xshift Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% F
set(sp_ax_X(6),'Position',[Xcorner+2*Xwidth+2*Xshift Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% A
set(sp_ax_X(1),'Position',[Xcorner Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

% B
set(sp_ax_X(2),'Position',[Xcorner+Xwidth+Xshift Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

% C
set(sp_ax_X(3),'Position',[Xcorner+2*Xwidth+2*Xshift Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

if (use_FracSig0_FracHGram1_SynInd2_FracLPpow3==1) && saveFig
    print([dirStruct.png 'Fig5_HGram_SpIQN_NHvPTS'], '-dpng',  '-r600');
end

