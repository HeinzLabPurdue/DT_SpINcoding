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


%% Load all AN data for Q374

anl.NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
anl.HIchins= [358 360 361 362 367 370];
demo.NH_chin_track_unit= [343, 1, 5]; % One more NH example = [322, 3, 2]
demo.HI_chin_track_unit= [361, 1, 3]; % HI example = [361, 1, 3]
all_chinDanishData= helper.load_danish_chin_data(anl);
nhChinInd= find(ismember([ [all_chinDanishData.chinID]', [all_chinDanishData.track]', [all_chinDanishData.unit]' ], demo.NH_chin_track_unit, 'rows'));
hiChinInd= find(ismember([ [all_chinDanishData.chinID]', [all_chinDanishData.track]', [all_chinDanishData.unit]' ], demo.HI_chin_track_unit, 'rows'));


%% Examples showing loss of onset

% Analysis parameters for /s/
fs= 5e3;
anl.tStart= 745e-3;
anl.tEnd= 835e-3;
anl.tZoom= 25e-3;
anl.winHardOnset= 10e-3;
anl.winSoftOnset= 00e-3;
anl.winSustained= 80e-3;
anl.tBinEdges= (anl.tStart-anl.tZoom):(1/fs):(anl.tEnd+anl.tZoom);
anl.tBinCenters= ( anl.tBinEdges(1:end-1) + anl.tBinEdges(2:end) )/2;
anl.tValidWindow= (anl.tBinCenters>anl.tStart) & (anl.tBinCenters<anl.tEnd);
anl.nw= 5;

demo.nh_data= all_chinDanishData(nhChinInd);
demo.hi_data= all_chinDanishData(hiChinInd);
demo.delay= 0; % Both 343 and 361 were without invFIR filtering

demo.nh_uRate_pos= histcounts( cell2mat(demo.nh_data.pos) - demo.delay, anl.tBinEdges) / length(demo.nh_data.pos);
demo.nh_uRate_neg= histcounts( cell2mat(demo.nh_data.neg) - demo.delay, anl.tBinEdges) / length(demo.nh_data.neg);
demo.hi_uRate_pos= histcounts( cell2mat(demo.hi_data.pos) - demo.delay, anl.tBinEdges) / length(demo.hi_data.pos);
demo.hi_uRate_neg= histcounts( cell2mat(demo.hi_data.neg) - demo.delay, anl.tBinEdges) / length(demo.hi_data.neg);

%% Read signal
[sig.data, sig.fs]= audioread(['mat_data' filesep 'FLN_Stim_S_P.wav']);
sig.data= gen_resample(sig.data, sig.fs, 20e3);
sig.data= gen_rescale(sig.data, 70);
sig.fs= 20e3;
sig.time= (1:length(sig.data))/sig.fs;
sig.fric_inds= (sig.time> anl.tStart) & (sig.time<anl.tEnd);
sig.fric_inds_zoom= (sig.time> (anl.tStart-anl.tZoom)) & (sig.time<(anl.tEnd+anl.tZoom));
sig.tWindow= sig.time(sig.fric_inds);

demo.latency_correction= 0e-3;

timeStruct.Onset.Start= anl.tStart;
timeStruct.Onset.HardEnd= timeStruct.Onset.Start + anl.winHardOnset; % Using values from Delgutte and Kiang 1984 (Fricative coding)
timeStruct.Onset.SoftEnd= timeStruct.Onset.HardEnd + anl.winSoftOnset;

timeStruct.Sustained.Start= anl.tEnd-anl.winSustained;
timeStruct.Sustained.End= anl.tEnd;

OnsetInds.Hard= (anl.tBinCenters> (timeStruct.Onset.Start+demo.latency_correction)) & (anl.tBinCenters<(timeStruct.Onset.HardEnd+demo.latency_correction)); % 10 ms in Delgutte's paper
OnsetInds.Soft= 0*OnsetInds.Hard;
temp_soft_start= find(anl.tBinCenters>=(timeStruct.Onset.HardEnd+demo.latency_correction), 1);
temp_soft_end= find(anl.tBinCenters<timeStruct.Onset.SoftEnd+demo.latency_correction, 1, 'last');
OnsetInds.Soft(temp_soft_start:temp_soft_end)= linspace(1, 0, temp_soft_end - temp_soft_start + 1);

anl.OnsetMask= OnsetInds.Hard + OnsetInds.Soft;
anl.SustainedMask= (anl.tBinCenters > (timeStruct.Sustained.Start+demo.latency_correction)) & (anl.tBinCenters < (timeStruct.Sustained.End+demo.latency_correction)) ; % 10 ms in Delgutte's paper

%% Plot
plt.lw1= 1;
plt.lw2= 1.5;
plt.mrkSize= 4;
plt.plotVar= true;
plt.octRange4Avging= .33;
plt.avgType= 'mean';
plt.MINpts= 1;
plt.xtick_freq_val= [1 2 4 8 10];
plt.xtick_freq_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_freq_val), 'uniformoutput', false);
plt.tick_len= [.02 .02];
plt.ytick_psd= -90:10:-60;

Ashift= 1.2*max(demo.nh_uRate_neg)+max(demo.hi_uRate_pos);
yTickScale= .2;
ytickValsDefault= [-yTickScale 0 yTickScale];
ytickVals= [ytickValsDefault Ashift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault, 1, 2)), 'UniformOutput', false);

clf;

sp_ax(1)= subplot(3, 2, 1);
yyaxis right;
hold on;
[Pxx,Freq_xx,lHan]= plot_dpss_psd(sig.data(sig.fric_inds), sig.fs, 'nw', anl.nw, 'xunit', 'khz');
set(lHan, 'color', get_color('gray'), 'linew', plt.lw1);
set(gca, 'YColor', 'k', 'YTick', plt.ytick_psd);
xlabel('');
ylabel('PSD (dB/Hz)');

yyaxis left;
plot(demo.nh_data.TC.freqkHz, demo.nh_data.TC.TCfit, '-', 'Color', 'b', 'linew', plt.lw2);
plot(demo.hi_data.TC.freqkHz, demo.hi_data.TC.TCfit, '-', 'Color', 'r', 'linew', plt.lw2);
set(gca, 'YColor', get_color('k'), 'XTick', plt.xtick_freq_val);
ylim([0 120]);
xlim([.75 10]);
ylabel('TC Thresh. (dB SPL)');
txtHan(1)= text(.05, .95, 'A', 'Units', 'normalized', 'FontWeight', 'bold');

sp_ax(2)= subplot(3, 2, 2);
hold on;
plot(anl.tBinCenters, Ashift+demo.nh_uRate_pos, 'Color', get_color('b'), 'linew', plt.lw1);
plot(anl.tBinCenters, Ashift+-demo.nh_uRate_neg, 'Color', get_color('lb'), 'linew', plt.lw1);
plot(anl.tBinCenters, demo.hi_uRate_pos, 'Color', get_color('r'), 'linew', plt.lw1);
plot(anl.tBinCenters, -demo.hi_uRate_neg, 'Color', get_color('lr'), 'linew', plt.lw1);
plot(sig.time(sig.fric_inds_zoom), -Ashift+sig.data(sig.fric_inds_zoom), 'Color', get_color('gray'), 'linew', plt.lw1);
ylabel('Rate (spikes/s)');
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs);
ylim([-0.7 .85]);

yyaxis right;
plot(anl.tBinCenters, anl.OnsetMask, '-', 'Color', get_color('c'), 'linew', plt.lw2);
plot(anl.tBinCenters, anl.SustainedMask, '-', 'Color', get_color('m'), 'linew', plt.lw2);
ylim([-0.075 5]);
set(gca, 'YTick', [0 1], 'YColor', 'k');
xlim([anl.tStart-anl.tZoom anl.tEnd+anl.tZoom]);
box off;
txtHan(2)= text(.05, .95, 'B', 'Units', 'normalized', 'FontWeight', 'bold');
ylabel('Onset & Sustained mask');

%% Onset and sustained rate calculation
onsetRates= nan(length(all_chinDanishData), 1);
sustainedRates= nan(length(all_chinDanishData), 1);

parfor unitVar= 1:length(all_chinDanishData)
    cur_data= all_chinDanishData(unitVar);
    if ~isnan(cur_data.CF_Hz)
        uR_pos= histcounts( cell2mat(cur_data.pos) - cur_data.delay, anl.tBinEdges) / length(cur_data.pos);
        uR_neg= histcounts( cell2mat(cur_data.neg) - cur_data.delay, anl.tBinEdges) / length(cur_data.neg);
        onsetRates(unitVar)= ( sum(uR_pos.*anl.OnsetMask) + sum(uR_neg.*anl.OnsetMask) ) / 2 / (anl.winHardOnset + .5*anl.winSoftOnset);
        sustainedRates(unitVar)= ( sum(uR_pos.*anl.SustainedMask) + sum(uR_neg.*anl.SustainedMask) ) / 2 / anl.winSustained;
    end
end

all_cf_kHz= [all_chinDanishData.CF_Hz]/1e3;
all_chinIDs= [all_chinDanishData.chinID];
all_SR= [all_chinDanishData.SR];
all_thresh_dB= [all_chinDanishData.thresh_dB];
all_Q10_local= [all_chinDanishData.Q10local];
all_TTR_dB= [all_chinDanishData.TTR_dB];
nhInds= ismember([all_chinDanishData.chinID], anl.NHchins)';
hiInds= ismember([all_chinDanishData.chinID], anl.HIchins)';
all_hearing= cell(length([all_chinDanishData.chinID]), 1);
all_hearing(nhInds)= {'NH'};
all_hearing(hiInds)= {'HI'};



cfRange_min_kHz= 3;
cfRange_max_kHz= 8;
valid_cf= (all_cf_kHz(:)>cfRange_min_kHz) & (all_cf_kHz(:)<cfRange_max_kHz);
[~, pOnset]= ttest2(onsetRates(nhInds&valid_cf), onsetRates(hiInds&valid_cf));
fprintf('Onset p=%.3f: %.1f < CF < %.1f kHz\n', pOnset, cfRange_min_kHz, cfRange_max_kHz);

[~, pSustained]= ttest2(sustainedRates(nhInds&valid_cf), sustainedRates(hiInds&valid_cf));
fprintf('Sustained p=%.3f: %.1f < CF < %.1f kHz\n', pOnset, cfRange_min_kHz, cfRange_max_kHz);

sp_ax(3)= subplot(3, 2, 3);
ylabel('PSD (dB/Hz)');

hold on;

plot(all_cf_kHz(nhInds), onsetRates(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), onsetRates(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');

[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');

plot(all_cf_kHz(hiInds), onsetRates(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), onsetRates(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');

[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');


set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
ylabel('Onset rate (spikes/s)');
txtHan(3)= text(.05, .95, 'C', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

sp_ax(5)= subplot(3, 2, 5);

hold on;

plot(all_cf_kHz(nhInds), sustainedRates(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), sustainedRates(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'b');

[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'b');

plot(all_cf_kHz(hiInds), sustainedRates(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), sustainedRates(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw2, 'color', 'r');

[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw1, 'color', 'r');

[~, pSus]= ttest2(sustainedRates(nhInds&valid_cf), sustainedRates(hiInds&valid_cf));
fprintf('Sustained p=%.3f: %.1f < CF < %.1f kHz\n', pSus, cfRange_min_kHz, cfRange_max_kHz);



set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
xlabel('Characteristic Frequency or Frequency (kHz)');
ylabel('Sustained rate (spikes/s)');
txtHan(4)= text(.05, .95, 'D', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

legLine(1)= plot(nan, nan, '-', 'color', 'b', 'linew', plt.lw2);
legLine(2)= plot(nan, nan, '--', 'color', 'b', 'linew', plt.lw1);
legLine(3)= plot(nan, nan, '-', 'color', 'r', 'linew', plt.lw2);
legLine(4)= plot(nan, nan, '--', 'color', 'r', 'linew', plt.lw1);
[legHan, icons]= legend(legLine, 'NH (Driven)', 'NH (Spont.)', 'HI (Driven)', 'HI (Spont.)', 'location', 'northwest', 'box', 'off');
legHan.FontSize= 9;
icons(5).XData= mean(icons(5).XData) + [-0.1 +.15];
icons(7).XData= mean(icons(7).XData) + [-0.1 +.15];
icons(9).XData= mean(icons(9).XData) + [-0.1 +.15];
icons(11).XData= mean(icons(11).XData) + [-0.1 +.15];
legHan.Position(1:2)= [.1 .23];


linkaxes(sp_ax([1 3 5]), 'x');
xlim(sp_ax(1), [.75 9]);

Fig6_tbl_s_coding_AN= table(all_chinIDs(:), all_cf_kHz(:), all_SR(:), all_hearing(:), all_thresh_dB(:), all_Q10_local(:), all_TTR_dB(:), sustainedRates(:), onsetRates(:), ...
    'VariableNames', {          'ChinID',    'CF_kHz',   'SpontRate', 'HearingStatus',  'thresh_dB',      'Q10local',        'TTR_dB',      'SusRate',          'OnRate'});

if saveStats
    writetable(Fig6_tbl_s_coding_AN, [dirStruct.stats 'Fig6_AN_s_all.txt'], 'Delimiter', 'tab');
end

%% FFR data
RootDataDir= ['FFR_Data' filesep];
allDirs= dir([RootDataDir '*SFR*']);

FFR_onset= nan(length(allDirs), 1);
FFR_onset_nf= nan(length(allDirs), 1);

nhDirInd= contains({allDirs.name}', 'NH');
hiDirInd= contains({allDirs.name}', 'PTS');
dirCellGroup= repmat('XX',length(allDirs), 1);

timeStruct.Onset.Start= 740e-3;
timeStruct.Onset.FFR_HardEnd= 20e-3 + timeStruct.Onset.Start;

ffr_demo.hi_ffr_chin= 367;
ffr_demo.nh_ffr_chin= 374;

for dirVar= 1:length(allDirs)
    curDir=  [RootDataDir allDirs(dirVar).name filesep];
    s_files= dir([curDir 'a*_S_*1*']); % choose dirs with attentuation=10 dB
    ChinID=cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp(curDir,'(-Q\d+_)','tokens'), 'UniformOutput', 0));
    
    if ~isempty(s_files)
        if contains(curDir, 'NH')
            dirCellGroup(dirVar, :)= 'NH';
        elseif contains(curDir, 'PTS')
            dirCellGroup(dirVar, :)= 'HI';
        end
        
        s_data_cell= cell(length(s_files), 2);
        nPairs_actual= nan(length(s_files), 1);
        
        for sfile_var=1:length(s_files)
            temp_data= load([curDir s_files(sfile_var).name]);
            temp_data = temp_data.data;
            s_data_cell{sfile_var, 1}= temp_data.AD_Data.AD_Avg_PO_V{1};
            s_data_cell{sfile_var, 2}= temp_data.AD_Data.AD_Avg_NP_V{1};
            
            nPairs_actual(sfile_var)= temp_data.Stimuli.RunLevels_params.nPairs_actual;
        end
        
        s_atten=temp_data.Stimuli.atten_dB;
        
        if s_atten<15
            
            s_data_pos= zeros(1, length(s_data_cell{sfile_var,1}));
            s_data_neg= zeros(1, length(s_data_cell{sfile_var,2}));
            fs_data= temp_data.Stimuli.RPsamprate_Hz;
            
            for i=1:length(s_files)
                s_data_pos= s_data_pos + s_data_cell{i, 1}*nPairs_actual(i)/sum(nPairs_actual);
                s_data_neg= s_data_neg + s_data_cell{i, 2}*nPairs_actual(i)/sum(nPairs_actual);
            end
            
            ffr_uV_scaleFactor= 1e6/temp_data.AD_Data.Gain;
            s_data_pos= ffr_uV_scaleFactor*gen_resample(s_data_pos, fs_data, fs);
            s_data_neg= ffr_uV_scaleFactor*gen_resample(s_data_neg, fs_data, fs);
            
            curFilt= helper.get_filter(fs);
            s_data_pos_filt= filtfilt(curFilt, s_data_pos);
            s_data_neg_filt= filtfilt(curFilt, s_data_neg);
            s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
            s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
            t_data= (1:length(s_data_env))/fs;
            
            
            if ChinID <= 370
                FFR_latency_correction= 0e-3;
            else
                FFR_latency_correction= 4.69e-3;
            end
            
            OnsetInds_Hard= (t_data> (timeStruct.Onset.Start+FFR_latency_correction)) & (t_data<(timeStruct.Onset.FFR_HardEnd+FFR_latency_correction)); % 10 ms in Delgutte's paper
            
            ffr_OnsetMask= (OnsetInds_Hard==1);
            ffr_SustainedMask= (t_data > (timeStruct.Sustained.Start+FFR_latency_correction)) & (t_data < (timeStruct.Sustained.End+FFR_latency_correction)) ; % 10 ms in Delgutte's paper
            
            FFR_onset(dirVar)= range(s_data_env(ffr_OnsetMask)) ;
            FFR_onset_nf(dirVar)= range(s_data_tfs(ffr_SustainedMask));
            
            figure(1);
            sp_ax(4)= subplot(3, 2, 4);
            if ChinID==ffr_demo.nh_ffr_chin
                yyaxis left
                hold on;
                Astim= 4; % uV
                plot(t_data, Astim+s_data_env, '-', 'Color', get_color('b'), 'linew', plt.lw1);
                plot(sig.time, -.75*Astim + sig.data / max(abs(sig.data)) * Astim / 3, '-', 'Color', get_color('gray'));
                set(gca, 'YColor', get_color('k'));
                ylabel('FFR_{ENV} Amp (\muV)');
                
                
                yyaxis right;
                hold on;
                plot(t_data, ffr_OnsetMask, '-', 'color', get_color('c'), 'linew', plt.lw2)
                ylim([-.5 5]);
                xlim([0.55 0.9]);
                set(gca, 'YColor', get_color('k'), 'YTick', [0 1]);
                ylabel('Onset mask');
            elseif ChinID==ffr_demo.hi_ffr_chin
                yyaxis left
                hold on;
                plot(t_data, s_data_env, '-', 'Color', get_color('r'), 'linew', plt.lw1);
                txtHan(5)= text(.05, .95, 'E', 'Units', 'normalized', 'FontWeight', 'bold');
                xlabel('Time (s)');
            end
        end
    end
end

yTickScale= 1.0;
ytickValsDefault= [-yTickScale 0 yTickScale];
ytickVals= [ytickValsDefault Astim+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault, 1, 2)), 'UniformOutput', false);

axes(sp_ax(4));
yyaxis left
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs);

figure(1);
sp_ax(6)= subplot(3, 2, 6);
hold on;
boxplot(FFR_onset, dirCellGroup, 'GroupOrder', {'NH', 'HI'});
plot([1 2], [5.2 5.2], 'k-', 'linew', 1);
plot(1.5, 5.65, 'kp', 'markersize', 6, 'linew', 1);
ylim([.5 5.9])

ylabel('FFR Onset Amp (\muV)');
txtHan(6)= text(.05, .95, 'F', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

nhChinsInds= strcmp(cellstr(dirCellGroup), 'NH');
[~, pFFR]= ttest2(FFR_onset(nhChinsInds), FFR_onset(~nhChinsInds));


set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len);
set(txtHan, 'FontSize', 11);



%% define new axes for AB
Xcorner= .07;
Xwidth= .37;
Xshift= .14;
Ycorner= .062;
Ywidth= .275;
Yshift= .05;

% D
set(sp_ax(5),'Position', [Xcorner Ycorner Xwidth Ywidth])
drawnow

% C
set(sp_ax(3),'Position', [Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% A
set(sp_ax(1),'Position', [Xcorner Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

% F
set(sp_ax(6),'Position', [Xcorner+Xwidth+Xshift Ycorner Xwidth .75*Ywidth])
drawnow

% E
set(sp_ax(4),'Position', [Xcorner+Xwidth+Xshift Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% B
set(sp_ax(2),'Position', [Xcorner+Xwidth+Xshift Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow



if saveFig
    print([dirStruct.png 'Fig6_fric_s'], '-dpng',  '-r600');
end

%% table for FFR 

var_ffr_chinID= cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp({allDirs.name}','(-Q\d+_)','tokens'), 'UniformOutput', 0);
var_ffr_hearingStatus= nhDirInd(:);
var_ffr_OnsetRate= FFR_onset(:);

% without NaN version 
Fig6_tbl_s_coding_FFR= table(var_ffr_chinID(:),   var_ffr_hearingStatus(:),   var_ffr_OnsetRate(:), ...
    'VariableNames', {          'ChinID',              'HearingStatus',          'OnRate'});
if saveStats
    writetable(Fig6_tbl_s_coding_FFR, [dirStruct.stats 'Fig6_FFR_s_all.txt'], 'Delimiter', 'tab');
end