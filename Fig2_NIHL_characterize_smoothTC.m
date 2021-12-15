clear;
clc;

saveFig= 1;
saveStats= 0;

figSize_cm= [10 5 17.6 10];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer', 'DefaultTextFontName', 'DefaultAxesFontName'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters', 'SansSerif', 'SansSerif'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];

if ~isfolder(dirStruct.png)
    mkdir(dirStruct.png);
    mkdir(dirStruct.stats);
end

useBootstrap= 1;
numBootstraps= 10;
bootstrap_frac= 0.75;

%%
NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
HIchins= [358 360 361 362 367 370];
prctPoint= [10 50];

%% ABRs
ABR_DataDir= ['ABR_Data' filesep];
freqs2use_Hz= [.5 1 2 4 8]*1e3;
freqs2use_kHz= freqs2use_Hz(:)/1e3;

abr_thresh_nh= nan(numel(freqs2use_Hz), numel(NHchins));
abr_thresh_preHI= nan(numel(freqs2use_Hz), numel(HIchins));
abr_thresh_postHI= nan(numel(freqs2use_Hz), numel(HIchins));

for nhVar= 1:length(NHchins)
    curChin= NHchins(nhVar);
    curDir= dir(sprintf('%sQ%d*NH*', ABR_DataDir, curChin));
    
    if isempty(curDir)
        fprintf('No baseline for Q%d\n', curChin);
    elseif numel(curDir)>1
        error('SSUp');
    else
        fprintf('Using %s\n', curDir.name);
        curData_fName= dir(sprintf('%s%s%s*.mat', ABR_DataDir, curDir.name, filesep));
        curData= load(sprintf('%s%s%s%s', ABR_DataDir, curDir.name, filesep, curData_fName.name));
        curData= curData.abrs;
        [~, inds]= ismember(freqs2use_Hz, curData.thresholds(:,1));
        abr_thresh_nh(:, nhVar)= curData.thresholds(inds,2);
    end
end

for nhVar= 1:length(HIchins)
    curChin= HIchins(nhVar);
    curDir= dir(sprintf('%sQ%d*NH*', ABR_DataDir, curChin));
    
    if isempty(curDir)
        fprintf('No baseline for Q%d\n', curChin);
    elseif numel(curDir)>1
        error('SSUp');
    else
        fprintf('Using %s\n', curDir.name);
        curData_fName= dir(sprintf('%s%s%s*.mat', ABR_DataDir, curDir.name, filesep));
        curData= load(sprintf('%s%s%s%s', ABR_DataDir, curDir.name, filesep, curData_fName.name));
        curData= curData.abrs;
        [~, inds]= ismember(freqs2use_Hz, curData.thresholds(:,1));
        abr_thresh_preHI(:, nhVar)= curData.thresholds(inds,2);
    end
end

for hiVar= 1:length(HIchins)
    curChin= HIchins(hiVar);
    curDir= dir(sprintf('%sQ%d*HI*', ABR_DataDir, curChin));
    
    if isempty(curDir)
        fprintf('No baseline for Q%d\n', curChin);
    elseif numel(curDir)>1
        error('SSUp');
    else
        fprintf('Using %s\n', curDir.name);
        curData_fName= dir(sprintf('%s%s%s*.mat', ABR_DataDir, curDir.name, filesep));
        curData= load(sprintf('%s%s%s%s', ABR_DataDir, curDir.name, filesep, curData_fName.name));
        curData= curData.abrs;
        [~, inds]= ismember(freqs2use_Hz, curData.thresholds(:,1));
        abr_thresh_postHI(:, hiVar)= curData.thresholds(inds,2);
    end
end

% Plot
plt.fontName='Arial';
plt.mrkSize= 4;
plt.mrkSize2= 4;
plt.lw1= 0.5;
plt.lw2= 1.5;
plt.lw3= 2.5;
plt.mrkLW= 1;
plt.fontSize= 9;
plt.ttlFontSize= 11;
plt.tick_len= [.02 .02];
plt.ax_lw= 1;
plt.Q10_ticks_vals= [.5 1 5 10];
plt.Q10_ticks_labs= cellfun(@(x) num2str(x), num2cell(plt.Q10_ticks_vals), 'UniformOutput', false);
plt.plotVar= true;
plt.step4Avging= 1;
plt.window4Avging= 1;
plt.avgType= 50;
plt.MINpts= 7;
plt.CFshift= 1.05;

sp_ax(1)= subplot(231);
hold on;
plot(freqs2use_kHz, abr_thresh_nh, '-o', 'color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw1);

plot(freqs2use_kHz, abr_thresh_preHI, '-o', 'color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw1);
plot(freqs2use_kHz, nanmean([abr_thresh_nh, abr_thresh_preHI], 2), '-', 'color', 'b', 'markersize', plt.mrkSize, 'linew', plt.lw3);

plot(freqs2use_kHz, abr_thresh_postHI, '-o', 'color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw1);
plot(freqs2use_kHz, nanmean(abr_thresh_postHI, 2), '-', 'color', 'r', 'markersize', plt.mrkSize, 'linew', plt.lw3);

set(gca, 'xscale', 'log', 'xtick', freqs2use_kHz, 'TickLength', plt.tick_len);
xlim([.4 10]);
ylim([-19 90]);

% xlab_han= xlabel('Frequency (kHz)');
ylabHan(1)= ylabel('ABR threshold (dB SPL)', 'Units', 'normalized');
set(gca, 'fontsize', plt.fontSize, 'box', 'off', 'linew', plt.ax_lw, 'FontName', plt.fontName);

abr_thresh_shift= nanmean(abr_thresh_postHI, 2) - nanmean([abr_thresh_nh, abr_thresh_preHI], 2);


% Table
abr_nh_inds= find(all(~isnan(abr_thresh_nh)));
NHchin_str= repmat(cellfun(@(x) ['Q' num2str(x)], num2cell(NHchins(abr_nh_inds)), 'UniformOutput', false), numel(freqs2use_kHz), 1);
NHchins_hearing= repmat({'NH'}, numel(freqs2use_Hz), numel(abr_nh_inds));
freqs_tbl_nh= repmat(freqs2use_kHz(:), 1, numel(abr_nh_inds));

HIchins_str= repmat(cellfun(@(x) ['Q' num2str(x)], num2cell(HIchins), 'UniformOutput', false), numel(freqs2use_kHz), 1);
freqs_tbl_hi= repmat(freqs2use_kHz(:), 1, numel(HIchins));

HIchins_hearing_pre= repmat({'NH'}, numel(freqs2use_Hz), numel(HIchins));
HIchins_hearing_post= repmat({'HI'}, numel(freqs2use_Hz), numel(HIchins));

var_abr.chinID= [NHchin_str, HIchins_str, HIchins_str];
var_abr.hearing= [NHchins_hearing, HIchins_hearing_pre, HIchins_hearing_post];
var_abr.freq_kHz= [freqs_tbl_nh, freqs_tbl_hi, freqs_tbl_hi];
var_abr.abr_thresh= [abr_thresh_nh(:, abr_nh_inds), abr_thresh_preHI, abr_thresh_postHI];


tbl_abr= table(var_abr.chinID(:), var_abr.hearing(:), var_abr.freq_kHz(:), var_abr.abr_thresh(:), ...
    'VariableNames',{'chinID', 'hearing', 'freq_kHz', 'abr_thresh'});

if saveStats
    writetable(tbl_abr, [dirStruct.stats 'Fig1_abr_data_all.txt'], 'Delimiter', 'tab');
end

%% DPOAE
DirStruct.Codes= pwd;
xTicks= freqs2use_Hz;

DPoae_rootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/Baselines/';

% each column is for one animal
freq27= load(['mat_data' filesep 'default_27freq_DPOAE.mat']);
freq27= freq27.default_27freq_DPOAE;

dp_summary_fName=['mat_data' filesep 'DPOAEdata.mat'];

forceRedoDP= false;
if ~exist(dp_summary_fName, 'file') || forceRedoDP
    addpath(DirStruct.Codes);
    dp_data_nh= nan(27, length(NHchins));
    dp_data_preHI= nan(27, length(HIchins));
    dp_data_postHI= nan(27, length(HIchins));
    
    for chinVar=1:length(NHchins)
        chinID=NHchins(chinVar);
        allChinDirs= dir([DPoae_rootDataDir '*' num2str(chinID) '*']);
        
        DPstruct.dirNum= find(contains(lower({allChinDirs.name}'), {'pre', 'nh'}));
        if ~isempty(DPstruct.dirNum)
            DPstruct.DataDir= allChinDirs(DPstruct.dirNum).name;
            DPstruct.dpFile= dir([DPoae_rootDataDir DPstruct.DataDir filesep '*dpoae*']);
            DPstruct.dpFile= [DPoae_rootDataDir DPstruct.DataDir filesep DPstruct.dpFile(1).name];
            DPstruct.calibFile= helper.get_lower_calibFile(DPstruct.dpFile);
            
            % Pre-exposure
            run(DPstruct.calibFile);
            DPstruct.calibData=ans;
            DPstruct.calibData=DPstruct.calibData.CalibData;
            %     run(NH.dpFile);
            %     NH.dpData= ans;
            out_DPOAE_data= helper.my_dpoae_analysis(DPstruct.dpFile);
            DPstruct.dpData=[out_DPOAE_data.dp_amp];
            DPstruct.freqs_Hz= [out_DPOAE_data.freq2];
            calib_at_freqs=0*DPstruct.freqs_Hz;
            
            for freqVar=1:length(DPstruct.freqs_Hz)
                calib_at_freqs(freqVar)= helper.CalibInterp(DPstruct.freqs_Hz(freqVar)/1e3, DPstruct.calibData);
            end
            
            if numel(DPstruct.dpData)~=27
                fprintf('chin %d\n', chinID);
            end
            
            dp_data_nh(:, chinVar)= interp1(DPstruct.freqs_Hz, DPstruct.dpData, freq27);
        end
    end
    
    
    for chinVar=1:length(HIchins)
        chinID=HIchins(chinVar);
        allChinDirs= dir([DPoae_rootDataDir '*' num2str(chinID) '*']);
        
        DPstruct.dirNum= find(contains(lower({allChinDirs.name}'), {'pre', 'nh'}));
        DPstruct.DataDir= allChinDirs(DPstruct.dirNum).name;
        DPstruct.dpFile= dir([DPoae_rootDataDir DPstruct.DataDir filesep '*dpoae*']);
        DPstruct.dpFile= [DPoae_rootDataDir DPstruct.DataDir filesep DPstruct.dpFile(1).name];
        DPstruct.calibFile= helper.get_lower_calibFile(DPstruct.dpFile);
        
        % Pre-exposure
        run(DPstruct.calibFile);
        DPstruct.calibData=ans;
        DPstruct.calibData=DPstruct.calibData.CalibData;
        %     run(NH.dpFile);
        %     NH.dpData= ans;
        out_DPOAE_data= helper.my_dpoae_analysis(DPstruct.dpFile);
        DPstruct.dpData=[out_DPOAE_data.dp_amp];
        DPstruct.freqs_Hz= [out_DPOAE_data.freq2];
        calib_at_freqs=0*DPstruct.freqs_Hz;
        
        for freqVar=1:length(DPstruct.freqs_Hz)
            calib_at_freqs(freqVar)= helper.CalibInterp(DPstruct.freqs_Hz(freqVar)/1e3, DPstruct.calibData);
        end
        
        if numel(DPstruct.dpData)~=27
            fprintf('chin %d\n', chinID);
        end
        
        dp_data_preHI(:, chinVar)= interp1(DPstruct.freqs_Hz, DPstruct.dpData, freq27);
    end
    
    
    for chinVar=1:length(HIchins)
        chinID=HIchins(chinVar);
        allChinDirs= dir([DPoae_rootDataDir '*' num2str(chinID) '*']);
        
        DPstruct.dirNum= find(contains(lower({allChinDirs.name}'), {'post', 'hi', 'pts'}));
        DPstruct.DataDir= allChinDirs(DPstruct.dirNum).name;
        DPstruct.dpFile= dir([DPoae_rootDataDir DPstruct.DataDir filesep '*dpoae*']);
        DPstruct.dpFile= [DPoae_rootDataDir DPstruct.DataDir filesep DPstruct.dpFile(1).name];
        DPstruct.calibFile= helper.get_lower_calibFile(DPstruct.dpFile);
        
        % Pre-exposure
        run(DPstruct.calibFile);
        DPstruct.calibData=ans;
        DPstruct.calibData=DPstruct.calibData.CalibData;
        %     run(NH.dpFile);
        %     NH.dpData= ans;
        out_DPOAE_data= helper.my_dpoae_analysis(DPstruct.dpFile);
        DPstruct.dpData=[out_DPOAE_data.dp_amp];
        DPstruct.freqs_Hz= [out_DPOAE_data.freq2];
        calib_at_freqs=0*DPstruct.freqs_Hz;
        
        for freqVar=1:length(DPstruct.freqs_Hz)
            calib_at_freqs(freqVar)= helper.CalibInterp(DPstruct.freqs_Hz(freqVar)/1e3, DPstruct.calibData);
        end
        
        if numel(DPstruct.dpData)~=27
            fprintf('chin %d\n', chinID);
        end
        
        dp_data_postHI(:, chinVar)= interp1(DPstruct.freqs_Hz, DPstruct.dpData, freq27);
    end
    save(dp_summary_fName, 'dp_data_nh', 'dp_data_preHI', 'dp_data_postHI');
    rmpath(DirStruct.Codes);
else
    load(dp_summary_fName);
end

sp_ax(4)= subplot(234);
hold on;
plot(freq27/1e3, dp_data_nh, '-o', 'color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw1);

plot(freq27/1e3, dp_data_preHI, '-o', 'color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw1);
% plot(freq27/1e3, nanmean(dp_data_nh, 2), '-', 'color', 'b', 'markersize', plt.mrkSize, 'linew', plt.lw3);
[~,~, lHan] = helper.octAVG_upd(repmat(freq27(:)/1e3, size(dp_data_nh, 2)+size(dp_data_preHI, 2), 1), [dp_data_nh(:); dp_data_preHI(:)], ...
    1, plt.step4Avging, plt.plotVar, plt.window4Avging);
set(lHan, 'linestyle', '-', 'color', 'b', 'markersize', plt.mrkSize, 'linew', plt.lw3);


plot(freq27/1e3, dp_data_postHI, '-o', 'color', helper.get_color('r'), 'markersize', plt.mrkSize, 'linew', plt.lw1);
% plot(freq27/1e3, nanmean(dp_data_postHI, 2), '-', 'color', 'r', 'markersize', plt.mrkSize, 'linew', plt.lw3);
[~,~, lHan] = helper.octAVG_upd(repmat(freq27(:)/1e3, size(dp_data_postHI, 2), 1), dp_data_postHI(:), 1, plt.step4Avging, plt.plotVar, plt.window4Avging);
set(lHan, 'linestyle', '-', 'color', 'r', 'markersize', plt.mrkSize, 'linew', plt.lw3);

set(gca, 'xscale', 'log', 'xtick', freqs2use_kHz, 'TickLength', plt.tick_len);
xlim([.4 10]);
xlab_han(1)= xlabel('Frequency (kHz)', 'Units', 'normalized');
ylabHan(2)= ylabel('DPOAE level (dB SPL)', 'Units', 'normalized');
set(gca, 'fontsize', plt.fontSize, 'box', 'off', 'linew', plt.ax_lw, 'FontName', plt.fontName);


% Table
dpoae_nh_inds= find(all(~isnan(dp_data_nh)));
NHchin_str= repmat(cellfun(@(x) ['Q' num2str(x)], num2cell(NHchins(dpoae_nh_inds)), 'UniformOutput', false), numel(freq27), 1);
NHchins_hearing= repmat({'NH'}, numel(freq27), numel(dpoae_nh_inds));
freqs_tbl_nh= repmat(freq27(:), 1, numel(dpoae_nh_inds));

HIchins_str= repmat(cellfun(@(x) ['Q' num2str(x)], num2cell(HIchins), 'UniformOutput', false), numel(freq27), 1);
freqs_tbl_hi= repmat(freq27(:), 1, numel(HIchins));

HIchins_hearing_pre= repmat({'NH'}, numel(freq27), numel(HIchins));
HIchins_hearing_post= repmat({'HI'}, numel(freq27), numel(HIchins));

var_dpoae.chinID= [NHchin_str, HIchins_str, HIchins_str];
var_dpoae.hearing= [NHchins_hearing, HIchins_hearing_pre, HIchins_hearing_post];
var_dpoae.freq_kHz= [freqs_tbl_nh, freqs_tbl_hi, freqs_tbl_hi];
var_dpoae.dpoae_amp= [dp_data_nh(:, dpoae_nh_inds), dp_data_preHI, dp_data_postHI];


tbl_dpoae= table(var_dpoae.chinID(:), var_dpoae.hearing(:), var_dpoae.freq_kHz(:), var_dpoae.dpoae_amp(:), ...
    'VariableNames',{'chinID', 'hearing', 'freq_kHz', 'dpoae_amp'});

if saveStats
    writetable(tbl_dpoae, [dirStruct.stats 'Fig1_dpoae_data_all.txt'], 'Delimiter', 'tab');
end


%% AN threshold data
anl.NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
anl.HIchins= [358 360 361 362 367 370];
uniqSpikeData= helper.load_danish_chin_data(anl);

nhInds= ismember([uniqSpikeData.chinID], NHchins);
hiInds= ismember([uniqSpikeData.chinID], HIchins);
cf_kHz= [uniqSpikeData.CF_Hz]/1e3;
thresh_dBSPL= [uniqSpikeData.thresh_dB];
TTR_dB= [uniqSpikeData.TTR_dB];
Q10_vals= [uniqSpikeData.Q10local];

save(['mat_data' filesep 'chin_thresh_data.mat'], 'nhInds', 'hiInds', 'cf_kHz', 'thresh_dBSPL', 'TTR_dB', 'Q10_vals');

minQ10= .5;
Q10_vals(Q10_vals<minQ10)= minQ10;

sp_ax(2)= subplot(232);
hold on;

plot(cf_kHz(nhInds), thresh_dBSPL(nhInds), 'x', 'color', helper.get_color('b'), 'MarkerSize', plt.mrkSize);
% [meanThreshNH,meanThFreqNH, lHan]= bf_trifilt(cf_kHz(nhInds), thresh_dBSPL(nhInds), plt.plotVar, 'b', plt.step4Avging, plt.avgType);
[meanThreshNH,meanThFreqNH, lHan] = helper.octAVG_upd(cf_kHz(nhInds), thresh_dBSPL(nhInds), plt.MINpts, plt.step4Avging, plt.plotVar, plt.window4Avging, plt.avgType);
set(lHan, 'linew', plt.lw3, 'LineStyle', '-', 'Color', 'b');
[~,~, lHan] = helper.octAVG_upd(cf_kHz(nhInds), thresh_dBSPL(nhInds), plt.MINpts, plt.step4Avging, plt.plotVar, plt.window4Avging, 10); % 10 ptile
set(lHan, 'linew', plt.lw2, 'LineStyle', '--', 'Color', 'b');


plot(cf_kHz(hiInds), thresh_dBSPL(hiInds), '+', 'color', helper.get_color('r'), 'MarkerSize', plt.mrkSize);
% [meanThreshHI,meanThFreqHI, lHan]= bf_trifilt(cf_kHz(hiInds), thresh_dBSPL(hiInds), plt.plotVar, 'r', plt.step4Avging, plt.avgType);
[meanThreshHI,meanThFreqHI, lHan] = helper.octAVG_upd(cf_kHz(hiInds), thresh_dBSPL(hiInds), plt.MINpts, plt.step4Avging, plt.plotVar, plt.window4Avging, plt.avgType);
set(lHan, 'linew', plt.lw3, 'LineStyle', '-', 'Color', 'r');
[~,~, lHan] = helper.octAVG_upd(cf_kHz(hiInds), thresh_dBSPL(hiInds), plt.MINpts, plt.step4Avging, plt.plotVar, plt.window4Avging, 10);
set(lHan, 'linew', plt.lw2, 'LineStyle', '--', 'Color', 'r');

% plot(cf_kHz(outlier_inds), thresh_dBSPL(outlier_inds), '+', 'color', helper.get_color('g'), 'MarkerSize', plt.mrkSize);
set(gca, 'xscale', 'log', 'xtick', freqs2use_kHz, 'TickLength', plt.tick_len);
xlim([.35 10]);
ylim([-19 90]);
ylabHan(3)= ylabel('ANF threshold (dB SPL)', 'Units', 'normalized');


sp_ax(5)= subplot(235);
hold on;
plot(cf_kHz(nhInds), Q10_vals(nhInds), 'x', 'color', helper.get_color('b'), 'MarkerSize', plt.mrkSize);
[meanQ10NH,meanQ10FreqNH, lHan] = helper.octAVG_upd(cf_kHz(nhInds), Q10_vals(nhInds), plt.MINpts, plt.step4Avging, plt.plotVar, plt.window4Avging, plt.avgType);
set(lHan, 'linew', plt.lw3, 'LineStyle', '-', 'Color', 'b');

plot(cf_kHz(hiInds), Q10_vals(hiInds), '+', 'color', helper.get_color('r'), 'MarkerSize', plt.mrkSize);
[meanQ10HI,meanQ10FreqHI, lHan] = helper.octAVG_upd(cf_kHz(hiInds), Q10_vals(hiInds), plt.MINpts, plt.step4Avging, plt.plotVar, plt.window4Avging, plt.avgType);
set(lHan, 'linew', plt.lw3, 'LineStyle', '-', 'Color', 'r');

set(gca, 'xscale', 'log', 'xtick', freqs2use_kHz, 'TickLength', plt.tick_len, 'yscale', 'log', 'YTick', plt.Q10_ticks_vals, 'YTickLabel', plt.Q10_ticks_labs);
xlim([.35 10]);
xlab_han(2)= xlabel('CF (kHz)', 'Units', 'normalized');
ylabHan(4)= ylabel('ANF Q_{10}', 'Units', 'normalized');


legHan_line(1)= plot(nan, nan, '-b', 'linew', plt.lw2);
legHan_line(2)= plot(nan, nan, '-r', 'linew', plt.lw2);
[legHan, icons]= legend(legHan_line, 'NH', 'HI', 'Location', 'northwest', 'box', 'off');
icons(3).XData= mean(icons(4).XData) + [0 +.25];
icons(5).XData= mean(icons(6).XData) + [0 +.25];
legHan.Position(1:2)= [0.41, 0.4];

var_AN.chinID= [uniqSpikeData(:).chinID];
var_AN.hearing= cell(length([uniqSpikeData.chinID]), 1);
var_AN.hearing(nhInds)= {'NH'};
var_AN.hearing(hiInds)= {'HI'};
var_AN.log_cf_Hz= log([uniqSpikeData.CF_Hz]);
var_AN.thresh_dB= [uniqSpikeData.thresh_dB];
var_AN.Q10local= [uniqSpikeData.Q10local];

tbl_AN= table(var_AN.chinID(:), var_AN.hearing(:), var_AN.log_cf_Hz(:), var_AN.thresh_dB(:), var_AN.Q10local(:), TTR_dB(:), ...
    'VariableNames',{'chinID', 'hearing', 'cf_Hz_log', 'an_thresh', 'Q10local', 'TTR_dB'});

writetable(tbl_AN, [dirStruct.stats 'Fig1_AN_data_all.txt'], 'Delimiter', 'tab');


%% Threshld distribution data
plt.ten_pTileMRK= 'v';
plt.ten_pTileCOL= helper.get_color('g');
plt.fifty_pTileMRK= 's';
plt.fifty_pTileCOL= helper.get_color('prp');
plt.abrShiftMRK= '*';

CFcenters= [.5 1 2 4 8];
lowPtilePoint= min(prctPoint);

% Cat data
load('CATdata/Fig4_workspace.mat');

mat_BF_NH_cat= cell2mat(BF_N);
mat_BF_Imild_cat= cell2mat(BF_Imild);
mat_BF_Isev_cat= cell2mat(BF_Isev);
mat_BF_HI_cat= [mat_BF_Imild_cat, mat_BF_Isev_cat];


mat_Thr_NH_cat= cell2mat(Thr_N);
mat_Thr_Imild_cat= cell2mat(Thr_Imild);
mat_Thr_Isev_cat= cell2mat(Thr_Isev);
mat_Thr_HI_cat= [mat_Thr_Imild_cat, mat_Thr_Isev_cat];

ten_prcTileData_NH_cat= nan(numBootstraps, length(CFcenters));
fifty_prcTileData_NH_cat= nan(numBootstraps, length(CFcenters));
ten_prcTileData_HI_cat= nan(numBootstraps, length(CFcenters));
fifty_prcTileData_HI_cat= nan(numBootstraps, length(CFcenters));

stdData_NH_cat= nan(length(CFcenters), 1);
stdData_HI_cat= nan(length(CFcenters), 1);

for cfVar=1:length(CFcenters)
    curCFcent= CFcenters(cfVar);
    curCFstart= curCFcent/sqrt(2);
    curCFend= curCFcent*sqrt(2);
    
    curCF_cat_data_NH= mat_Thr_NH_cat((mat_BF_NH_cat>curCFstart) & (mat_BF_NH_cat<curCFend));
    curCF_cat_data_HI= mat_Thr_HI_cat((mat_BF_HI_cat>curCFstart) & (mat_BF_HI_cat<curCFend));
    
    for bootVar=1:numBootstraps
        nPoints_NH= round(numel(curCF_cat_data_NH)*bootstrap_frac);
        ten_prcTileData_NH_cat(bootVar, cfVar)= prctile(randsample(curCF_cat_data_NH, nPoints_NH, true), min(prctPoint));
        fifty_prcTileData_NH_cat(bootVar, cfVar)= prctile(randsample(curCF_cat_data_NH, nPoints_NH, true), max(prctPoint));

        nPoints_HI= round(numel(curCF_cat_data_HI)*bootstrap_frac);        
        ten_prcTileData_HI_cat(bootVar, cfVar)= prctile(randsample(curCF_cat_data_HI, nPoints_HI, true), min(prctPoint));
        fifty_prcTileData_HI_cat(bootVar, cfVar)= prctile(randsample(curCF_cat_data_HI, nPoints_HI, true), max(prctPoint));
    end
    
%     for prcVar=1:length(prctPoint)
%         curPrctVal= prctPoint(prcVar);
%         prcTileData_NH_cat(cfVar, prcVar)= prctile(, curPrctVal);
%         prcTileData_HI_cat(cfVar, prcVar)= prctile(, curPrctVal);
%     end

    stdData_NH_cat(cfVar)= std(mat_Thr_NH_cat((mat_BF_NH_cat>curCFstart) & (mat_BF_NH_cat<curCFend)));
    stdData_HI_cat(cfVar)= std(mat_Thr_HI_cat((mat_BF_HI_cat>curCFstart) & (mat_BF_HI_cat<curCFend)));
end

% tenPtile_Shift_HI_cat= prcTileData_HI_cat(:, prctPoint==lowPtilePoint) - prcTileData_NH_cat(:, prctPoint==lowPtilePoint);
% fiftyPtile_Shift_HI_cat= prcTileData_HI_cat(:, prctPoint==50) - prcTileData_NH_cat(:, prctPoint==50);
tenPtile_Shift_HI_cat= ten_prcTileData_HI_cat-ten_prcTileData_NH_cat;
fiftyPtile_Shift_HI_cat= fifty_prcTileData_HI_cat-fifty_prcTileData_NH_cat;

%  Chin data
chin_data= load(['mat_data' filesep 'chin_thresh_data.mat']);

mat_BF_NH_chin= chin_data.cf_kHz(chin_data.nhInds);
mat_BF_HI_chin= chin_data.cf_kHz(chin_data.hiInds);

mat_Thr_NH_chin= chin_data.thresh_dBSPL(chin_data.nhInds);
mat_Thr_HI_chin= chin_data.thresh_dBSPL(chin_data.hiInds);

% prcTileData_NH_chin= nan(length(CFcenters), length(prctPoint));
% prcTileData_Imild_chin= nan(length(CFcenters), length(prctPoint));

ten_prcTileData_NH_chin= nan(numBootstraps, length(CFcenters));
fifty_prcTileData_NH_chin= nan(numBootstraps, length(CFcenters));
ten_prcTileData_Imild_chin= nan(numBootstraps, length(CFcenters));
fifty_prcTileData_Imild_chin= nan(numBootstraps, length(CFcenters));

stdData_NH_chin= nan(length(CFcenters), 1);
stdData_HI_chin= nan(length(CFcenters), 1);

for cfVar=1:length(CFcenters)
    curCFcent= CFcenters(cfVar);
    curCFstart= curCFcent/sqrt(2);
    curCFend= curCFcent*sqrt(2);
    
    curCF_chin_data_NH= mat_Thr_NH_chin((mat_BF_NH_chin>curCFstart) & (mat_BF_NH_chin<curCFend));
    curCF_chin_data_HI= mat_Thr_HI_chin( (mat_BF_HI_chin>curCFstart) & (mat_BF_HI_chin<curCFend) );
    
    
    for bootVar=1:numBootstraps
        nPoints_NH= round(numel(curCF_chin_data_NH)*bootstrap_frac);
        ten_prcTileData_NH_chin(bootVar, cfVar)= prctile(randsample(curCF_chin_data_NH, nPoints_NH, true), min(prctPoint));
        fifty_prcTileData_NH_chin(bootVar, cfVar)= prctile(randsample(curCF_chin_data_NH, nPoints_NH, true), max(prctPoint));
        
        nPoints_HI= round(numel(curCF_chin_data_HI)*bootstrap_frac);
        ten_prcTileData_Imild_chin(bootVar, cfVar)= prctile(randsample(curCF_chin_data_HI, nPoints_HI, true), min(prctPoint));
        fifty_prcTileData_Imild_chin(bootVar, cfVar)= prctile(randsample(curCF_chin_data_HI, nPoints_HI, true), max(prctPoint));
    end

%     for prcVar=1:length(prctPoint)
%         curPrctVal= prctPoint(prcVar);
%         
%         if numel(curCF_chin_data_NH)>=plt.MINpts
%             low_prcTileData_NH_chin(cfVar, prcVar)= prctile(curCF_chin_data_NH, curPrctVal);
%         end
%         
%         if numel(curCF_chin_data_HI)>=plt.MINpts
%             low_prcTileData_Imild_chin(cfVar, prcVar)= prctile(curCF_chin_data_HI, curPrctVal);
%         end
%     end
    stdData_NH_chin(cfVar)= std(mat_Thr_NH_chin((mat_BF_NH_chin>curCFstart) & (mat_BF_NH_chin<curCFend)));
    stdData_HI_chin(cfVar)= std(mat_Thr_HI_chin((mat_BF_HI_chin>curCFstart) & (mat_BF_HI_chin<curCFend)));
end

% tenPtile_Shift_chin= low_prcTileData_Imild_chin(:, prctPoint==lowPtilePoint) - low_prcTileData_NH_chin(:, prctPoint==lowPtilePoint);
% fiftyPtile_Shift_chin= low_prcTileData_Imild_chin(:, prctPoint==50) - low_prcTileData_NH_chin(:, prctPoint==50);
tenPtile_Shift_chin= ten_prcTileData_Imild_chin-ten_prcTileData_NH_chin;
fiftyPtile_Shift_chin= fifty_prcTileData_Imild_chin-fifty_prcTileData_NH_chin;

valid_cf_inds_chin= 2:numel(CFcenters);

sp_ax(3)= subplot(233);
hold on;
% plot(CFcenters, fiftyPtile_Shift_chin, plt.fifty_pTileMRK, 'color', plt.fifty_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
% plot(CFcenters, tenPtile_Shift_chin, plt.ten_pTileMRK, 'color', plt.ten_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
errorbar(CFcenters(valid_cf_inds_chin)/plt.CFshift, mean(fiftyPtile_Shift_chin(:, valid_cf_inds_chin), 1), ...
std(fiftyPtile_Shift_chin(:, valid_cf_inds_chin), 0, 1)/sqrt(numBootstraps), plt.fifty_pTileMRK, 'color', plt.fifty_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
errorbar(CFcenters(valid_cf_inds_chin)*plt.CFshift, mean(tenPtile_Shift_chin(:, valid_cf_inds_chin), 1), ...
    std(tenPtile_Shift_chin(:, valid_cf_inds_chin), 0, 1)/sqrt(numBootstraps), plt.ten_pTileMRK, 'color', plt.ten_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
plot(freqs2use_kHz(valid_cf_inds_chin), abr_thresh_shift(valid_cf_inds_chin), '+', 'color', helper.get_color('k'), 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
set(gca, 'xscale', 'log', 'xtick', CFcenters, 'TickLength', plt.tick_len);

leg_line_han(1)= plot(nan, nan, plt.fifty_pTileMRK, 'color', plt.fifty_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
leg_line_han(2)= plot(nan, nan, plt.ten_pTileMRK, 'color', plt.ten_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
leg_line_han(3)= plot(nan, nan, '+', 'color', helper.get_color('k'), 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
[leg3_han, icons3]= legend(leg_line_han, '50%', '10%', 'ABR', 'box', 'off');
count= 1;
icons3(count).Position(1)= icons3(count).Position(1) - .12; count= count+1;
icons3(count).Position(1)= icons3(count).Position(1) - .12; count= count+1;
icons3(count).Position(1)= icons3(count).Position(1) - .12;
leg3_han.Position = [0.76 0.58 0.12 0.11];



ttleHan(1)= title('Chinchilla', 'Units', 'normalized');
ttleHan(1).Position(2)= .925;


sp_ax(6)= subplot(236);
hold on;
% lHan(1)= plot(CFcenters, fiftyPtile_Shift_Isev_cat, plt.fifty_pTileMRK, 'color', helper.get_color('m'), 'MarkerSize', plt.mrkSize2);
% lHan(2)= plot(CFcenters, tenPtile_Shift_Isev_cat, plt.ten_pTileMRK, 'color', helper.get_color('m'), 'MarkerSize', plt.mrkSize2);
% lHan(3)= plot(CFcenters, fiftyPtile_Shift_HI_cat, plt.fifty_pTileMRK, 'color', plt.fifty_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
% lHan(4)= plot(CFcenters, tenPtile_Shift_HI_cat, plt.ten_pTileMRK, 'color', plt.ten_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
lHan(3)= errorbar(CFcenters/plt.CFshift, mean(fiftyPtile_Shift_HI_cat, 1), std(fiftyPtile_Shift_HI_cat, 0, 1)/sqrt(numBootstraps), plt.fifty_pTileMRK, ...
    'color', plt.fifty_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);
lHan(4)= errorbar(CFcenters*plt.CFshift, mean(tenPtile_Shift_HI_cat, 1), std(tenPtile_Shift_HI_cat, 0, 1)/sqrt(numBootstraps), plt.ten_pTileMRK, ...
    'color', plt.ten_pTileCOL, 'MarkerSize', plt.mrkSize2, 'LineWidth', plt.mrkLW);

set(gca, 'xscale', 'log', 'xtick', CFcenters, 'TickLength', plt.tick_len);
xlab_han(3)= xlabel('Band center frequency (kHz)', 'Units', 'normalized');
ylabHan_cat_chin= ylabel('Shift in impaired rel. normal (dB)', 'Units', 'normalized');
ylabHan_cat_chin.Position(1)= -.13;
ylabHan_cat_chin.Position(2)= 1.15;
ttleHan(2)= title('Cat', 'Units', 'normalized');
ttleHan(2).Position(2)= .925;

linkaxes(sp_ax, 'x')
linkaxes(sp_ax([3 6]), 'y')
ylim(sp_ax(3), [-10 max(ylim(sp_ax(3)))]);
ylim(sp_ax(3), [0 53]);

set(findall(gcf,'-property','FontSize'),'FontSize', plt.fontSize);
set(ttleHan,'FontSize', plt.ttlFontSize);
helper.add_subplot_letter_txp(2, 3, 11, .05, .95);

%% define new axes for AB
Xcorner_AB= .057;
Xwidth_AB= .263;
Xshift_AB= .074;
Ycorner_AB= .092;
Ywidth_AB= .40;
Yshift_AB= .085;

% B
set(sp_ax(4),'Position',[Xcorner_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% A
set(sp_ax(1),'Position',[Xcorner_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow

% D
set(sp_ax(5),'Position',[Xcorner_AB+Xwidth_AB+Xshift_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% C
set(sp_ax(2),'Position',[Xcorner_AB+Xwidth_AB+Xshift_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow

% F
set(sp_ax(6),'Position',[Xcorner_AB+2*Xwidth_AB+2*Xshift_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% E
set(sp_ax(3),'Position',[Xcorner_AB+2*Xwidth_AB+2*Xshift_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow


plt.ylab_x12= -.13;
ylabHan(1).Position(1)= plt.ylab_x12;
ylabHan(2).Position(1)= plt.ylab_x12;

plt.ylab_x34= -.135;
ylabHan(3).Position(1)= plt.ylab_x34;
ylabHan(4).Position(1)= plt.ylab_x34;
drawnow;

xlabY= -.11;
xlab_han(1).Position(2)= xlabY;
xlab_han(2).Position(2)= xlabY;
xlab_han(3).Position(2)= xlabY;

%%
if saveFig
    fName= sprintf('Fig1_nihl_model_bstrap');
    print([dirStruct.png fName], '-dpng',  '-r600');
    saveas(gcf, [dirStruct.eps fName], 'epsc');
end
