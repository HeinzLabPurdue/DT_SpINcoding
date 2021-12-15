clear;
clc;

saveFig= 0;
saveStats= 0;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];

subtractSR= 1;

figSize_cm= [15 5 17 18];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);


%% All other NH chins
anl.NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
anl.HIchins= [358 360 361 362 367 370];
all_chinID= [anl.NHchins, anl.HIchins];

demo.NH_chin_track_unit= [379, 2, 6]; % One more NH example = [322, 3, 2] [338, 3, 2]
demo.HI_chin_track_unit= [370, 4, 1]; % good HI example 
all_chinDanishData= helper.load_danish_chin_data(anl);
nhChinInd= find(ismember([ [all_chinDanishData.chinID]', [all_chinDanishData.track]', [all_chinDanishData.unit]' ], demo.NH_chin_track_unit, 'rows'));
hiChinInd= find(ismember([ [all_chinDanishData.chinID]', [all_chinDanishData.track]', [all_chinDanishData.unit]' ], demo.HI_chin_track_unit, 'rows'));


%% Examples showing loss of onset
fs= 5e3;

% Analysis parameters for /d/
anl_d.tStart= 8e-3;
anl_d.tEnd= 32e-3;
anl_d.tZoom= 10e-3;
anl_d.winHardOnset= 8e-3;
anl_d.winSoftOnset= 0e-3;
anl_d.winSustained= 15e-3;

% Analysis parameters for /d/
anl_g.tStart= 197e-3;
anl_g.tEnd= 221e-3;
anl_g.tZoom= 10e-3;
anl_g.winHardOnset= 8e-3;
anl_g.winSoftOnset= 0e-3;
anl_g.winSustained= 15e-3;

% Analysis parameters for both
anl.tBinEdges= (anl_d.tStart-anl_d.tZoom):(1/fs):(anl_g.tEnd+anl_d.tZoom);
anl.tBinCenters= ( anl.tBinEdges(1:end-1) + anl.tBinEdges(2:end) )/2;
anl.nw= 5;

% Time structure for /d/
timeStruct.Onset.Start.d= anl_d.tStart;
timeStruct.Onset.HardEnd.d= timeStruct.Onset.Start.d + anl_d.winHardOnset; % Using values from Delgutte and Kiang 1984 (Fricative coding)
timeStruct.Onset.SoftEnd.d= timeStruct.Onset.HardEnd.d + anl_d.winSoftOnset;

timeStruct.Sustained.Start.d= anl_d.tEnd-anl_d.winSustained;
timeStruct.Sustained.End.d= anl_d.tEnd;

OnsetInds_d.Hard= (anl.tBinCenters> timeStruct.Onset.Start.d) & (anl.tBinCenters<(timeStruct.Onset.HardEnd.d)); % 10 ms in Delgutte's paper
OnsetInds_d.Soft= 0*OnsetInds_d.Hard;
temp_soft_start= find(anl.tBinCenters>=(timeStruct.Onset.HardEnd.d), 1);
temp_soft_end= find(anl.tBinCenters<timeStruct.Onset.SoftEnd.d, 1, 'last');
OnsetInds_d.Soft(temp_soft_start:temp_soft_end)= linspace(1, 0, temp_soft_end - temp_soft_start + 1);

anl_d.OnsetMask= OnsetInds_d.Hard + OnsetInds_d.Soft;
anl_d.SustainedMask= (anl.tBinCenters > (timeStruct.Sustained.Start.d)) & (anl.tBinCenters < (timeStruct.Sustained.End.d)) ; % 10 ms in Delgutte's paper

% Time structure for /g/
timeStruct.Onset.Start.g= anl_g.tStart;
timeStruct.Onset.HardEnd.g= timeStruct.Onset.Start.g + anl_g.winHardOnset; % Using values from Delgutte and Kiang 1984 (Fricative coding)
timeStruct.Onset.SoftEnd.g= timeStruct.Onset.HardEnd.g + anl_g.winSoftOnset;

timeStruct.Sustained.Start.g= anl_g.tEnd-anl_g.winSustained;
timeStruct.Sustained.End.g= anl_g.tEnd;

OnsetInds_g.Hard= (anl.tBinCenters> timeStruct.Onset.Start.g) & (anl.tBinCenters<(timeStruct.Onset.HardEnd.g)); % 10 ms in Delgutte's paper
OnsetInds_g.Soft= 0*OnsetInds_g.Hard;
temp_soft_start= find(anl.tBinCenters>=(timeStruct.Onset.HardEnd.g), 1);
temp_soft_end= find(anl.tBinCenters<timeStruct.Onset.SoftEnd.g, 1, 'last');
OnsetInds_g.Soft(temp_soft_start:temp_soft_end)= linspace(1, 0, temp_soft_end - temp_soft_start + 1);

anl_g.OnsetMask= OnsetInds_g.Hard + OnsetInds_g.Soft;
anl_g.SustainedMask= (anl.tBinCenters > (timeStruct.Sustained.Start.g)) & (anl.tBinCenters < (timeStruct.Sustained.End.g)) ; % 10 ms in Delgutte's paper

% damo data
demo.nh_data= all_chinDanishData(nhChinInd);
demo.hi_data= all_chinDanishData(hiChinInd);
demo.delay= 0;

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

sig.d_inds= (sig.time> anl_d.tStart) & (sig.time<anl_d.tEnd);
sig.d_inds_zoom= (sig.time> (anl_d.tStart-anl_d.tZoom)) & (sig.time<(anl_d.tEnd+anl_d.tZoom));

sig.g_inds= (sig.time> anl_g.tStart) & (sig.time<anl_g.tEnd);
sig.g_inds_zoom= (sig.time> (anl_g.tStart-anl_g.tZoom)) & (sig.time<(anl_g.tEnd+anl_g.tZoom));


%% Plot
plt.lw1= 1;
plt.lw2= 1.5;
plt.lw3= 2.0;
plt.mrkSize= 4;
plt.plotVar= true;
plt.octRange4Avging= .33;
plt.avgType= 'mean';
plt.MINpts= 1;
plt.xtick_freq_val= [.5 1 2 4 8];
plt.xtick_freq_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_freq_val), 'uniformoutput', false);
plt.tick_len= [.02 .02];
plt.ytick_psd= -90:10:-70;
plt.sigAdjust= 2.5e-3; % because AN onset is slightly delayed wrt sig (travelling wave + other delay)
plt.g_time_ticks= [.19:.02:.23];

clf;

sp_ax(1)= subplot(4, 3, 1);
yyaxis right;
hold on;
[~,~,lHan_d]= plot_dpss_psd(sig.data(sig.d_inds), sig.fs, 'nw', anl.nw, 'xunit', 'khz');
set(lHan_d, 'color', get_color('k'), 'linew', plt.lw1);
[~,~,lHan_g]= plot_dpss_psd(sig.data(sig.g_inds), sig.fs, 'nw', anl.nw, 'xunit', 'khz');
set(lHan_g, 'color', get_color('gray'), 'linew', plt.lw1);
set(gca, 'YColor', 'k', 'YTick', plt.ytick_psd);
xlabel('Frequency (kHz)');
ylabel('PSD (dB/Hz)');
ylim([-99 -69]);

yyaxis left;
plot(demo.nh_data.TC.freqkHz, demo.nh_data.TC.TCfit, '-', 'Color', 'b', 'linew', plt.lw2);
plot(demo.hi_data.TC.freqkHz, demo.hi_data.TC.TCfit, '-', 'Color', 'r', 'linew', plt.lw2);
set(gca, 'YColor', get_color('k'), 'XTick', plt.xtick_freq_val);
ylim([0 120]);
xlim([.495 10.5]);
Ashift= 0.75*max(demo.nh_uRate_neg)+max(demo.hi_uRate_pos);
ylabel('TC Thresh. (dB SPL)');
txtHan(1)= text(.05, 1.05, 'A', 'Units', 'normalized', 'FontWeight', 'bold');

plosive_line_han(1)= plot(nan, nan, '-', 'color', get_color('k'), 'linew', plt.lw1);
plosive_line_han(2)= plot(nan, nan, '-', 'color', get_color('gray'), 'linew', plt.lw1);
plosive_line_han(3)= plot(nan, nan, '-', 'color', get_color('b'), 'linew', plt.lw2);
plosive_line_han(4)= plot(nan, nan, '-', 'color', get_color('r'), 'linew', plt.lw2);
[legHan_plosive, icons_plosive]= legend(plosive_line_han, '/d/', '/g/', 'NH', 'HI', 'box', 'off', 'Location', 'northeast');
icons_plosive(5).XData= mean(icons_plosive(5).XData) + [0 +.25];
icons_plosive(7).XData= mean(icons_plosive(7).XData) + [0 +.25];
icons_plosive(9).XData= mean(icons_plosive(9).XData) + [0 +.25];
icons_plosive(11).XData= mean(icons_plosive(11).XData) + [0 +.25];
legHan_plosive.Position(1:2)= [.205 .905];

yTickScale= 0.1;
ytickValsDefault= [-yTickScale 0 yTickScale];
ytickVals= [ytickValsDefault Ashift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault, 1, 2)), 'UniformOutput', false);

sp_ax(2)= subplot(4, 3, 2);
hold on;
plot(anl.tBinCenters, Ashift+demo.nh_uRate_pos, 'Color', get_color('b'), 'linew', plt.lw1);
plot(anl.tBinCenters, Ashift+-demo.nh_uRate_neg, 'Color', get_color('lb'), 'linew', plt.lw1);
plot(anl.tBinCenters, demo.hi_uRate_pos, 'Color', get_color('r'), 'linew', plt.lw1);
plot(anl.tBinCenters, -demo.hi_uRate_neg, 'Color', get_color('lr'), 'linew', plt.lw1);
plot(sig.time(sig.d_inds_zoom)+plt.sigAdjust, -Ashift+sig.data(sig.d_inds_zoom), 'Color', get_color('k'), 'linew', plt.lw1);
ylab_rate_han= ylabel('Rate (spikes/bin)', 'Units', 'normalized');
ylab_rate_han.Position(1)= -.17;
xlab_time_han= xlabel('Time(sec)', 'units', 'normalized');
xlab_time_han.Position(1)= 1.2;
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs);


yyaxis right;
plot(anl.tBinCenters, anl_d.OnsetMask, '-', 'Color', get_color('g'), 'linew', plt.lw1);
plot(anl.tBinCenters, anl_d.SustainedMask, '-', 'Color', get_color('m'), 'linew', plt.lw1);
ylim([-0.25 5]);
set(gca, 'YTick', [0 1], 'YColor', 'k');
xlim([anl_d.tStart-anl_d.tZoom anl_d.tEnd+anl_d.tZoom]);
box off;
txtHan(2)= text(.05, 1.05, 'B. /d/', 'Units', 'normalized', 'FontWeight', 'bold');
xlim([max(0, anl_d.tStart-anl_d.tZoom), anl_d.tEnd+anl_d.tZoom])

sp_ax(3)= subplot(4, 3, 3);
yyaxis left;
hold on;
plot(anl.tBinCenters, Ashift+demo.nh_uRate_pos, '-', 'Color', get_color('b'), 'linew', plt.lw1);
plot(anl.tBinCenters, Ashift+-demo.nh_uRate_neg, '-', 'Color', get_color('lb'), 'linew', plt.lw1);
plot(anl.tBinCenters, demo.hi_uRate_pos, '-', 'Color', get_color('r'), 'linew', plt.lw1);
plot(anl.tBinCenters, -demo.hi_uRate_neg, '-', 'Color', get_color('lr'), 'linew', plt.lw1);
plot(sig.time(sig.g_inds_zoom)+plt.sigAdjust, -Ashift+sig.data(sig.g_inds_zoom), '-', 'Color', get_color('gray'), 'linew', plt.lw1);
set(gca, 'xtick', plt.g_time_ticks);
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs, 'YColor', 'k');

yyaxis right;
plot(anl.tBinCenters, anl_g.OnsetMask, '-', 'Color', get_color('g'), 'linew', plt.lw1);
plot(anl.tBinCenters, anl_g.SustainedMask, '-', 'Color', get_color('m'), 'linew', plt.lw1);
ylim([-0.25 5]);
set(gca, 'YTick', [0 1], 'YColor', 'k');
box off;
txtHan(3)= text(.05, 1.05, 'C. /g/', 'Units', 'normalized', 'FontWeight', 'bold');
xlim([max(0, anl_g.tStart-anl_g.tZoom), anl_g.tEnd+anl_g.tZoom])
ylabel('Onset & Sustained mask');


%% Onset and sustained rate calculation
onsetRates_d= nan(length(all_chinDanishData), 1);
sustainedRates_d= nan(length(all_chinDanishData), 1);
onsetRates_g= nan(length(all_chinDanishData), 1);
sustainedRates_g= nan(length(all_chinDanishData), 1);

parfor unitVar= 1:length(all_chinDanishData)
    cur_data= all_chinDanishData(unitVar);
    uR_pos= histcounts( cell2mat(cur_data.pos) - cur_data.delay, anl.tBinEdges) / length(cur_data.pos);
    uR_neg= histcounts( cell2mat(cur_data.neg) - cur_data.delay, anl.tBinEdges) / length(cur_data.neg);
    
    onsetRates_d(unitVar)= ( sum(uR_pos.*anl_d.OnsetMask) + sum(uR_neg.*anl_d.OnsetMask) ) / 2 / (anl_d.winHardOnset + .5*anl_d.winSoftOnset);
    sustainedRates_d(unitVar)= ( sum(uR_pos.*anl_d.SustainedMask) + sum(uR_neg.*anl_d.SustainedMask) ) / 2 / anl_d.winSustained;
    
    onsetRates_g(unitVar)= ( sum(uR_pos.*anl_g.OnsetMask) + sum(uR_neg.*anl_g.OnsetMask) ) / 2 / (anl_g.winHardOnset + .5*anl_g.winSoftOnset);
    sustainedRates_g(unitVar)= ( sum(uR_pos.*anl_g.SustainedMask) + sum(uR_neg.*anl_g.SustainedMask) ) / 2 / anl_g.winSustained;
end

all_cf_kHz= [all_chinDanishData.CF_Hz]/1e3;
all_chinIDs= [all_chinDanishData.chinID];
nhInds= ismember([all_chinDanishData.chinID], anl.NHchins)';
hiInds= ismember([all_chinDanishData.chinID], anl.HIchins)';

all_SR= [all_chinDanishData.SR];
all_thresh_dB= [all_chinDanishData.thresh_dB];
all_Q10_local= [all_chinDanishData.Q10local];
all_TTR_dB= [all_chinDanishData.TTR_dB];
all_hearing= cell(length([all_chinDanishData.chinID]), 1);
all_hearing(nhInds)= {'NH'};
all_hearing(hiInds)= {'HI'};

if subtractSR
    % \d\ 
    onsetRates_d= onsetRates_d(:)-all_SR(:);
    sustainedRates_d= sustainedRates_d(:)-all_SR(:);
    
    % \g\
    onsetRates_g= onsetRates_g(:)-all_SR(:);
    sustainedRates_g= sustainedRates_g(:)-all_SR(:);
end


% /d/ Onset
sp_ax(4)= subplot(4, 2, 3);

hold on;
plot(all_cf_kHz(nhInds), onsetRates_d(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), onsetRates_d(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'b');
end

plot(all_cf_kHz(hiInds), onsetRates_d(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), onsetRates_d(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'r');
end

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
ylabel('Onset rate (spikes/s)');
txtHan(4)= text(.05, .95, 'D', 'Units', 'normalized', 'FontWeight', 'bold');
ttlHan(1)= text(.48, 1.05, '/d/', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

% /d/ Sustained
sp_ax(6)= subplot(4, 2, 5);

hold on;
plot(all_cf_kHz(nhInds), sustainedRates_d(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), sustainedRates_d(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'b');
end

plot(all_cf_kHz(hiInds), sustainedRates_d(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), sustainedRates_d(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'r');
end

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
xlab_cf_han= xlabel('Characteristic Frequency (kHz)', 'Units', 'normalized');
xlab_cf_han.Position(1)= 1.15;
xlab_cf_han.Position(2)= -.18;


ylabel('Sustained rate (spikes/s)');
txtHan(6)= text(.05, .95, 'E', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

% /g/ Onset
sp_ax(5)= subplot(4, 2, 4);

hold on;

plot(all_cf_kHz(nhInds), onsetRates_g(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), onsetRates_g(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'b');
end

plot(all_cf_kHz(hiInds), onsetRates_g(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), onsetRates_g(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'r');
end

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
ylabel('Onset rate (spikes/s)');
txtHan(5)= text(.05, .95, 'F', 'Units', 'normalized', 'FontWeight', 'bold');
ttlHan(2)= text(.48, 1.05, '/g/', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

% /g/ Sustained
sp_ax(7)= subplot(4, 2, 6);

hold on;
plot(all_cf_kHz(nhInds), sustainedRates_g(nhInds), 'x', 'Color', get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), sustainedRates_g(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'b');
end

plot(all_cf_kHz(hiInds), sustainedRates_g(hiInds), '+', 'Color', get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), sustainedRates_g(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'r');
end

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
ylabel('Sustained rate (spikes/s)');
txtHan(7)= text(.05, .95, 'G', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

%% add legend
axes(sp_ax(6));
if ~subtractSR
    legLine(1)= plot(nan, nan, '-', 'color', 'b', 'linew', plt.lw3);
    legLine(2)= plot(nan, nan, '--', 'color', 'b', 'linew', plt.lw2);
    legLine(3)= plot(nan, nan, '-', 'color', 'r', 'linew', plt.lw3);
    legLine(4)= plot(nan, nan, '--', 'color', 'r', 'linew', plt.lw2);
    [legHan, icons]= legend(legLine, 'NH (Driven)', 'NH (Spont.)', 'HI (Driven)', 'HI (Spont.)', 'location', 'northwest', 'box', 'off');
    legHan.FontSize= 9;
    icons(5).XData= mean(icons(5).XData) + [-0.1 +.1];
    icons(7).XData= mean(icons(7).XData) + [-0.1 +.1];
    icons(9).XData= mean(icons(9).XData) + [-0.1 +.1];
    icons(11).XData= mean(icons(11).XData) + [-0.1 +.1];
    legHan.Position(1:2)= [.1 .39];
    
else
    legLine(1)= plot(nan, nan, '-', 'color', 'b', 'linew', plt.lw2);
    legLine(2)= plot(nan, nan, '-', 'color', 'r', 'linew', plt.lw2);
    [legHan, icons]= legend(legLine, 'NH', 'HI', 'location', 'northwest', 'box', 'off');
    legHan.FontSize= 9;
    icons(3).XData= mean(icons(3).XData) + [-0.08 +.2];
    icons(5).XData= mean(icons(5).XData) + [-0.08 +.2];
    legHan.Position(1:2)= [.08 .4];
    axis tight;    
end

%%

linkaxes(sp_ax(4:5), 'y');
linkaxes(sp_ax(6:7), 'y');
drawnow;
linkaxes(sp_ax(4:7), 'x');
xlim(sp_ax(4), [.75 10]);

% without NaN version 
Fig7_tbl_plosive_coding_AN= table(all_chinIDs(:), all_cf_kHz(:), all_SR(:), all_hearing(:), all_thresh_dB(:), all_Q10_local(:), all_TTR_dB(:), ...
    sustainedRates_g(:), onsetRates_g(:), sustainedRates_d(:), onsetRates_d, ...
    'VariableNames', {          'ChinID',   'CF_kHz',    'SpontRate',  'HearingStatus',      'thresh_dB',        'Q10local',      'TTR_dB', ...
        'SusRate_g',     'OnRate_g',         'SusRate_d',      'OnRate_d'});
if saveStats
    writetable(Fig7_tbl_plosive_coding_AN, [dirStruct.stats 'Fig7_AN_plosive_all.txt'], 'Delimiter', 'tab');
end

%% Stats 
CFsplit_kHz= 2;
CFrange_low_kHz= [0.5 CFsplit_kHz];
CFrange_high_kHz= [CFsplit_kHz 8];
valid_cf_low= (all_cf_kHz(:)>CFrange_low_kHz(1)) & (all_cf_kHz(:)<CFrange_low_kHz(2));
valid_cf_high= (all_cf_kHz(:)>CFrange_high_kHz(1)) & (all_cf_kHz(:)<CFrange_high_kHz(2));


% /d/ 

% ------ Onset 

% ------ ------ Low Freq
[~, pOnset_dLF]= ttest2(onsetRates_d(nhInds&valid_cf_low), onsetRates_d(hiInds&valid_cf_low));
fprintf('Onset (LowFreq | /d/) p=%.3f: %.1f < CF < %.1f kHz\n', pOnset_dLF, CFrange_low_kHz(1), CFrange_low_kHz(2));

% ------ ------ High Freq
[~, pOnset_dHF]= ttest2(onsetRates_d(nhInds&valid_cf_high), onsetRates_d(hiInds&valid_cf_high));
fprintf('Onset (HighFreq | /d/) p=%.3f: %.1f < CF < %.1f kHz\n', pOnset_dHF, CFrange_high_kHz(1), CFrange_high_kHz(2));


% ------ Sustained

% ------ ------ Low Freq
[~, pSus_dLF]= ttest2(sustainedRates_d(nhInds&valid_cf_low), sustainedRates_d(hiInds&valid_cf_low));
fprintf('Sustained (LowFreq | /d/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_dLF, CFrange_low_kHz(1), CFrange_low_kHz(2));

% ------ ------ High Freq
[~, pSus_dHF]= ttest2(sustainedRates_d(nhInds&valid_cf_high), sustainedRates_d(hiInds&valid_cf_high));
fprintf('Sustained (HighFreq | /d/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_dHF, CFrange_high_kHz(1), CFrange_high_kHz(2));


% /g/ 

% ------ Onset 

% ------ ------ Low Freq
[~, pOnset_gLF]= ttest2(onsetRates_g(nhInds&valid_cf_low), onsetRates_g(hiInds&valid_cf_low));
fprintf('Onset (LowFreq | /g/) p=%.3f: %.1f < CF < %.1f kHz\n', pOnset_gLF, CFrange_low_kHz(1), CFrange_low_kHz(2));

% ------ ------ High Freq
[~, pOnset_gHF]= ttest2(onsetRates_g(nhInds&valid_cf_high), onsetRates_g(hiInds&valid_cf_high));
fprintf('Onset (HighFreq | /g/) p=%.3f: %.1f < CF < %.1f kHz\n', pOnset_gHF, CFrange_high_kHz(1), CFrange_high_kHz(2));


% ------ Sustained

% ------ ------ Low Freq
[~, pSus_gLF]= ttest2(sustainedRates_g(nhInds&valid_cf_low), sustainedRates_g(hiInds&valid_cf_low));
fprintf('Sustained (LowFreq | /g/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_gLF, CFrange_low_kHz(1), CFrange_low_kHz(2));

% ------ ------ High Freq
[~, pSus_gHF]= ttest2(sustainedRates_g(nhInds&valid_cf_high), sustainedRates_g(hiInds&valid_cf_high));
fprintf('Sustained (HighFreq | /g/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_gHF, CFrange_high_kHz(1), CFrange_high_kHz(2));


%% FFR data
RootDataDir= ['FFR_Data' filesep];
allDirs= dir([RootDataDir '*SFR*']);

FFR_onset_d= nan(length(allDirs), 1);
FFR_onset_d_nf= nan(length(allDirs), 1);

FFR_onset_g= nan(length(allDirs), 1);
FFR_onset_g_nf= nan(length(allDirs), 1);

nhDirInd= contains({allDirs.name}', 'NH');
hiDirInd= contains({allDirs.name}', 'PTS');
dirCellGroup= repmat('XX',length(allDirs), 1);

timeStruct.Onset.Start.d= 7e-3;
timeStruct.Onset.FFR_HardEnd.d= 15e-3 + timeStruct.Onset.Start.d;

timeStruct.Onset.Start.g= 194e-3;
timeStruct.Onset.FFR_HardEnd.g= 20e-3 + timeStruct.Onset.Start.g;

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
                FFR_latency_correction= 4.59e-3;
            end
            
            % /d/
            OnsetInds_Hard_d= (t_data> (timeStruct.Onset.Start.d+FFR_latency_correction)) & (t_data<(timeStruct.Onset.FFR_HardEnd.d+FFR_latency_correction)); % 10 ms in Delgutte's paper
            ffr_OnsetMask_d= (OnsetInds_Hard_d==1);
            FFR_onset_d(dirVar)= range(s_data_env(ffr_OnsetMask_d)) ;
            
            % /g/
            OnsetInds_Hard_g= (t_data> (timeStruct.Onset.Start.g+FFR_latency_correction)) & (t_data<(timeStruct.Onset.FFR_HardEnd.g+FFR_latency_correction)); % 10 ms in Delgutte's paper
            
            ffr_OnsetMask_g= (OnsetInds_Hard_g==1);
            
            FFR_onset_g(dirVar)= range(s_data_env(ffr_OnsetMask_g)) ;
            
        end
    end
end

figure(1);
sp_ax(8)= subplot(4, 2, 7);
hold on 
boxplot(FFR_onset_d, dirCellGroup, 'GroupOrder', {'NH', 'HI'});
ylabel('FFR_{ENV} Onset (\muV)');
txtHan(8)= text(.05, .95, 'H', 'Units', 'normalized', 'FontWeight', 'bold');
box off;
plot([1 2], [5 5], 'k-', 'linew', 1);
plot(1.5, 5.5, 'kp', 'markersize', 6, 'linew', 1);
ylim([1.2 5.8]);
xlim([.5 2.5])

nhChinsInds= strcmp(cellstr(dirCellGroup), 'NH');
[~, pFFR_d]= ttest2(FFR_onset_d(nhChinsInds), FFR_onset_d(~nhChinsInds));

sp_ax(9)= subplot(4, 2, 8);
hold on;
boxplot(FFR_onset_g, dirCellGroup, 'GroupOrder', {'NH', 'HI'});
ylabel('FFR_{ENV} Onset (\muV)');
txtHan(9)= text(.05, .95, 'I', 'Units', 'normalized', 'FontWeight', 'bold');
box off;
plot([1 2], [5 5], 'k-', 'linew', 1);
text(1.5, 5.5, 'ns', 'FontSize', 8, 'HorizontalAlignment', 'center');
ylim([1.2 5.8]);
xlim([.5 2.5])

[~, pFFR_g]= ttest2(FFR_onset_g(nhChinsInds), FFR_onset_g(~nhChinsInds));

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len);
set(sp_ax(1:3),'TickLength', 1.5*plt.tick_len);

set(ttlHan, 'FontSize', 13);

%% define new axes for AB
Xcorner= .07;
Xwidth= .39;
Xshift= .1;
Ycorner= .03;
Ywidth= .1925;
Yshift= .05;

Xwidth_top= .24;
Xshift_top= .08;


% H
set(sp_ax(8),'Position', [Xcorner Ycorner Xwidth .75*Ywidth])
drawnow

% E
set(sp_ax(6),'Position', [Xcorner Ycorner+Ywidth+1*Yshift Xwidth Ywidth])
drawnow

% D
set(sp_ax(4),'Position', [Xcorner Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

% I
set(sp_ax(9),'Position', [Xcorner+Xwidth+Xshift Ycorner Xwidth .75*Ywidth])
drawnow

% G
set(sp_ax(7),'Position', [Xcorner+Xwidth+Xshift Ycorner+Ywidth+1*Yshift Xwidth Ywidth])
drawnow

% F
set(sp_ax(5),'Position', [Xcorner+Xwidth+Xshift Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow


% A
set(sp_ax(1),'Position', [Xcorner Ycorner+3*Ywidth+3.5*Yshift Xwidth_top Ywidth])
drawnow

% B
set(sp_ax(2),'Position', [Xcorner+Xwidth_top+1.65*Xshift_top Ycorner+3*Ywidth+3.5*Yshift .95*Xwidth_top Ywidth])
drawnow

% C
set(sp_ax(3),'Position', [Xcorner+1.9*Xwidth_top+2.5*Xshift_top Ycorner+3*Ywidth+3.5*Yshift .95*Xwidth_top Ywidth])
drawnow


if saveFig
    print([dirStruct.png 'Fig7_plosive'], '-dpng',  '-r600');
end


%%
var_ffr_chinID= cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp({allDirs.name}','(-Q\d+_)','tokens'), 'UniformOutput', 0);
var_ffr_hearingStatus= nhDirInd(:);
var_ffr_OnsetRate_d= FFR_onset_d(:);
var_ffr_OnsetRate_g= FFR_onset_g(:);

Fig7_tbl_plosive_coding_FFR= table(var_ffr_chinID(:),   var_ffr_hearingStatus(:),   var_ffr_OnsetRate_d(:),   var_ffr_OnsetRate_g(:), ...
    'VariableNames', {          'ChinID',              'HearingStatus',             'OnRate_d',             'OnRate_g'});
if saveStats
    writetable(Fig7_tbl_plosive_coding_FFR, [dirStruct.stats 'Fig7_FFR_plosive_all.txt'], 'Delimiter', 'tab');
end