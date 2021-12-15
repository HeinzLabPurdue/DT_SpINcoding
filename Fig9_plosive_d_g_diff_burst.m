clear;
clc;

saveFig= 0;
saveStats= 0;
dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];

subtractSR= 1;

figSize_cm= [15 5 17.6 15];
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
anl_d.tStart_sec= 8e-3;
anl_d.tEnd_sec= 33e-3;
anl_d.tZoom_sec= 10e-3;

% Analysis parameters for /d/
anl_g.tStart_sec= 197e-3;
anl_g.tEnd_sec= 222e-3;
anl_g.tZoom_sec= 10e-3;

% Analysis parameters for both
anl.tBinEdges_sec= (anl_d.tStart_sec-anl_d.tZoom_sec):(1/fs):(anl_g.tEnd_sec+anl_g.tZoom_sec);
anl.tBinCenters_ms= ( anl.tBinEdges_sec(1:end-1) + anl.tBinEdges_sec(2:end) )/2 * 1e3;
anl.nw= 8;

anl_d.resp_mask= (anl.tBinCenters_ms/1e3> anl_d.tStart_sec) & (anl.tBinCenters_ms/1e3<anl_d.tEnd_sec); 
anl_d.valid_mask= (anl.tBinCenters_ms/1e3> max(0, (anl_d.tStart_sec-anl_d.tZoom_sec))) & (anl.tBinCenters_ms/1e3<(anl_d.tEnd_sec+anl_d.tZoom_sec)); 
anl_g.resp_mask= (anl.tBinCenters_ms/1e3> anl_g.tStart_sec) & (anl.tBinCenters_ms/1e3<anl_g.tEnd_sec); 
anl_g.valid_mask= (anl.tBinCenters_ms/1e3> (anl_g.tStart_sec-anl_g.tZoom_sec)) & (anl.tBinCenters_ms/1e3<(anl_g.tEnd_sec+anl_g.tZoom_sec)); 

% demo data
demo.nh_data= all_chinDanishData(nhChinInd);
demo.hi_data= all_chinDanishData(hiChinInd);
demo.delay= 0;

demo.nh_uRate_pos= histcounts( cell2mat(demo.nh_data.pos) - demo.delay, anl.tBinEdges_sec) / length(demo.nh_data.pos) * fs;
demo.nh_uRate_neg= histcounts( cell2mat(demo.nh_data.neg) - demo.delay, anl.tBinEdges_sec) / length(demo.nh_data.neg) * fs;
demo.hi_uRate_pos= histcounts( cell2mat(demo.hi_data.pos) - demo.delay, anl.tBinEdges_sec) / length(demo.hi_data.pos) * fs;
demo.hi_uRate_neg= histcounts( cell2mat(demo.hi_data.neg) - demo.delay, anl.tBinEdges_sec) / length(demo.hi_data.neg) * fs;

%% Read signal

[sig.data, sig.fs]= audioread(['mat_data' filesep 'FLN_Stim_S_P.wav']);
sig.data= helper.gen_resample(sig.data, sig.fs, 20e3);
sig.data= helper.gen_rescale(sig.data, 70);
sig.fs= 20e3;
sig.time_sec= (1:length(sig.data))/sig.fs;

sig.d_inds= (sig.time_sec> anl_d.tStart_sec) & (sig.time_sec<anl_d.tEnd_sec);
sig.d_inds_zoom= (sig.time_sec> (anl_d.tStart_sec-anl_d.tZoom_sec)) & (sig.time_sec<(anl_d.tEnd_sec+anl_d.tZoom_sec));

sig.g_inds= (sig.time_sec> anl_g.tStart_sec) & (sig.time_sec<anl_g.tEnd_sec);
sig.g_inds_zoom= (sig.time_sec> (anl_g.tStart_sec-anl_g.tZoom_sec)) & (sig.time_sec<(anl_g.tEnd_sec+anl_g.tZoom_sec));

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
plt.g_time_ticks= 190:20:230;

clf;

sp_ax(1)= subplot(4, 3, 1);
yyaxis right;
hold on;
[StimPSD_d_dB,freq_psd_Hz_d,lHan_d]= helper.plot_dpss_psd(sig.data(sig.d_inds), sig.fs, 'nw', anl.nw, 'xunit', 'khz');
set(lHan_d, 'color', helper.get_color('k'), 'linew', plt.lw1);
[StimPSD_g_dB, freq_psd_Hz_g,lHan_g]= helper.plot_dpss_psd(sig.data(sig.g_inds), sig.fs, 'nw', anl.nw, 'xunit', 'khz');
set(lHan_g, 'color', helper.get_color('gray'), 'linew', plt.lw1);
set(gca, 'YColor', 'k', 'YTick', plt.ytick_psd);
xlabel('Frequency (kHz)');
ylabel('PSD (dB/Hz)');
ylim([-99 -69]);

yyaxis left;
plot(demo.nh_data.TC.freqkHz, demo.nh_data.TC.TCfit, '-', 'Color', 'b', 'linew', plt.lw2);
plot(demo.hi_data.TC.freqkHz, demo.hi_data.TC.TCfit, '-', 'Color', 'r', 'linew', plt.lw2);
set(gca, 'YColor', helper.get_color('k'), 'XTick', plt.xtick_freq_val);
ylim([0 120]);
xlim([.495 10.5]);
Ashift= 0.75*max(demo.nh_uRate_neg)+max(demo.hi_uRate_pos);
ylabel('TC thresh. (dB SPL)');
txtHan(1)= text(.05, 1.05, 'A', 'Units', 'normalized', 'FontWeight', 'bold');

plosive_line_han(1)= plot(nan, nan, '-', 'color', helper.get_color('k'), 'linew', plt.lw1);
plosive_line_han(2)= plot(nan, nan, '-', 'color', helper.get_color('gray'), 'linew', plt.lw1);
plosive_line_han(3)= plot(nan, nan, '-', 'color', helper.get_color('b'), 'linew', plt.lw2);
plosive_line_han(4)= plot(nan, nan, '-', 'color', helper.get_color('r'), 'linew', plt.lw2);
[legHan_plosive, icons_plosive]= legend(plosive_line_han, '/d/', '/g/', 'NH', 'HI', 'box', 'off', 'Location', 'northeast');
icons_plosive(5).XData= mean(icons_plosive(5).XData) + [0 +.25];
icons_plosive(7).XData= mean(icons_plosive(7).XData) + [0 +.25];
icons_plosive(9).XData= mean(icons_plosive(9).XData) + [0 +.25];
icons_plosive(11).XData= mean(icons_plosive(11).XData) + [0 +.25];
legHan_plosive.Position(1:2)= [.205 .890];

yTickScale= 0.2 * fs;
ytickValsDefault= [0 yTickScale];
ytickVals= [ytickValsDefault Ashift+ytickValsDefault];
ytickLabs= cellfun(@(x) num2str(x), num2cell(repmat(ytickValsDefault/1e3, 1, 2)), 'UniformOutput', false);

sp_ax(2)= subplot(4, 3, 2);
hold on;
plot(anl.tBinCenters_ms, Ashift+demo.nh_uRate_pos, 'Color', helper.get_color('b'), 'linew', plt.lw1);
plot(anl.tBinCenters_ms, Ashift+-demo.nh_uRate_neg, 'Color', helper.get_color('lb'), 'linew', plt.lw1);
plot(anl.tBinCenters_ms, demo.hi_uRate_pos, 'Color', helper.get_color('r'), 'linew', plt.lw1);
plot(anl.tBinCenters_ms, -demo.hi_uRate_neg, 'Color', helper.get_color('lr'), 'linew', plt.lw1);
plot(1e3*(sig.time_sec(sig.d_inds_zoom)+plt.sigAdjust), -Ashift+sig.data(sig.d_inds_zoom)*fs*2, 'Color', helper.get_color('k'), 'linew', plt.lw1);
ylab_rate_han= ylabel('Rate (spikes/s)', 'Units', 'normalized');
ylab_rate_han.Position(1)= -.12;
xlab_time_han= xlabel('Time (ms)', 'units', 'normalized');
xlab_time_han.Position(1)= 1.2;
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs);
text(-.15, 1.05, 'x10^3', 'Units', 'normalized');
ylim(Ashift*[-1.5 2]);

yyaxis right;
plot(anl.tBinCenters_ms(anl_d.valid_mask), anl_d.resp_mask(anl_d.valid_mask), '-', 'Color', helper.get_color('g'), 'linew', plt.lw1);
ylim([-0.25 5]);
set(gca, 'YTick', [0 1], 'YColor', 'k');
xlim([anl_d.tStart_sec-anl_d.tZoom_sec anl_d.tEnd_sec+anl_d.tZoom_sec]);
box off;
txtHan(2)= text(.05, 1.05, 'B. /d/', 'Units', 'normalized', 'FontWeight', 'bold');
xlim([max(0, anl_d.tStart_sec-anl_d.tZoom_sec), anl_d.tEnd_sec+anl_d.tZoom_sec]*1e3)

sp_ax(3)= subplot(4, 3, 3);
yyaxis left;
hold on;
plot(anl.tBinCenters_ms, Ashift+demo.nh_uRate_pos, '-', 'Color', helper.get_color('b'), 'linew', plt.lw1);
plot(anl.tBinCenters_ms, Ashift+-demo.nh_uRate_neg, '-', 'Color', helper.get_color('lb'), 'linew', plt.lw1);
plot(anl.tBinCenters_ms, demo.hi_uRate_pos, '-', 'Color', helper.get_color('r'), 'linew', plt.lw1);
plot(anl.tBinCenters_ms, -demo.hi_uRate_neg, '-', 'Color', helper.get_color('lr'), 'linew', plt.lw1);
plot(1e3*(sig.time_sec(sig.g_inds_zoom)+plt.sigAdjust), -Ashift+sig.data(sig.g_inds_zoom)*fs*2, '-', 'Color', helper.get_color('gray'), 'linew', plt.lw1);
set(gca, 'xtick', plt.g_time_ticks);
set(gca, 'YTick', ytickVals, 'YTickLabel', ytickLabs, 'YColor', 'k');
text(-.15, 1.05, 'x10^3', 'Units', 'normalized');
ylim(Ashift*[-1.5 2]);

yyaxis right;
plot(anl.tBinCenters_ms(anl_g.valid_mask), anl_g.resp_mask(anl_g.valid_mask), '-', 'Color', helper.get_color('g'), 'linew', plt.lw1);
ylim([-0.25 5]);
set(gca, 'YTick', [0 1], 'YColor', 'k');
box off;
txtHan(3)= text(.05, 1.05, 'C. /g/', 'Units', 'normalized', 'FontWeight', 'bold');
xlim([max(0, anl_g.tStart_sec-anl_g.tZoom_sec), anl_g.tEnd_sec+anl_g.tZoom_sec]*1e3)
ylabel('Burst mask');

%% Onset and sustained rate calculation
BurstRates_d= nan(length(all_chinDanishData), 1);
BurstRates_g= nan(length(all_chinDanishData), 1);

parfor unitVar= 1:length(all_chinDanishData)
    cur_data= all_chinDanishData(unitVar);
    uR_pos= histcounts( cell2mat(cur_data.pos) - cur_data.delay, anl.tBinEdges_sec) / length(cur_data.pos);
    uR_neg= histcounts( cell2mat(cur_data.neg) - cur_data.delay, anl.tBinEdges_sec) / length(cur_data.neg);
    
    BurstRates_d(unitVar)= ( sum(uR_pos.*anl_d.resp_mask) + sum(uR_neg.*anl_d.resp_mask) ) / 2 / (anl_d.tEnd_sec - anl_d.tStart_sec);
    BurstRates_g(unitVar)= ( sum(uR_pos.*anl_g.resp_mask) + sum(uR_neg.*anl_g.resp_mask) ) / 2 / (anl_g.tEnd_sec - anl_g.tStart_sec);
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
    BurstRates_d= BurstRates_d(:)-all_SR(:);
    
    % \g\
    BurstRates_g= BurstRates_g(:)-all_SR(:);
end

% /d/ Rate
sp_ax(4)= subplot(4, 2, 3);

hold on;
plot(all_cf_kHz(nhInds), BurstRates_d(nhInds), 'x', 'Color', helper.get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), BurstRates_d(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'b');
end

plot(all_cf_kHz(hiInds), BurstRates_d(hiInds), '+', 'Color', helper.get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), BurstRates_d(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'r');
end

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
ylabel('Driven rate (spikes/s)');
txtHan(4)= text(.05, 1.05, 'D', 'Units', 'normalized', 'FontWeight', 'bold');
ttlHan(1)= text(.48, 1.1, '/d/', 'Units', 'normalized', 'FontWeight', 'bold');
box off;
xlab_cf_han= xlabel('CF (kHz)', 'Units', 'normalized');
% xlab_cf_han.Position(1)= 1.15;
xlab_cf_han.Position(2)= -.15;

% /g/ Rate
sp_ax(5)= subplot(4, 2, 4);
hold on;

plot(all_cf_kHz(nhInds), BurstRates_g(nhInds), 'x', 'Color', helper.get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), BurstRates_g(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), all_SR(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'b');
end

plot(all_cf_kHz(hiInds), BurstRates_g(hiInds), '+', 'Color', helper.get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), BurstRates_g(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

if ~subtractSR
    [~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), all_SR(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
    set(lHan, 'LineStyle', '--', 'LineWidth', plt.lw2, 'color', 'r');
end

xlab_cf_han= xlabel('CF (kHz)', 'Units', 'normalized');
% xlab_cf_han.Position(1)= 1.15;
xlab_cf_han.Position(2)= -.15;

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
% ylabel('Onset rate (spikes/s)');
txtHan(5)= text(.05, 1.05, 'E', 'Units', 'normalized', 'FontWeight', 'bold');
ttlHan(2)= text(.48, 1.1, '/g/', 'Units', 'normalized', 'FontWeight', 'bold');
box off;

%% add legend
axes(sp_ax(5));
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
    legHan.Position(1:2)= [.86 .63];
    axis tight;    
end


%% comapring /d/ vs /g/ 
Rates_d_minus_g= BurstRates_d - BurstRates_g;

sp_ax(6)= subplot(4, 2, 5);
hold on;

plot(freq_psd_Hz_d, StimPSD_d_dB-StimPSD_g_dB, 'k', 'LineWidth', plt.lw3);
plot([min(all_cf_kHz) max(all_cf_kHz)], zeros(1,2), '--', 'color', helper.get_color('gray'), 'LineWidth', plt.lw3);
set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val);
txtHan(6)= text(.05, 0.95, {'F. /d/ rel. /g/'; '  (stimulus)'}, 'Units', 'normalized', 'FontWeight', 'bold');
xlabel('Frequency (kHz)');
ylabel('Rel. power (dB)');

sp_ax(7)= subplot(4, 2, 6);
hold on;

plot(all_cf_kHz(nhInds), Rates_d_minus_g(nhInds), 'x', 'Color', helper.get_color('b'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(nhInds), Rates_d_minus_g(nhInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'b');

plot(all_cf_kHz(hiInds), Rates_d_minus_g(hiInds), '+', 'Color', helper.get_color('r'), 'MarkerSize', plt.mrkSize);
[~, ~, lHan] = helper.octAVG(all_cf_kHz(hiInds), Rates_d_minus_g(hiInds), plt.MINpts, plt.octRange4Avging, plt.plotVar);
set(lHan, 'LineStyle', '-', 'LineWidth', plt.lw3, 'color', 'r');

plot([min(all_cf_kHz) max(all_cf_kHz)], zeros(1,2), '--', 'color', helper.get_color('gray'), 'LineWidth', plt.lw3);
txtHan(7)= text(.05, 0.95, {'G. /d/ rel. /g/'; '  (response)'}, 'Units', 'normalized', 'FontWeight', 'bold');
ylabel('Rel. rate (spikes/s)');

set(gca, 'XScale', 'log', 'YColor', 'k', 'XTick', plt.xtick_freq_val, 'ytick', -25:25:50);
xlabel('CF (kHz)');
ylim([-40 60]);

%%
linkaxes(sp_ax(4:7), 'x');
xlim(sp_ax(4:7), [.75 10]);
ylim(sp_ax(4:5), [0 145]);

% without NaN version 
Fig7_tbl_plosive_coding_AN= table(all_chinIDs(:), all_cf_kHz(:), all_SR(:), all_hearing(:), all_thresh_dB(:), all_Q10_local(:), all_TTR_dB(:), ...
    BurstRates_g(:), BurstRates_d(:), ...
    'VariableNames', {          'ChinID',   'CF_kHz',    'SpontRate',  'HearingStatus',      'thresh_dB',        'Q10local',      'TTR_dB', ...
        'BurstRates_g',     'BurstRates_d'});
if saveStats
    writetable(Fig7_tbl_plosive_coding_AN, [dirStruct.stats 'Fig7_AN_plosive_burst.txt'], 'Delimiter', 'tab');
end

%% Stats 
CFsplit_kHz= 2;
CFrange_low_kHz= [0.5 CFsplit_kHz];
CFrange_high_kHz= [CFsplit_kHz 8];
valid_cf_low= (all_cf_kHz(:)>CFrange_low_kHz(1)) & (all_cf_kHz(:)<CFrange_low_kHz(2));
valid_cf_high= (all_cf_kHz(:)>CFrange_high_kHz(1)) & (all_cf_kHz(:)<CFrange_high_kHz(2));

% /d/ 


% ------ Burst rate

% ------ ------ Low Freq
[~, pSus_dLF]= ttest2(BurstRates_d(nhInds&valid_cf_low), BurstRates_d(hiInds&valid_cf_low));
fprintf('Sustained (LowFreq | /d/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_dLF, CFrange_low_kHz(1), CFrange_low_kHz(2));

% ------ ------ High Freq
[~, pSus_dHF]= ttest2(BurstRates_d(nhInds&valid_cf_high), BurstRates_d(hiInds&valid_cf_high));
fprintf('Sustained (HighFreq | /d/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_dHF, CFrange_high_kHz(1), CFrange_high_kHz(2));


% /g/ 

% ------ Sustained

% ------ ------ Low Freq
[~, pSus_gLF]= ttest2(BurstRates_g(nhInds&valid_cf_low), BurstRates_g(hiInds&valid_cf_low));
fprintf('Sustained (LowFreq | /g/) p=%.3f: %.1f < CF < %.1f kHz\n', pSus_gLF, CFrange_low_kHz(1), CFrange_low_kHz(2));

% ------ ------ High Freq
[~, pSus_gHF]= ttest2(BurstRates_g(nhInds&valid_cf_high), BurstRates_g(hiInds&valid_cf_high));
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
            s_data_pos= ffr_uV_scaleFactor*helper.gen_resample(s_data_pos, fs_data, fs);
            s_data_neg= ffr_uV_scaleFactor*helper.gen_resample(s_data_neg, fs_data, fs);
            
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
ylabel('SBR_{ENV} onset (\muV)');
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
% ylabel('FFR_{ENV} Onset (\muV)');
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
set(txtHan, 'FontSize', 11);
set(ttlHan, 'FontSize', 14);

%% define new axes for AB
Xcorner= .07;
Xwidth= .39;
Xshift= .1;

Ycorner= .035;
Ywidth3= .21;
Ywidth2= .17;
Ywidth1= .13;
Yshift= .075;

Xwidth_top= .245;
Xshift_top= .075;

% H
set(sp_ax(8),'Position', [Xcorner Ycorner Xwidth Ywidth1])
drawnow

% E
set(sp_ax(6),'Position', [Xcorner Ycorner+Ywidth1+Yshift Xwidth Ywidth2])
drawnow

% D
set(sp_ax(4),'Position', [Xcorner Ycorner+Ywidth1+Ywidth2+2*Yshift Xwidth Ywidth2])
drawnow

% I
set(sp_ax(9),'Position', [Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth1])
drawnow

% G
set(sp_ax(7),'Position', [Xcorner+Xwidth+Xshift Ycorner+Ywidth1+Yshift Xwidth Ywidth2])
drawnow

% F
set(sp_ax(5),'Position', [Xcorner+Xwidth+Xshift Ycorner+Ywidth1+Ywidth2+2*Yshift Xwidth Ywidth2])
drawnow


% A
set(sp_ax(1),'Position', [Xcorner Ycorner+Ywidth1+2*Ywidth2+3.5*Yshift Xwidth_top Ywidth3])
drawnow

% B
set(sp_ax(2),'Position', [Xcorner+Xwidth_top+1.65*Xshift_top Ycorner+Ywidth1+2*Ywidth2+3.5*Yshift .95*Xwidth_top Ywidth3])
drawnow

% C
set(sp_ax(3),'Position', [Xcorner+1.9*Xwidth_top+2.5*Xshift_top Ycorner+Ywidth1+2*Ywidth2+3.5*Yshift .95*Xwidth_top Ywidth3])
drawnow

if saveFig
    print([dirStruct.png 'Fig7_plosive_diff_burst'], '-dpng',  '-r600');
    saveas(gcf, [dirStruct.eps 'Fig7_plosive_diff_burst'], 'epsc');
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
