clear;
clc;

saveFig= 0;
saveStats= 0;

dirStruct.png= [pwd filesep 'final_figs' filesep];
dirStruct.stats= [pwd filesep 'tables_for_stats' filesep];
dirStruct.eps= [pwd filesep 'final_figs_eps' filesep];

figSize_cm= [15 5 11.6 5];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

plotDist0_Prec1= 1;
plotX_Rate0_CF1= 0;

%% All other NH chins
anl.NHchins= [321 322 325 338 341 343 346 347 354 355 373 374 379 394 395];
anl.HIchins= [358 360 361 362 367 370];
all_chinDanishData= helper.load_danish_chin_data(anl);

uniq_chin_track_unit= [ [all_chinDanishData.chinID]', [all_chinDanishData.track]', [all_chinDanishData.unit]'];
uniq_chin_track_unit_Cell= num2cell(uniq_chin_track_unit, 2);
unitID= cellfun(@(x) sprintf('Q%d_t%d_u%d', x), uniq_chin_track_unit_Cell, 'UniformOutput', false);

%% Plot variables
nhChinsInds= ismember(uniq_chin_track_unit(:,1), anl.NHchins);
hiChinsInds= ismember(uniq_chin_track_unit(:,1), anl.HIchins);

plt.lw= 1;
plt.markerSize= 4;
plt.fSize= 9;
plt.txt_fSize= 11;
plt.ttl_q_X= .5;
plt.ttl_q_Y= 1.08;
plt.ttl_sp_X= .01;
plt.ttl_sp_Y= plt.ttl_q_Y;

%% Analysis parameters
tWindow= 1300e-3;
wind_cost_combo= {[tWindow, 8], [tWindow, 200], [tWindow, 3000]}; % Stress [envelope, pitch, and TFS]
spLetters= 'ABCD';
txtHan= nan(length(wind_cost_combo), 1);
ttlHan= nan(length(wind_cost_combo), 1);
sp_ax= nan(length(wind_cost_combo), 1);
AllTable= [];

for comboVar=1:length(wind_cost_combo)
    curCombo= wind_cost_combo{comboVar};
    
    anl.tWindows= curCombo(1); 
    anl.VPcosts= curCombo(2); % shifting cost: per second
    
    cell_VPdist_struct= cell(length(all_chinDanishData), 1);
    all_cfs= nan(length(all_chinDanishData), 1);
    all_SRps= nan(length(all_chinDanishData), 1);
    inloop_VPdist= nan(length(all_chinDanishData), 1);
    nSpikes= nan(length(all_chinDanishData), 1);
    
    fName= [dirStruct.stats 'all_VP_data' strrep(num2str(anl.tWindows*1e3), ' ', '_') 'ms'];
    while contains(fName, '__')
        fName= strrep(fName, '__', '_');
    end
    if ~exist([fName '.mat'], 'file')
        parfor unitVar= 1:length(all_chinDanishData)
            curData= all_chinDanishData(unitVar);
            if ~isnan(curData.CF_Hz) & ~isempty(curData.pos) %#ok<AND2>
                [cell_VPdist_struct{unitVar}, nSpikes(unitVar)]= helper.analyze_SpkTrn_precision_VPnorm(curData, anl);
                all_cfs(unitVar)= curData.CF_Hz;
                all_SRps(unitVar)= curData.SR;
                if ~isnan(nSpikes(unitVar))
                    inloop_VPdist(unitVar, :)= nanmean(cell_VPdist_struct{unitVar}{1}, 1);
                end
            end
        end
    else
    end
    
    inloop_Precision= log(1./inloop_VPdist);
    sp_ax(comboVar)= subplot(1, length(wind_cost_combo), comboVar);
    hold on;
    if plotDist0_Prec1
        if ~plotX_Rate0_CF1
            plot(nSpikes(nhChinsInds), inloop_Precision(nhChinsInds), 'x', 'linew', plt.lw, 'markersize', plt.markerSize);
            plot(nSpikes(hiChinsInds), inloop_Precision(hiChinsInds), '+', 'linew', plt.lw, 'markersize', plt.markerSize);
        elseif plotX_Rate0_CF1
            plot(all_cfs(nhChinsInds), inloop_Precision(nhChinsInds), 'x', 'linew', plt.lw, 'markersize', plt.markerSize);
            plot(all_cfs(hiChinsInds), inloop_Precision(hiChinsInds), '+', 'linew', plt.lw, 'markersize', plt.markerSize);
        end
    else
        if ~plotX_Rate0_CF1
            plot(nSpikes(nhChinsInds), inloop_VPdist(nhChinsInds), 'x', 'linew', plt.lw, 'markersize', plt.markerSize);
            plot(nSpikes(hiChinsInds), inloop_VPdist(hiChinsInds), '+', 'linew', plt.lw, 'markersize', plt.markerSize);
            ylim([-0.02 1.02]);
        elseif plotX_Rate0_CF1
            plot(all_cfs(nhChinsInds), inloop_VPdist(nhChinsInds), 'x', 'linew', plt.lw, 'markersize', plt.markerSize);
            plot(all_cfs(hiChinsInds), inloop_VPdist(hiChinsInds), '+', 'linew', plt.lw, 'markersize', plt.markerSize);
            ylim([-0.02 1.02]);
            set(gca, 'XScale', 'log');
        end
    end

    if rem(1/anl.VPcosts*1e3, 1)==0
        txtHan(comboVar)= text(plt.ttl_q_X, plt.ttl_q_Y, sprintf('2/q= %.0f ms', 2/anl.VPcosts*1e3), 'units', 'normalized', 'Interpreter', 'tex', 'HorizontalAlignment', 'center');
    elseif rem(1/anl.VPcosts*1e3*10, 1)==0
         txtHan(comboVar)= text(plt.ttl_q_X, plt.ttl_q_Y, sprintf('2/q= %.1f ms', 2/anl.VPcosts*1e3), 'units', 'normalized', 'Interpreter', 'tex', 'HorizontalAlignment', 'center');
    else 
         txtHan(comboVar)= text(plt.ttl_q_X, plt.ttl_q_Y, sprintf('2/q= %.2f ms', 2/anl.VPcosts*1e3), 'units', 'normalized', 'Interpreter', 'tex', 'HorizontalAlignment', 'center');
    end
    ttlHan(comboVar)= text(plt.ttl_sp_X, plt.ttl_sp_Y, spLetters(comboVar), 'units', 'normalized', 'Interpreter', 'tex', 'HorizontalAlignment', 'center');
    
    non_NAN= ~isnan(inloop_VPdist);

    nPoints= sum(non_NAN);
    tbl = table(anl.tWindows*1e3*ones(nPoints, 1), anl.VPcosts*ones(nPoints,1), unitID(non_NAN)  , nSpikes(non_NAN), all_SRps(non_NAN), ...
        log10(all_cfs(non_NAN)), nhChinsInds(non_NAN), inloop_VPdist(non_NAN), inloop_Precision(non_NAN), ...
        'VariableNames',{'WindowDur_ms', 'VPcost', 'UnitID', 'Rate', 'SRps', 'CF_Hz_log', 'HearingStatus', 'VPdistance', 'Precision'});
    AllTable= [AllTable; tbl]; %#ok<AGROW>
end

if saveStats
    writetable(AllTable, [dirStruct.stats 'Fig2_VPprecision_Table.txt'],'Delimiter', 'tab');
end

if rem(numel(sp_ax), 2)==1
    axes(sp_ax((numel(sp_ax)+1)/2));
    if ~plotX_Rate0_CF1
        xLab_han= xlabel('Discharge rate (spikes/s)', 'units', 'normalized');
    else
        xLab_han= xlabel('CF (kHz)', 'units', 'normalized');
    end
else
    axes(sp_ax(numel(sp_ax)/2));
    if ~plotX_Rate0_CF1
        xLab_han= xlabel('Discharge Rate (spikes/s)', 'units', 'normalized');
    else
        xLab_han= xlabel('CF (kHz)', 'units', 'normalized');
    end
    xLab_han.Position(1)= 1.1;
end

axes(sp_ax(1));
if plotDist0_Prec1
    ylabel('Precision');
else
    ylabel('VPdistance');
end

legend('NH', 'HI', 'box', 'on', 'location', 'southeast', 'box', 'off');
linkaxes(sp_ax, 'x');
set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);
set(txtHan, 'fontsize', plt.txt_fSize);
set(ttlHan, 'fontsize', plt.txt_fSize);
set(findall(gcf,'-property','TickLength'),'TickLength', [.03 .01]);

%% Set axes placement/size
Xcorner=.085;
Xwidth=.25;
Xshift=.07;
Ycorner=.195;
Ywidth=.72;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow

% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow
 
% C
set(sp_ax(3),'Position',[Xcorner+2*Xwidth+2*Xshift Ycorner Xwidth Ywidth])
drawnow

%%
if saveFig

    if plotDist0_Prec1
        postFix1= '_Prec';
    else
        postFix1= '_VPdist_norm'; %#ok<*UNRCH>
    end
    
    if plotX_Rate0_CF1
        postFix2= '_vCF';
    else
        postFix2= '_vRate';
    end
    print([dirStruct.png 'Fig2' postFix1 postFix2], '-dpng',  '-r600');
    saveas(gcf, [dirStruct.eps 'Fig2' postFix1 postFix2], 'epsc');
end