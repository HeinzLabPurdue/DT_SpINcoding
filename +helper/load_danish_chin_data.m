function chinDanishData= load_danish_chin_data(anl)

ReqFields= {'NHchins', 'HIchins'};
if ~all(ismember(ReqFields, fieldnames(anl)))
    error('Required fields for anl are %s and %s\n', ReqFields{1}, ReqFields{2});
end
anl.ALLchins= [anl.NHchins, anl.HIchins];

chinDanishData= [];

%% Q374 (synchrony capture chin) is processed differently from other chins
if ismember(374, anl.ALLchins)
    ANdataDir= ['ANF_Data' filesep];
    Chin_Q374= 374;
    
    SynCapData= load(sprintf('%sQ%d_allSyncCap.mat', ANdataDir, Chin_Q374));
    SynCapData=SynCapData.SynCapData;
    
    for unitVar= 1:length(SynCapData)
        cur_unit_data= SynCapData(unitVar);
        
        cur_delay= 4.59e-3; % This is the additional delay of using RP2-inv-calib + HB7 (TDT)
        
        chinDanishData(unitVar).chinID= 374;  %#ok<*AGROW>
        chinDanishData(unitVar).track= cur_unit_data.track;  %#ok<*AGROW>
        chinDanishData(unitVar).unit= cur_unit_data.unit;  %#ok<*AGROW>
        chinDanishData(unitVar).CF_Hz= cur_unit_data.CF_Hz;  %#ok<*AGROW>
        chinDanishData(unitVar).thresh_dB= cur_unit_data.thresh_dB;  %#ok<*AGROW>
        chinDanishData(unitVar).TTR_dB= cur_unit_data.TTR_dB;  %#ok<*AGROW>
        chinDanishData(unitVar).Q10local= cur_unit_data.Q10local;  %#ok<*AGROW>
        chinDanishData(unitVar).Q10global= cur_unit_data.Q10global;  %#ok<*AGROW>
        chinDanishData(unitVar).SR= cur_unit_data.SR;
        chinDanishData(unitVar).TC= cur_unit_data.TC;
        chinDanishData(unitVar).pos= cur_unit_data.danish.pos;
        chinDanishData(unitVar).neg= cur_unit_data.danish.neg;
        chinDanishData(unitVar).delay= cur_delay;
        chinDanishData(unitVar).hearing= 'NH';
    end
end

%% Load Q394 
if ismember(394, anl.ALLchins)
    Q394_data= load(['ANF_Data' filesep 'Q394_all_UnitData_65v80.mat']);
    Q394_data= Q394_data.all_UnitData_65v80;
    
    count= numel(chinDanishData);
    for unitVar= 1:length(Q394_data)
        count= count+1;
        cur_unit_data= Q394_data(unitVar);
        
        cur_delay= 4.59e-3; %This is the additional delay of using RP2-inv-calib
        
        chinDanishData(count).chinID= 394;  %#ok<*AGROW>
        chinDanishData(count).track= cur_unit_data.track;  %#ok<*AGROW>
        chinDanishData(count).unit= cur_unit_data.unit;  %#ok<*AGROW>
        chinDanishData(count).CF_Hz= cur_unit_data.CF_Hz;  %#ok<*AGROW>
        chinDanishData(count).thresh_dB= cur_unit_data.thresh_dB;  %#ok<*AGROW>
        chinDanishData(count).TTR_dB= cur_unit_data.TTR_dB;  %#ok<*AGROW>
        chinDanishData(count).Q10local= cur_unit_data.Q10local;  %#ok<*AGROW>
        chinDanishData(count).Q10global= cur_unit_data.Q10global;  %#ok<*AGROW>
        chinDanishData(count).SR= cur_unit_data.SR;
        chinDanishData(count).TC= cur_unit_data.allTCdata;
        chinDanishData(count).pos= cur_unit_data.Spikes_65dB_Pos;
        chinDanishData(count).neg= cur_unit_data.Spikes_65dB_Neg;
        chinDanishData(count).delay= cur_delay;
        chinDanishData(count).hearing= 'NH';
    end
end


%% Load Q395 
if ismember(395, anl.ALLchins)
    Q395_data= load(['ANF_Data' filesep 'Q395_all_UnitData_50v65v80.mat']);
    Q395_data= Q395_data.all_UnitData_50v65v80;
    
    count= numel(chinDanishData);
    for unitVar= 1:length(Q395_data)
        count= count+1;
        cur_unit_data= Q395_data(unitVar);
        
        cur_delay= 4.59e-3; %This is the additional delay of using RP2-inv-calib
        
        chinDanishData(count).chinID= 395;  %#ok<*AGROW>
        chinDanishData(count).track= cur_unit_data.track;  %#ok<*AGROW>
        chinDanishData(count).unit= cur_unit_data.unit;  %#ok<*AGROW>
        chinDanishData(count).CF_Hz= cur_unit_data.CF_Hz;  %#ok<*AGROW>
        chinDanishData(count).thresh_dB= cur_unit_data.thresh_dB;  %#ok<*AGROW>
        chinDanishData(count).TTR_dB= cur_unit_data.TTR_dB;  %#ok<*AGROW>
        chinDanishData(count).Q10local= cur_unit_data.Q10local;  %#ok<*AGROW>
        chinDanishData(count).Q10global= cur_unit_data.Q10global;  %#ok<*AGROW>
        chinDanishData(count).SR= cur_unit_data.SR;
        chinDanishData(count).TC= cur_unit_data.allTCdata;
        chinDanishData(count).pos= cur_unit_data.Spikes_65dB_Pos;
        chinDanishData(count).neg= cur_unit_data.Spikes_65dB_Neg;
        chinDanishData(count).delay= cur_delay;
        chinDanishData(count).hearing= 'NH';
    end
end

%% Load data for all other chins
dirStruct.loading_dir= ['ANF_Data' filesep];
allfiles=dir([dirStruct.loading_dir '*.mat']);

allChinSpikeData = [];
for chinVar=1:length(allfiles)
    curChinID= sscanf(allfiles(chinVar).name, 'Q%d_*');
    if ismember(curChinID, setxor(anl.ALLchins, [374, 394, 395]))
        temp = load([dirStruct.loading_dir allfiles(chinVar).name]);
        allChinSpikeData = [allChinSpikeData; temp.spike_data'];
    end
end

allChinSpikeData= allChinSpikeData(strcmp({allChinSpikeData.noise}', 'SSN')); % Only keep SSN data since SSN and FLN data for the same unit/snr combo is the same
all_chin_track_unit= [ [allChinSpikeData.chinID]', [allChinSpikeData.track]', [allChinSpikeData.unit]'];
[uniq_chin_track_unit, uniqInds]= unique(all_chin_track_unit, 'rows');
uniqSpikeData= allChinSpikeData(uniqInds);
% clear allChinSpikeData;

%% Combine Q374 and other chins
count= length(chinDanishData);

for unitVar= 1:length(uniqSpikeData)
    cur_chin_track_unit= uniq_chin_track_unit(unitVar, :);
    valid_inds= find(ismember(all_chin_track_unit, cur_chin_track_unit, 'rows'));
    
    count= count+1;
    
    cur_unit_data= uniqSpikeData(unitVar);
    if cur_unit_data.chinID>368
        cur_delay= 4.59e-3;
    else
        cur_delay= 0;
    end
    
    chinDanishData(count).pos= [];
    chinDanishData(count).neg= [];
    
    for snrVar= 1:length(valid_inds)
        
        chinDanishData(count).pos= [chinDanishData(count).pos; cur_unit_data.SpikeTrains{1,1}];
        chinDanishData(count).neg= [chinDanishData(count).neg; cur_unit_data.SpikeTrains{1,2}];
    end
    
    
    chinDanishData(count).chinID= cur_unit_data.chinID;
    chinDanishData(count).track= cur_unit_data.track;  %#ok<*AGROW>
    chinDanishData(count).unit= cur_unit_data.unit;  %#ok<*AGROW>
    chinDanishData(count).CF_Hz= cur_unit_data.CF_Hz;  %#ok<*AGROW>
    chinDanishData(count).thresh_dB= cur_unit_data.thresh_dB;  %#ok<*AGROW>
    chinDanishData(count).TTR_dB= cur_unit_data.TTR_dB;  %#ok<*AGROW>
    chinDanishData(count).Q10local= cur_unit_data.Q10local;  %#ok<*AGROW>
    chinDanishData(count).Q10global= cur_unit_data.Q10global;  %#ok<*AGROW>
    chinDanishData(count).SR= cur_unit_data.SR;
    chinDanishData(count).TC= cur_unit_data.TC;
    chinDanishData(count).delay= cur_delay;
    if ismember(cur_chin_track_unit(1), anl.NHchins)
        chinDanishData(count).hearing= 'NH';
    elseif ismember(cur_chin_track_unit(1), anl.HIchins)
        chinDanishData(count).hearing= 'HI';
    end
end