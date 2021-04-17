function [VPdistance, nSpikes_avg]= analyze_SpkTrn_precision_VP(curData, anl)

curChinID= curData(1).chinID;
curTrack= curData(1).track;
curUnit= curData(1).unit;

fs= 20e3;
[stim, fsOrg]=audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_09_09-Q321_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_N.wav');
stim= gen_resample(stim, fsOrg, fs);

durStim=length(stim)/fs;

if length(curData)>1
    error('asdsad');
end

if ~isempty(curData.pos)
    
    clean_speech_spikes_pos= curData.pos;
    clean_speech_spikes_neg= curData.neg; %#ok<*AGROW>
    
    nSpikes_pos= nanmean(cell2mat(cellfun(@(x,t2) sum(x<=t2), clean_speech_spikes_pos, repmat({durStim}, size(clean_speech_spikes_pos)), 'UniformOutput', false)));
    nSpikes_neg= nanmean(cell2mat(cellfun(@(x,t2) sum(x<=t2), clean_speech_spikes_neg, repmat({durStim}, size(clean_speech_spikes_neg)), 'UniformOutput', false)));
    nSpikes_avg= (nSpikes_pos+nSpikes_neg)/2;
    
    tWindows= anl.tWindows;
    
    VPdistance= cell(numel(tWindows), 1);
    
    for winVar= 1:length(tWindows)
        curWindowSize= tWindows(winVar);
        numSegments= floor(durStim/curWindowSize);
        
        % Initialize 
        curWin_VPdistance= nan(numSegments, numel(anl.VPcosts)); % DC = difcor
        
        for segVar= 1:numSegments
            % Get spike corresponding to the current window
            tStart= (segVar-1)*curWindowSize;
            tEnd= segVar*curWindowSize;
            cur_spikes_pos= cellfun(@(x,t1,t2) x(x>t1 & x<=t2), clean_speech_spikes_pos, repmat({tStart}, size(clean_speech_spikes_pos)), repmat({tEnd}, size(clean_speech_spikes_pos)), ...
                'UniformOutput', false);
            cur_spikes_neg= cellfun(@(x,t1,t2) x(x>t1 & x<=t2), clean_speech_spikes_neg, repmat({tStart}, size(clean_speech_spikes_neg)), repmat({tEnd}, size(clean_speech_spikes_neg)), ...
                'UniformOutput', false);
            
            % having issues with empty cells. Using temporary fix. 
            nSpikes_pos= cellfun(@(x) numel(x), cur_spikes_pos);
            NZpos= sum(nSpikes_pos==0);
            nSpikes_neg= cellfun(@(x) numel(x), cur_spikes_neg);
            NZneg= sum(nSpikes_neg==0);
            
            for costVar=1:length(anl.VPcosts)
                curCost= anl.VPcosts(costVar);
                tempDist_2d_pos= helper.VPspkdlPW(cur_spikes_pos(nSpikes_pos~=0), cur_spikes_pos(nSpikes_pos~=0), curCost);
                
                tempDist_2d_neg= helper.VPspkdlPW(cur_spikes_neg(nSpikes_neg~=0), cur_spikes_neg(nSpikes_neg~=0), curCost);
                
                up_tri_Mask_pos= triu(ones(size(tempDist_2d_pos)), 1)==1;
                dist_pos= [tempDist_2d_pos(up_tri_Mask_pos); repmat(nSpikes_pos(nSpikes_pos~=0), NZpos, 1); zeros(NZpos*(NZpos-1)/2, 1)];

                up_tri_Mask_neg= triu(ones(size(tempDist_2d_neg)), 1)==1;
                dist_neg= [tempDist_2d_neg(up_tri_Mask_neg); repmat(nSpikes_neg(nSpikes_neg~=0), NZneg, 1); zeros(NZneg*(NZneg-1)/2, 1)];
                curWin_VPdistance(segVar, costVar)= (nanmean(dist_pos) + nanmean(dist_neg))/2;
            end
            
        end
        
        VPdistance{winVar}= curWin_VPdistance;
    end
    
else
   VPdistance=  nan;
   nSpikes_avg= nan;
   warning('No match for Q%d/t%d/u%d', curChinID, curTrack, curUnit);
end