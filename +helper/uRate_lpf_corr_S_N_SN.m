% in this function, plus and pos are synonymous. Same with minus and neg
function [covStruct, rmsStruct]=...
    uRate_lpf_corr_S_N_SN(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, ...
    fsOrg, modFreq, N_lp, tStart, tEnd)


debugModeOn1=0;
fs=2e3;

SetOnsetZeroDuration=10e-3;
OnsetInds2SetZero=round(fs*SetOnsetZeroDuration);

corrCoef0Corr1= 1;

%% find p and q for resampling
temp= rats(fs/fsOrg);
temp(temp==' ')= [];
if ~isempty(find(temp=='/', 1))
    pqEndPartition= find(temp=='/');
    p= str2double(temp(1:pqEndPartition-1));
    q= str2double(temp(pqEndPartition+1:end));
else
    p= str2double(temp);
    q= 1;
end

S_rate_plus=resample(S_rate_plus_Org, p, q);
S_rate_minus=resample(S_rate_minus_Org, p, q);
% S_rate_tfs=S_rate_plus; %(S_rate_plus-S_rate_minus)/2;

N_rate_plus=resample(N_rate_plus_Org, p, q);
N_rate_minus=resample(N_rate_minus_Org, p, q);
% N_rate_tfs=N_rate_plus; %(N_rate_plus-N_rate_minus)/2;

SN_rate_plus=resample(SN_rate_plus_Org, p, q);
SN_rate_minus=resample(SN_rate_minus_Org, p, q);

curFilt= designfilt('lowpassiir','FilterOrder',N_lp, ...
    'PassbandFrequency',modFreq,'PassbandRipple',0.2, ...
    'SampleRate',fs);

if debugModeOn1
    fvtool(curFilt);
    xlim([0 2*modFreq]);
end

% for +ve polatiry
S_rate_pos_filt=filtfilt(curFilt, S_rate_plus);
S_rate_pos_filt(1:OnsetInds2SetZero)=0;
SN_rate_pos_filt=filtfilt(curFilt, SN_rate_plus);
SN_rate_pos_filt(1:OnsetInds2SetZero)=0;
N_rate_pos_filt=filtfilt(curFilt, N_rate_plus);
N_rate_pos_filt(1:OnsetInds2SetZero)=0;

% for -ve polatiry
S_rate_neg_filt=filtfilt(curFilt, S_rate_minus);
S_rate_neg_filt(1:OnsetInds2SetZero)=0;
SN_rate_neg_filt=filtfilt(curFilt, SN_rate_minus);
SN_rate_neg_filt(1:OnsetInds2SetZero)=0;
N_rate_neg_filt=filtfilt(curFilt, N_rate_minus);
N_rate_neg_filt(1:OnsetInds2SetZero)=0;

indStart=max(1, round(tStart*fs));
indEnd=min(length(S_rate_plus), round(tEnd*fs));
validINDs=indStart:indEnd;

if corrCoef0Corr1
    covStruct.s_sn.pos= S_rate_pos_filt(validINDs) * SN_rate_pos_filt(validINDs)';
    covStruct.sn_n.pos= SN_rate_pos_filt(validINDs) * N_rate_pos_filt(validINDs)';
    covStruct.s_n.pos= S_rate_pos_filt(validINDs) * N_rate_pos_filt(validINDs)';
    
    covStruct.s_sn.neg= S_rate_neg_filt(validINDs) * SN_rate_neg_filt(validINDs)';
    covStruct.sn_n.neg= SN_rate_neg_filt(validINDs) * N_rate_neg_filt(validINDs)';
    covStruct.s_n.neg= S_rate_neg_filt(validINDs) * N_rate_neg_filt(validINDs)';
    
else
    covStruct.s_sn.pos= corr2(S_rate_pos_filt(validINDs), SN_rate_pos_filt(validINDs));
    covStruct.sn_n.pos=corr2(SN_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
    covStruct.s_n.pos=corr2(S_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
    
    covStruct.s_sn.neg=corr2(S_rate_neg_filt(validINDs), SN_rate_neg_filt(validINDs));
    covStruct.sn_n.neg=corr2(SN_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
    covStruct.s_n.neg=corr2(S_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
end

% Create rms strcture in case normalization is required outside this
% function
rmsStruct.s.pos= rms(S_rate_pos_filt(validINDs));
rmsStruct.sn.pos= rms(SN_rate_pos_filt(validINDs));
rmsStruct.n.pos= rms(N_rate_pos_filt(validINDs));

rmsStruct.s.neg= rms(S_rate_neg_filt(validINDs));
rmsStruct.sn.neg= rms(SN_rate_neg_filt(validINDs));
rmsStruct.n.neg= rms(N_rate_neg_filt(validINDs));
