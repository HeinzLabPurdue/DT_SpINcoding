% Returns fractional power in dB
function [outPower_dB, D]= psd_cf_energy(CF_Hz, PSD_dB, freq, fs)

if ~exist('fs', 'var')
    fs= max(freq)*2;
end

if size(PSD_dB, 2)~=1 % PSD_dB should be a column vector
    PSD_dB= PSD_dB';
end

if size(freq, 1)~=1 % make freq a row vector
    freq= freq';
end
% freq= 1:10e3;
% CF_Hz= .2e3;


CF_kHz= CF_Hz/1e3;
%% FINAL EQUATIONS for Kale Heinz 2010 Q10 data
if CF_kHz<.1 | CF_kHz>10 %#ok<OR2>
    warning('May not be the right Q10');
end

Q10_5th_pctile=(CF_kHz.^0.4019)*10^0.0218;
Q10_95th_pctile=(CF_kHz.^0.4019)*10^0.6405;

Q10_mean= (Q10_5th_pctile+Q10_95th_pctile)/2;
Q10_BW_Hz= CF_Hz/Q10_mean;
Q3_BW_Hz= Q10_BW_Hz* 10^(-7/80)/sqrt(2); % assuming a Fourth Order filter, see Jorgensen Dau 2011
% sqrt(2) is adhoc- but looks like it works great! Need to think why

Fc2= roots([1, -Q3_BW_Hz, -CF_Hz^2]);
Fc2(Fc2<=0)= [];
Fc1= CF_Hz^2/Fc2;
filtOrder= 4;

if Fc2<fs/2
    D = designfilt('bandpassiir', 'FilterOrder', filtOrder, ...
        'HalfPowerFrequency1', Fc1, 'HalfPowerFrequency2', Fc2,...
        'SampleRate', fs);
elseif Fc1<fs/2
    D = designfilt('highpassiir','FilterOrder',filtOrder, ...
        'PassbandFrequency',Fc1,'PassbandRipple',0.2, ...
        'SampleRate',fs);
else
    D= false;
    outPower_dB= nan;
    return;
end

H= abs(freqz(D, freq, fs));
% H_dB= db(abs(H));

PSD_pow= db2mag(2*PSD_dB);
outPower_dB= .5*db((H*PSD_pow)/sum(PSD_pow));

