% Returns fractional power in dB
function outPower_dB= psd_f0_energy(F0_co_Hz, PSD_dB, freq, fs)

if ~exist('fs', 'var')
    fs= max(freq)*2;
end

if size(PSD_dB, 2)~=1 % PSD_dB should be a column vector
    PSD_dB= PSD_dB';
end

if size(freq, 1)~=1 % make freq a row vector
    freq= freq';
end

filtOrder= 10;
D= designfilt('lowpassiir','FilterOrder', filtOrder, ...
         'PassbandFrequency',F0_co_Hz,'PassbandRipple',0.2, ...
         'SampleRate',fs);

H= abs(freqz(D, freq, fs));
% H_dB= db(abs(H));

PSD_pow= db2pow(PSD_dB); % db2mag(2*PSD_dB);
outPower_dB= .5*db((H*PSD_pow)/sum(PSD_pow));