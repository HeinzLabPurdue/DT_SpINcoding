function curFilt= get_filter(fs_data, FreqWindow)
if ~exist('FreqWindow', 'var')
    HalfPowerFrequency1=50;
    HalfPowerFrequency2=1e3;
else
    HalfPowerFrequency1=FreqWindow(1);
    HalfPowerFrequency2=FreqWindow(2);
end

N_bp_half= 2; 

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end
