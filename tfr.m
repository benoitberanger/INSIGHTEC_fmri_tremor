%%

INPUT = ECR_R_reginfo;

% ACC_ZZZ_reginfo.time_in = (0:length(ACC_Z)-1)/data.fsample;
% [ACC_ZZZ_reginfo.in, ACC_ZZZ_reginfo.time_in] = farm.resample( ACC_Z, ACC_ZZZ_reginfo.time_in, data.fsample, 1/10 );
% INPUT = ACC_ZZZ_reginfo;


data_TFR = [];
data_TFR.trial{1}   = INPUT.in;
data_TFR.time {1}   = INPUT.time_in;
data_TFR.fsample    = 500;
data_TFR.label      = {'ECR_L'};
data_TFR.sampleinfo = [1 length(data_TFR.trial{1})];

cfg_TFR            = [];
cfg_TFR.channel    = 'ECR_L';
cfg_TFR.method     = 'mtmconvol';               % Select method (choose 'mtmconvol' for regressor)
cfg_TFR.output     = 'pow';                     % Select output ('pow'=power)  
cfg_TFR.taper      = 'hanning';                 % Windowing ('hanning'=hanning taper)
cfg_TFR.foi        = 2:0.1:10;                   % Frequency range you are interested in (usually 2:0.5:8, make sure you at least include 3-8 Hz)   
cfg_TFR.t_ftimwin  = 10*ones(1,length(cfg_TFR.foi));          % Wavelet length (seconds; 1 wavelet per frequency). For practical reasons usually take 2 second (which will contain enough time to detect tremor frequency)
% cfg_TFR.toi        = '50%';                 % Temporal resolution of you time-frequency representation (resolution in seconds) ('orig': original resolution; 'timedat': one step specified under conf.prepemg.timedat;)
cfg_TFR.toi        = [data_TFR.time{1}(1) : 0.1 : data_TFR.time{1}(end)]
cfg_TFR.pad        = 'maxperlen';               % Padding (use 'maxperlen' for default)

TFRhann = ft_freqanalysis(cfg_TFR, data_TFR)
powspctrm = squeeze(TFRhann.powspctrm);


%% ft_singleplotTFR

cfg_spTFR = [];
ft_singleplotTFR(cfg_spTFR, TFRhann);

%% Surface

figure
surf(TFRhann.time,TFRhann.freq,squeeze(TFRhann.powspctrm(4,:,:)), 'EdgeColor', 'none')
xlabel('time (s)')
ylabel('freq (Hz)')


%% _mean_ plot

figure
mean_powspctrm = nanmean( powspctrm , 2);
plot( TFRhann.freq, mean_powspctrm  )


%% max(_mean_) plot

figure
[Y,I] = max(mean_powspctrm);
d = 0;
selected_power = nanmean( powspctrm(I-d:I+d,:) , 1);
plot( TFRhann.time, selected_power  )











return