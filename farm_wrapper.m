%% Init

clear
clc

assert( ~isempty(which('ft_preprocessing')), 'FieldTrip library not detected. Check your MATLAB paths, or get : https://github.com/fieldtrip/fieldtrip' )
assert( ~isempty(which('farm_rootdir'))    ,      'FARM library not detected. Check your MATLAB paths, or get : https://github.com/benoitberanger/FARM' )


%% Get file & sequence paramters

% data_path = '/network/lustre/iss01/cenir/analyse/irm/users/benoit.beranger/INSIGHTEC_fmri_tremor/nifti/2020_02_07_DEV_382_INISGHTEC_PiloteACC01/electrophy';
% fname     = '2020_02_07_INSIGHTEC_pil_001_01_EM_V1_run1';
% sequence.nVol   = 841;  % integer or NaN, if [] it means use all volumes

% data_path = '/network/lustre/iss01/cenir/analyse/irm/users/benoit.beranger/INSIGHTEC_fmri_tremor/nifti/2020_02_18_DEV_349_INSIGHTEC_PiloteEMGACC_02/electrophy';
% fname     = '2020_02_18_DEV_349_INSIGHTEC_PiloteEMGACC_02_run01';
% sequence.nVol   = 831;  % integer or NaN, if [] it means use all volumes

data_path = fullfile(pwd,'nifti/2020_02_25_ULTRABRAIN_01_002_FJ_INCLUSION/electrophy');
fname     = '2020_02_25_Patient002_FJ_run01';
sequence.nVol   = [];  % integer or NaN, if [] it means use all volumes


fname_eeg = fullfile(data_path, [fname '.eeg' ]);
fname_hdr = fullfile(data_path, [fname '.vhdr']);
fname_mrk = fullfile(data_path, [fname '.vmrk']);

sequence.TR     = 1.000; % in seconds
sequence.nSlice = 72;
sequence.MB     = 6;   % multiband factor

MRI_trigger_message = 'R128';

% In this sample dataset, channels are { 'EXT_D' 'FLE_D' 'EXT_G' 'FLE_G' }
% FARM will be performed on all 4 channels, so I create a regex that will fetch them :
channel_regex = 'FCR|ECR|DEL|BIC|TRI';


%% Load data
% Optimal length for a dataset is a bunch of seconds before the start of
% the fmri sequence, and a bunch of seconds after the end of the fmri
% sequence, before any other sequence.

% Read header & events
cfg           = [];
cfg.dataset   = fname_hdr;
raw_event     = ft_read_event (fname_mrk);
event         = farm_change_marker_value(raw_event, MRI_trigger_message, 'V'); % rename volume marker, just for comfort
event         = farm_delete_marker(event, 'Sync On');                          % not useful for FARM, this marker comes from the clock synchronization device

% Load data
data                    = ft_preprocessing(cfg); % load data
data.cfg.event          = event;                 % store events
data.sequence           = sequence;              % store sequence parameters
data.volume_marker_name = 'V';                   % name of the volume event in data.cfg.event


% Some paramters tuning
data.cfg.intermediate_results_overwrite = false; % don't overwrite files
data.cfg.intermediate_results_save      = true;  % write on disk intermediate results
data.cfg.intermediate_results_load      = true;  % if intermediate result file is detected, to not re-do step and load file

% Plot
% ft_databrowser(data.cfg, data)


%% ------------------------------------------------------------------------
%% FARM
% Main FARM functions are below.

% A lot of functions use what is called "regular expressions" (regex). It allows to recognize patterns in strings of characters
% This a powerfull tool, which is common to almost all programing languages. Open some documentation with : doc regular-expressions


%% Check input data
farm_check_data( data )


%% Channel selection
% In your dataset, you might have different nature of signal, for exemple EMG + Accelerometer.
% To perform FARM pipeline only on EMG, you need to select the corresponding channels.

% Select channel for the next processing steps
data = farm_select_channel( data, channel_regex );

fprintf('channel selected : %s \n', data.selected_channels_name{:})


%% Initial HPF @ 30Hz

data = farm_initial_hpf( data );


%% Which channel with greater artifacts ?

data = farm_detect_channel_with_greater_artifact( data );
fprintf('channel with greater artifacts : %s \n', data.label{data.target_channel})


%% Add slice markers : initialize sdur & dtime

data = farm_add_slice_marker( data );


%% Prepare slice candidates for the template generation

data = farm_pick_slice_for_template( data );


%% Optimize slice markers : optimize sdur & dtime
% with an unconstrained non-linear optimization

data = farm_optimize_sdur_dtime( data );


%% Slice correction : compute slice template using best candidates

data = farm_compute_slice_template( data );


%% Volume correction : replace volume-segment (dtime) by 0
% In the FARM article, this method is more advanced, and overwrite less points
% But I didn't succed to code it properly, so I used a "zero filling"

data = farm_volume_correction( data );


%% Revove noise residuals using PCA
% Here, the templates will be substracted, then PCA will be perform on the residuals.
% PCs will bi fitted to theses residials, and substracted.

data = farm_optimize_slice_template_using_PCA( data );


%% Revove noise residuals using ANC
% ANC will remove the last residuals not fitted by the PCs

% Don't know why ANC diverges in this dataset
% Clue : in Niazy et al., they think the filtering diverges when the amplitude is large,
% which is the case for EMG burst compared to EEG.

% data = farm_adaptive_noise_cancellation( data );


%% Remove slice markers
% More convenient

data = farm_remove_slice_marker( data );


%% Time-frequencey analysis (TFA)

% Fetch EMG data, filter, envelope
cleanEMG    = farm_get_timeseries( data, channel_regex, 'pca_clean', +[30 250]   );
envelopeEMG = farm_emg_envelope( cleanEMG, data.fsample );

% Fetch ACC data, filter
cleanACC = farm_get_timeseries( data,         'ACC',       'raw', +[ 2  15],2 );

% Concat
clean = [envelopeEMG;cleanACC]; % concatenate in a single variable so the TFA will be perform over all channels

% Downsample for faster convolution
time           = (0:length(clean)-1)/data.fsample;
new_fsample    = 500; % Hz
[ new_timeseries, new_time ] = farm.resample( clean', time, data.fsample, new_fsample/data.fsample ); new_timeseries = new_timeseries';

new_timeseries = farm.normalize_range( new_timeseries );

dF = 0.1; % Hz
dT = 0.1; % s

for W = 2%[1 2 4 8 16]
    
    % Prepare TFA
    data_TFA = [];
    % data_TFA.trial{1}   = log(new_timeseries+2);
    data_TFA.trial{1}   = new_timeseries;
    data_TFA.time {1}   = new_time;
    data_TFA.fsample    = new_fsample;
    data_TFA.label      = data.label(2:end);
    data_TFA.sampleinfo = [1 size(data_TFA.trial{1},2)];
    
    cfg_TFA            = [];
    cfg_TFA.method     = 'mtmconvol';               % Select method (choose 'mtmconvol' for regressor)
    cfg_TFA.output     = 'pow';                     % Select output ('pow'=power)
    cfg_TFA.taper      = 'hanning';                 % Windowing ('hanning'=hanning taper)
    cfg_TFA.foi        = 2:dF:10;                  % Frequency range you are interested in (usually 2:0.5:8, make sure you at least include 3-8 Hz)
    cfg_TFA.t_ftimwin  = W*ones(1,length(cfg_TFA.foi));          % Wavelet length (seconds; 1 wavelet per frequency). For practical reasons usually take 2 second (which will contain enough time to detect tremor frequency)
%     cfg_TFA.t_ftimwin  = W./cfg_TFA.foi;
%     cfg_TFA.toi        = '50%';                 % Temporal resolution of you time-frequency representation (resolution in seconds) ('orig': original resolution; 'timedat': one step specified under conf.prepemg.timedat;)
    cfg_TFA.toi        = [data_TFA.time{1}(1) : dT : data_TFA.time{1}(end)]
%     cfg_TFA.toi        = 'all';
    cfg_TFA.pad        = 'maxperlen';               % Padding (use 'maxperlen' for default)
    
    TFRhann = ft_freqanalysis(cfg_TFA, data_TFA)
    
    
    % ft_singleplotTFR
    % figure
    % cfg_spTFR = [];
    % cfg_spTFR.channel = 'ACC_X'
    % ft_singleplotTFR(cfg_spTFR, TFRhann);
    
    
    rangeF = 1; % Hz
    clear h
    
    for chan = 1 : length(data_TFA.label)
        
        powspctrm = squeeze( TFRhann.powspctrm(chan,:,:) );
        powspctrm( isnan( powspctrm ) ) = 0;
        mean_powspctrm = mean( powspctrm , 2);
        mean_powspctrm = farm.normalize_range( mean_powspctrm' );
        
%         figure(W)
%         hold on
%         plot(TFRhann.freq,mean_powspctrm)
        
        [Y,I] = max(mean_powspctrm);
        h(chan) = TFRhann.freq(I);
        
        dN = rangeF / dF / 2;
        selected_power = mean( powspctrm(I-dN:I+dN,:) , 1);
        selected_power = farm.normalize_range( selected_power );
        
%         figure(W)
%         hold on
%         plot(selected_power)
        
        % normal
%         reginfo = farm_make_regressor( selected_power, 1/dT, data.sequence.TR );
%         farm_save_regressor( data, reginfo, ['power_' TFRhann.label{chan}] )
        
        % log
        log_selected_power = selected_power;
        log_selected_power(log_selected_power==0) = min(log_selected_power(log_selected_power~=0));
        log_selected_power = log(log_selected_power);
        log_selected_power = farm.normalize_range( log_selected_power ) + 1;
        reginfo = farm_make_regressor( log_selected_power, 1/dT, data.sequence.TR );
        farm_save_regressor( data, reginfo, ['logpower_' TFRhann.label{chan}] )
        
        % log modulator (unconvolved)
        log_selected_power_TR      = log_selected_power(1 : round(data.sequence.TR*1/dT) : end);
        log_selected_power_TR_diff = [0 diff(log_selected_power_TR)]; % first derivative
        reginfo_ = reginfo;
        reginfo_. reg = log_selected_power_TR     ;
        reginfo_.dreg = log_selected_power_TR_diff;
        farm_save_regressor( data, reginfo_, ['modulator_' TFRhann.label{chan}] )
        
        figure(W)
        hold on
        plot(reginfo.reg)
        
    end
    
%     figure(W)
%     set(gcf,'name','mean_powspctrm')
%     legend(TFRhann.label)
%     figure(W)
%     set(gcf,'name','selected_power')
%     legend(TFRhann.label)
    figure(W)
    set(gcf,'name','log_selected_power')
    legend(TFRhann.label)
    
end


return


%% Plot

% % Raw
% farm_plot_carpet     (data, 'ECR_R', 'raw'      , +[30 250])
% farm_plot_FFT        (data, 'ECR_R', 'raw'      , +[30 250])
% farm_plot_spectrogram(data, 'ECR_R', 'raw'      , +[30 250])
%
% % After processing
% farm_plot_carpet     (data, 'ECR_R', 'pca_clean', +[30 250])
% farm_plot_FFT        (data, 'ECR_R', 'pca_clean', +[30 250])
% farm_plot_spectrogram(data, 'ECR_R', 'pca_clean', +[30 250])


%% Convert clean EMG to regrssors

% farm_plot_FFT( data, {'ECR','FCR'}, 'pca_clean', +[30 250] )

ECR_R         = farm_get_timeseries( data, 'ECR_R', 'pca_clean', +[30 250] );                  % (1 x nSamples)
ECR_R_reginfo = farm_emg_regressor ( data,  ECR_R );
% farm_plot_regressor(ECR_R_reginfo,'ECR_R')

FCR_R         = farm_get_timeseries( data, 'FCR_R', 'pca_clean', +[30 250] );                  % (1 x nSamples)
FCR_R_reginfo = farm_emg_regressor ( data,  FCR_R );
% farm_plot_regressor(FCR_R_reginfo,'FCR_R')

ECR_L         = farm_get_timeseries( data, 'ECR_L', 'pca_clean', +[30 250] );                  % (1 x nSamples)
ECR_L_reginfo = farm_emg_regressor ( data,  ECR_L );
% farm_plot_regressor(ECR_L_reginfo,'ECR_L')

FCR_L         = farm_get_timeseries( data, 'FCR_L', 'pca_clean', +[30 250] );                  % (1 x nSamples)
FCR_L_reginfo = farm_emg_regressor ( data,  FCR_L );
% farm_plot_regressor(FCR_L_reginfo,'FCR_L')

ALL         = farm_get_timeseries( data, {'ECR','FCR'}, 'pca_clean', +[30 250] );                  % (1 x nSamples)
ALL_reginfo = farm_emg_regressor ( data,  ALL );
% farm_plot_regressor(ALL_reginfo,'ALL')


%% Convert Acceleromter to regressors

% Diagnostic
% farm_plot_FFT( data, 'ACC', 'raw', +[1 15] , 2 )
% farm_plot_FFT( data, 'ACC', 'raw', +[4 6] , 2 )

ACC_X  = farm_get_timeseries( data, 'ACC_X', 'raw', +[1 15] , 2 );
ACC_Y  = farm_get_timeseries( data, 'ACC_Y', 'raw', +[1 15] , 2 );
ACC_Z  = farm_get_timeseries( data, 'ACC_Z', 'raw', +[1 15] , 2 );

ACC_X_reginfo = farm_acc_regressor( data, ACC_X );
ACC_Y_reginfo = farm_acc_regressor( data, ACC_Y );
ACC_Z_reginfo = farm_acc_regressor( data, ACC_Z );

% farm_plot_regressor(ACC_X_reginfo,'ACC_X')
% farm_plot_regressor(ACC_Y_reginfo,'ACC_Y')
% farm_plot_regressor(ACC_Z_reginfo,'ACC_Z')

ACC_XYZ = farm_get_timeseries( data, 'ACC', 'raw', +[1 15] , 2 );
ACC_XYZ_reginfo_mean = farm_acc_regressor( data, ACC_XYZ, 'mean' );
ACC_XYZ_reginfo_pca  = farm_acc_regressor( data, ACC_XYZ, 'pca'  );
% farm_plot_regressor(ACC_XYZ_reginfo_mean,'ACC_XYZ_mean')
% farm_plot_regressor(ACC_XYZ_reginfo_pca ,'ACC_XYZ_pca' )


%% Save EMG regressor

farm_save_regressor( data, ECR_R_reginfo,   'ECR_R' )
farm_save_regressor( data, FCR_R_reginfo,   'FCR_R' )
farm_save_regressor( data, ECR_L_reginfo,   'ECR_L' )
farm_save_regressor( data, FCR_L_reginfo,   'FCR_L' )

farm_save_regressor( data,   ALL_reginfo,   'EMG'   )


%% Save ACC regressor

farm_save_regressor( data, ACC_X_reginfo,   'ACC_X' )
farm_save_regressor( data, ACC_Y_reginfo,   'ACC_Y' )
farm_save_regressor( data, ACC_Z_reginfo,   'ACC_Z' )

farm_save_regressor( data, ACC_XYZ_reginfo_mean, 'ACC_XYZ_mean' )
farm_save_regressor( data, ACC_XYZ_reginfo_pca , 'ACC_XYZ_pca'  )

