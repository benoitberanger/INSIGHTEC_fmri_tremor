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

%%

% timeseries = farm_get_timeseries( data, 'ACC_X','raw',+1,4);
%
% figure(10)
% clf
% hold on
%
% plot((0:size(timeseries,2)-1)/data.fsample, farm.normalize_range( timeseries ))
%
% posture_event = ft_filter_event( data.cfg.event, 'value', 'S 11' );
% volume_event = farm.sequence.get_volume_event( data );
%
% t_posture = [posture_event.sample]-volume_event(1).sample;
% plot(t_posture/data.fsample, 0*ones(size(t_posture)),'sr');


%%

% Fetch ACC data, filter
cleanACC = farm_get_timeseries( data,         'ACC_X',       'raw', +[ 2  100],2 );

% Concat
clean = [cleanACC]; % concatenate in a single variable so the TFA will be perform over all channels

% Downsample for faster convolution
time           = (0:length(clean)-1)/data.fsample;
new_fsample    = 500; % Hz
[ new_timeseries, new_time ] = farm.resample( clean', time, data.fsample, new_fsample/data.fsample ); new_timeseries = new_timeseries';

new_timeseries = farm.normalize_range( new_timeseries' );

dF = 0.1; % Hz
dT = 0.1; % s



W = 2%[1 2 4 8 16]

% Prepare TFA
data_TFA = [];
% data_TFA.trial{1}   = log(new_timeseries+2);
data_TFA.trial{1}   = new_timeseries;
data_TFA.time {1}   = new_time;
data_TFA.fsample    = new_fsample;
data_TFA.label      = {'ACC_X'};
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

%%

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
    %
    %     figure(W)
    %     hold on
    %     plot(selected_power)
    
    %         % normal
    %         reginfo = farm_make_regressor( selected_power, 1/dT, data.sequence.TR );
    %         farm_save_regressor( data, reginfo, ['power_' TFRhann.label{chan}] )
    
            % log
            log_selected_power = selected_power;
            log_selected_power(log_selected_power==0) = min(log_selected_power(log_selected_power~=0));
            log_selected_power = log(log_selected_power);
            log_selected_power = farm.normalize_range( log_selected_power ) + 1;
    %         reginfo = farm_make_regressor( log_selected_power, 1/dT, data.sequence.TR );
    %         farm_save_regressor( data, reginfo, ['logpower_' TFRhann.label{chan}] )
    
    %         % log modulator (unconvolved)
    %         selected_power_TR      = selected_power(1 : round(data.sequence.TR*1/dT) : end);
    %         selected_power_TR_diff = [0 diff(selected_power_TR)]; % first derivative
    %         reginfo_ = reginfo;
    %         reginfo_. reg = selected_power_TR     ;
    %         reginfo_.dreg = selected_power_TR_diff;
    %         farm_save_regressor( data, reginfo_, ['modulator_' TFRhann.label{chan}] )
    
    %         figure(W)
    %         hold on
    %         plot(reginfo_.reg)
    
end

%     figure(W)
%     set(gcf,'name','mean_powspctrm')
%     legend(TFRhann.label)
%     figure(W)
%     set(gcf,'name','selected_power')
%     legend(TFRhann.label)
%     figure(W)
%     set(gcf,'name','selected_power')
%     legend(TFRhann.label)

volume_event = farm.sequence.get_volume_event( data );
first_volume_sammple = volume_event(1).sample;

posture_event_10 = ft_filter_event( data.cfg.event, 'value', 'S 10' );
posture_event_11 = ft_filter_event( data.cfg.event, 'value', 'S 11' );
posture_event = [posture_event_10 posture_event_11];
posture_sample = [posture_event.sample] - first_volume_sammple;
posture_time = posture_sample/data.fsample;
relax_event   = ft_filter_event( data.cfg.event, 'value', 'S 12' );
relax_sample  = [relax_event.sample] - first_volume_sammple;
relax_time = relax_sample/data.fsample;

figure(10)
clf

ax(1) = subplot(2,1,1);
hold on
plot(new_time,new_timeseries)
stem(posture_time,ones(size(posture_time)),'color','red'  ,'marker','none')
stem(  relax_time,ones(size(  relax_time)),'color','black','marker','none')

ax(2) = subplot(2,1,2);
hold on
plot(TFRhann.time,selected_power    ,'color','blue' ,'displayname',   'power@peakfreq')
plot(TFRhann.time,log_selected_power,'color','black','displayname','logpower@peakfreq')
stem(posture_time,ones(size(posture_time)),'displayname','posure','color','red'  ,'marker','none')
stem(  relax_time,ones(size(  relax_time)),'displayname','relax','color','black','marker','none')
legend

linkaxes(ax,'x')
