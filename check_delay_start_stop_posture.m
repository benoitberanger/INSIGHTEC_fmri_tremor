clear
clc
close all

% data = load('/network/lustre/iss01/cenir/analyse/irm/users/benoit.beranger/INSIGHTEC_fmri_tremor/nifti/2021_04_06_ULTRABRAIN_001_006_CJ_V1_CONTROL/electrophy/2021_04_06_ULTRABRAIN_001_006_CJ_run01_pca_clean.mat');
% acc  = farm_get_timeseries(data,'(^ACC)|(R_ACC)','raw');

data = load('/network/lustre/iss01/cenir/analyse/irm/users/benoit.beranger/INSIGHTEC_fmri_tremor/nifti/2021_03_30_ULTRABRAIN_001_005_BD_V3_CONTROLE/electrophy/2021_03_30_INSIGHTEC_005_BD_V3_CONTROL_run01_pca_clean.mat');
acc  = farm_get_timeseries(data,'(^ACC)|(L_ACC)','raw');

comb = farm_combine_timeseries( acc, 'euclidian_norm' );

timeseries = farm.filter(comb, data.fsample, [1 15], 2);
timeseries = farm.normalize_range(timeseries);
new_fsample = 500;

% plotFFT(timeseries, data.fsample, [0 20])

%% 

time           = (0:size(timeseries,2)-1)/data.fsample;
[ new_timeseries, new_time ] = farm.resample( timeseries, time, data.fsample, new_fsample/data.fsample );

data_emg_acc            = [];
data_emg_acc.trial{1}   = new_timeseries;
data_emg_acc.time {1}   = new_time;
data_emg_acc.fsample    = new_fsample;
data_emg_acc.label      = {'ACC'};
data_emg_acc.sampleinfo = [1 size(data_emg_acc.trial{1},2)];
data_emg_acc.info.channel_idx = 1;

cfg = [];
cfg.minmax_foi = [1 20];
TFA = farm.tfa.perform_time_frequency_analysis( data_emg_acc, [] );
TFA = farm.tfa.postprocessing( TFA, [] );
figH = farm_plot_TFA( data, TFA );


%%

volume_event = farm.sequence.get_volume_event(data);
first_volume_sample = volume_event(1).sample;

posture_start = ft_filter_event( data.cfg.event, 'value', 'S 11' );
posture_stop  = ft_filter_event( data.cfg.event, 'value', 'S 12' );

posture_start_sample = [posture_start.sample] - first_volume_sample;
posture_stop_sample  = [posture_stop .sample] - first_volume_sample;

posture_start_onset = posture_start_sample/data.fsample;
posture_stop_onset  = posture_stop_sample /data.fsample;

posture_start_curve = ones(size(posture_start_onset));
posture_stop_curve  = ones(size(posture_stop_onset ));

ax = figH.Children.SelectedTab.Children(3);
lim = ax.YLim(2);
hold(ax, 'on')

stem(ax, posture_start_onset, posture_start_curve*lim, 'red' )
stem(ax, posture_stop_onset , posture_stop_curve *lim, 'blue')
