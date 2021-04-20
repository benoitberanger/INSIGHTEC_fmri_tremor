%% Init

clear
clc

assert( ~isempty(which('ft_preprocessing')), 'FieldTrip library not detected. Check your MATLAB paths, or get : https://github.com/fieldtrip/fieldtrip' )
assert( ~isempty(which('farm_rootdir'))    ,      'FARM library not detected. Check your MATLAB paths, or get : https://github.com/benoitberanger/FARM' )

load e
e = e(end);

fmri_volume = e.getSerie('run_nm').getVolume('sw').removeEmpty.getPath;
nVol = zeros(size(fmri_volume));
for i = 1 : length(fmri_volume)
    nii     = nifti(fmri_volume{i});
    nVol(i) = size(nii.dat,4);
end

subjdir       = e.getPath;
electrophydir = fullfile(subjdir,'electrophy');
filepath      = gfile(electrophydir,'run\d{2}.eeg$');
filepath      = cellstr(char(filepath));

fname_eeg = filepath;
fname_hdr = spm_file(filepath, 'ext', '.vhdr');
fname_mrk = spm_file(filepath, 'ext', '.vmrk');

for iRun = 1 : length(filepath)
    %% Get file & sequence paramters
    
    sequence.TR     = 1.000; % in seconds
    sequence.nSlice = 72;
    sequence.MB     = 6;     % multiband factor
    sequence.nVol   = nVol(iRun);
    
    MRI_trigger_message = 'R128';
    
    emg_channel_regex = 'FCR|ECR|DEL|BIC|TRI';
%     emg_channel_regex = 'TRI';
    
    %% Load data
    % Optimal length for a dataset is a bunch of seconds before the start of
    % the fmri sequence, and a bunch of seconds after the end of the fmri
    % sequence, before any other sequence.
    
    % Read header & events
    cfg           = [];
    cfg.dataset   = fname_hdr{iRun};
    raw_event     = ft_read_event(fname_mrk{iRun});
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
    
    % Output directory
    % If no outdir is defined, use the same as inputdir
    outdir = '/network/lustre/iss01/cenir/analyse/irm/users/benoit.beranger/INSIGHTEC_fmri_tremor/test';
    data.cfg.outdir.intermediate = fullfile( outdir ); % intermediate results
    data.cfg.outdir.BVAexport    = fullfile( outdir ); % export final results in {.eeg, .vhdr, .vmrk}
    data.cfg.outdir.MATexport    = fullfile( outdir ); % export final results in .mat
    data.cfg.outdir.png          = fullfile( outdir ); % write PNG here, for visual quick check
    data.cfg.outdir.regressor    = fullfile( outdir ); % write regressor here, in .mat

    
    % Plot
    % ft_databrowser(data.cfg, data)
    
    
    %% ------------------------------------------------------------------------
    %% FARM
    % Main FARM functions are below.
    
    % A lot of functions use what is called "regular expressions" (regex). It allows to recognize patterns in strings of characters
    % This a powerfull tool, which is common to almost all programing languages. Open some documentation with : doc regular-expressions
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [ datapoints, channel_idx, channel_name, stage ] = farm.plot.get_datapoints( data, emg_channel_regex, 'raw' );
%     filter = -500;
%     datapoints = farm.filter(datapoints, data.fsample, filter);
%     data.trial{1}(channel_idx,:) = datapoints;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     data = farm_main_workflow( data, emg_channel_regex );
%     farm_plot_FFT(data, emg_channel_regex, 'pca_clean', [30 250]  );
%     return
    


end