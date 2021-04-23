%% Init

clear
clc

assert( ~isempty(which('ft_preprocessing')), 'FieldTrip library not detected. Check your MATLAB paths, or get : https://github.com/fieldtrip/fieldtrip' )
assert( ~isempty(which('farm_rootdir'))    ,      'FARM library not detected. Check your MATLAB paths, or get : https://github.com/benoitberanger/FARM' )

load e

fmri_volume = e.getSerie('run_nm').getVolume('^v').removeEmpty.getPath;
fmri_volume = sort(fmri_volume);
nVol = zeros(size(fmri_volume));
for i = 1 : length(fmri_volume)
    nii     = nifti(fmri_volume{i});
    nVol(i) = size(nii.dat,4);
end

subjdir       = e.getPath;
electrophydir = fullfile(subjdir,'electrophy');
filepath      = gfile(electrophydir,'run\d{2}.eeg$');
filepath      = cellstr(char(filepath));
filepath      = sort(filepath);

fname_eeg = filepath;
fname_hdr = spm_file(filepath, 'ext', '.vhdr');
fname_mrk = spm_file(filepath, 'ext', '.vmrk');

%%

t0 = tic;

for iRun = 1 : length(filepath)
    %% Get file & sequence paramters
    
    sequence        = [];
    sequence.TR     = 1.000; % in seconds
    sequence.nSlice = 72;
    sequence.MB     = 6;     % multiband factor
    sequence.nVol   = nVol(iRun);
    
    MRI_trigger_message = 'R128';
    
    emg_channel_regex = 'FCR|ECR|DEL|BIC|TRI';
    
    
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
    
    assert( nVol(iRun) <= numel(farm.sequence.get_volume_event(data)) )
    
    % Some paramters tuning
    data.cfg.intermediate_results_overwrite = 0; % don't overwrite files
    data.cfg.intermediate_results_save      = 1; % write on disk intermediate results
    data.cfg.intermediate_results_load      = 1; % if intermediate result file is detected, to not re-do step and load file
    
    % Plot
    % ft_databrowser(data.cfg, data)
    %     cfg.dataset
    %     [data.label num2cell(data.hdr.orig.impedances.channels)]
    % continue
    
    %% ------------------------------------------------------------------------
    %% FARM
    % Main FARM functions are below.
    
    % A lot of functions use what is called "regular expressions" (regex). It allows to recognize patterns in strings of characters
    % This a powerfull tool, which is common to almost all programing languages. Open some documentation with : doc regular-expressions
    
    fname = farm.io.mat.get_fname(data, 'pca_clean');
    if exist(fname, 'file')
        fprintf('existing file %s \n', fname)
        data = load(fname);
    else
        data = farm_main_workflow( data, emg_channel_regex );
        farm_export_BVA(data)
        farm_export_mat(data,[],1)
    end
    
    
    %% Some plots
    
    figH = farm_plot_FFT(data, emg_channel_regex, 'pca_clean', [30 250]  ); farm_print_figure( data, figH ); close(figH);
    figH = farm_plot_FFT(data,             'ACC',       'raw', [ 2   8],2); farm_print_figure( data, figH ); close(figH);
    
    
    %% 1 ACC or 1 ACC ?
    
    is_1_acc = any(strcmp( data.label, 'ACC_X'));
    is_2_acc = any(strcmp( data.label, 'R_ACC_X'));
    assert(xor(is_1_acc,is_2_acc),'wtf ?') % just to check
    
    if is_1_acc
        %% Time-Frequency Analysis
        
        cfg_TFA = [];
        cfg_TFA.emg_regex = emg_channel_regex;
        cfg_TFA.acc_regex = '^ACC';
        
        TFA = farm_time_frequency_analysis_emg_acc( data, cfg_TFA );
        figH = farm_plot_TFA( data, TFA ); farm_print_figure( data, figH ); close(figH);
        
        
        %% Coherence Analysis
        
        cfg_coh = [];
        cfg_coh.emg_regex = emg_channel_regex;
        cfg_coh.acc_regex = '^ACC';
        
        coh = farm_coherence_analysis_emg_acc( data, cfg_coh );
        % ft_connectivityplot([], coh);
        figH = farm_plot_coherence( data, coh ); farm_print_figure( data, figH ); close(figH);
        
        
        %% Select best EMG channel, that matches ACC using coherence
        
        cfg_select_emg = [];
        cfg_select_emg.emg_regex = emg_channel_regex;
        cfg_select_emg.acc_regex = '^ACC';
        
        best_emg = farm_select_best_emg_using_acc_coherence( data, cfg_select_emg );
        
        
        %% Generate regressors
        
        reginfo      = farm_make_regressor( data, best_emg.peakpower, best_emg.fsample);
        reginfo.name = ['peakpower@bestemg==' best_emg.label];
        figH         = farm_plot_regressor( data, reginfo ); farm_print_figure( data, figH ); close(figH);
        farm_save_regressor(data, reginfo)
        
        
        %% Accelerometer : this regressor will be a backup in case of bad EMG
        
        acc          = farm_get_timeseries(data,'^ACC','raw', [2 8],2);
        reginfo      = farm_acc_regressor(data, acc);
        reginfo.name = 'euclidiannorm@ACCXYZ_R';
        figH         = farm_plot_regressor( data, reginfo ); farm_print_figure( data, figH ); close(figH);
        farm_save_regressor(data, reginfo)
        
        
    else
        %% Time-Frequency Analysis
        
        cfg_TFA = [];
        cfg_TFA.emg_regex = emg_channel_regex;
        cfg_TFA.acc_regex = 'ACC';
        
        TFA = farm_time_frequency_analysis_emg_acc( data, cfg_TFA );
        figH = farm_plot_TFA( data, TFA ); farm_print_figure( data, figH ); close(figH);
        
        
        %% Coherence Analysis
        
        cfg_coh = [];
        cfg_coh.emg_regex = emg_channel_regex;
        cfg_coh.acc_regex = 'ACC';
        
        coh = farm_coherence_analysis_emg_acc( data, cfg_coh );
        % ft_connectivityplot([], coh);
        figH = farm_plot_coherence( data, coh ); farm_print_figure( data, figH ); close(figH);
        
        side = {'L', 'R'};
        
        for s = side
            %% Select best EMG channel, that matches ACC using coherence
            
            LR = char(s);
            
            cfg_select_emg = [];
            cfg_select_emg.emg_regex = sprintf('%s_(%s)',LR,emg_channel_regex);
            cfg_select_emg.acc_regex = sprintf('%s_ACC' ,LR                  );
            
            best_emg = farm_select_best_emg_using_acc_coherence( data, cfg_select_emg );
            
            
            %% Generate regressors
            
            reginfo      = farm_make_regressor( data, best_emg.peakpower, best_emg.fsample);
            reginfo.name = ['peakpower@bestemg==' best_emg.label];
            figH         = farm_plot_regressor( data, reginfo ); farm_print_figure( data, figH ); close(figH);
            farm_save_regressor(data, reginfo)
            
            
            %% Accelerometer : this regressor will be a backup in case of bad EMG
            
            acc          = farm_get_timeseries(data,sprintf('%s_ACC' ,LR),'raw', [2 8],2);
            reginfo      = farm_acc_regressor(data, acc);
            reginfo.name = sprintf('euclidiannorm@ACCXYZ_%s' ,LR);
            figH         = farm_plot_regressor( data, reginfo ); farm_print_figure( data, figH ); close(figH);
            farm_save_regressor(data, reginfo)
            
            
        end % LR
        
    end % 1 or 2 acc ?
    
    toc(t0);
    
end % run file

toc(t0);
