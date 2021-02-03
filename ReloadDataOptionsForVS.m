function [tlck_actv,tlck_bsln] = ReloadDataOptionsForVS(mycfg)

CleanFieldtripFromPaths;



dataset = [mycfg.subject_dir '/' mycfg.MEG];

% Full path to the fieldtrip toolbox
%--------------------------------------------------------------------------
fieldtrip_path = '/cubric/software/MEG/fieldtrip-20190219/';
addpath(fieldtrip_path); cd(fieldtrip_path); ft_defaults

    
%  loading the data, defining triggers ... 
%==========================================================================


% Covariance frequencies and time-range
%--------------------------------------------------------------------------
cov_foi = mycfg.cov_fwin; %in Hz
cov_toi = mycfg.cov_toi;  %in seconds

% Stimulus and baseline time-range (for beamformer contrast)
%--------------------------------------------------------------------------
bsln_toi = mycfg.b_toi; %baseline
actv_toi = mycfg.a_toi; %stimulus


%This does a check of what triggers are in the dataset
%--------------------------------------------------------------------------
cfg                     = [];
cfg.dataset             = dataset;
cfg.trialdef.eventtype  = '?';
trigs                   = ft_definetrial(cfg);

trig_names = unique({trigs.event.type}); % Filter list

if ismember(mycfg.trigger,trig_names)
    fprintf('Found trigger: %s',mycfg.trigger);
end


% Now we know the triggers exist, define trials
%--------------------------------------------------------------------------
cfg                     = [];
cfg.dataset             = dataset;
cfg.trialdef.eventtype  = mycfg.trigger;
cfg.trialdef.prestim    = abs(bsln_toi(1));   ...  * ^
cfg.trialdef.poststim   = abs(actv_toi(2));   ...  * ^
cfg_data                = ft_definetrial(cfg);

% *^Ok, maybe this could be a different parameter, but by deault I assume
% you'd want the length of a given trial to be beginning-of-baseline to end
% -of-active period... seems reasonable.


% NEW: ARTIFACT REJECT
if mycfg.Artifact
    
    % jump artifacts
    %--------------------------------------------------------------------------
    cfg            = [];
    cfg.trl        = cfg_data.trl;
    cfg.datafile   = dataset;
    cfg.headerfile = dataset;
    cfg.continuous = 'yes';
    
    % channel selection, cutoff and padding
    %cfg.artfctdef.zvalue.channel = 'MEG';
    cfg.artfctdef.zvalue.channel = {'MEG','-MRCNT*','-*STAT*','-MP*','-MM*','-MRSYN*'};
    cfg.artfctdef.zvalue.cutoff  = 35;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative = 'yes';
    cfg.artfctdef.zvalue.medianfilter = 'yes';
    cfg.artfctdef.zvalue.medianfiltord = 9;
    cfg.artfctdef.zvalue.absdiff = 'yes';
    
    % identify
    cfg.artfctdef.zvalue.interactive = 'no';
    [cfg, artifact_jump] = ft_artifact_zvalue(cfg);
    
    % EOG artifacts
    %--------------------------------------------------------------------------
    cfg            = [];
    cfg.trl        = cfg_data.trl;
    cfg.datafile   = dataset;
    cfg.headerfile = dataset;
    cfg.continuous = 'yes';
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel     = 'EOG';
    cfg.artfctdef.zvalue.cutoff      = 5;
    cfg.artfctdef.zvalue.trlpadding  = 0;
    cfg.artfctdef.zvalue.artpadding  = 0.1;
    cfg.artfctdef.zvalue.fltpadding  = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter   = 'yes';
    cfg.artfctdef.zvalue.bpfilttype = 'but';
    cfg.artfctdef.zvalue.bpfreq     = [2 20];
    cfg.artfctdef.zvalue.bpfiltord  = 4;
    cfg.artfctdef.zvalue.hilbert    = 'yes';
    
    % identify
    cfg.artfctdef.zvalue.interactive = 'no';
    [cfg, artifact_EOG] = ft_artifact_zvalue(cfg);
    
    % now remove (nan)
    cfg_data.artfctdef.reject          = 'complete';
    cfg_data.artfctdef.feedback        =  'no' ;
    cfg_data.artfctdef.eog.artifact    = artifact_EOG;
    cfg_data.artfctdef.jump.artifact   = artifact_jump;
    
    [cfg_data] = ft_rejectartifact(cfg_data);
    
end

% Pre-process the data
%--------------------------------------------------------------------------
cfg_data.channel    = {'meg'};%{'-MRCNT*','-*STAT*','-MP*','-MM*','-MRSYN*'}';%{'MEG'};
cfg_data.precision  = 'single';
cfg_data.demean     = 'yes';
cfg_data.bpfilter   = 'yes';
cfg_data.bpfreq     = cov_foi;
cfg_data.bpfilttype = 'but';
cfg_data.bpfiltdir  = 'twopass';
data_preproc        = ft_preprocessing(cfg_data);

% Compute the covariance
%--------------------------------------------------------------------------
cfg = [];
cfg.channel = {'MEG'};
cfg.removemean = 'no';
cfg.covariance = 'yes';
cfg.covariancewindow = cov_toi;% [-1.5 1.5];
data_tlck = ft_timelockanalysis(cfg, data_preproc);



% Define baseline and stimulus epochs  -
%   we will then beam these using the common filter computed above
%--------------------------------------------------------------------------
cfg        = [];
cfg.toilim = bsln_toi;
data_bsln  = ft_redefinetrial(cfg, data_preproc);
cfg        = [];
cfg.toilim = actv_toi;
data_actv  = ft_redefinetrial(cfg, data_preproc); 

% Timelock analysis of bandpass-filtered baseline and active windows
cfg             = [];
cfg.channel     = {'MEG'};
cfg.keeptrials  = 'yes';
cfg.covariance  = 'yes'; 
tlck_bsln       = ft_timelockanalysis(cfg, data_bsln);
tlck_actv       = ft_timelockanalysis(cfg, data_actv);



