function [tlck_actv,tlck_bsln,nt] = ReloadDataOptionsForVS(mycfg)

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


if strcmp(mycfg.trigger,lower('emg'))
    % Reading EMG onset to set trial definition - this is a bit long
    % winded....
    %--------------------------------------------------------------------
    hdr  = ft_read_header(dataset);
    chanindx    = strmatch('EMG', hdr.label);
    WindowToSearch = [-600 1200]; % Search around the EMG Marker for EMG Onsets
    sdThreshold = 2.5;
    
    if length(chanindx)>1
        error('only one EMG channel supported');
    end
    
    cfg = [];
    cfg.dataset = dataset;
    cfg.trialdef.eventtype      = 'stim_off';
    event=ft_definetrial(cfg);
    
    % read all data of the EMG channel, assume continuous file format
    emg = ft_read_data(cfg.dataset, 'header', hdr, ...
        'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials, ...
        'chanindx', chanindx, 'checkboundary', false);
    
    % apply filtering, hilbert transformation and boxcar convolution (for smoothing)
    emgflt      = ft_preproc_bandpassfilter(emg, hdr.Fs, [15 150]); % bandpassfilter
    emgnth      = ft_preproc_bandstopfilter(emgflt,hdr.Fs, [48 52]);% notch filter
    emgrect     = abs(emgnth);
    
    mark = [];
    for i = 1 : length(event.trl)
        ss = event.trl(i,1) + (WindowToSearch (1)); %Get the first sample to search
        es = event.trl(i,1) + (WindowToSearch (2)); %Get the Last sample to search
        RTdiff(i) = inf;
        
        noise = mean(emgrect(ss:es));
        noisesd = std(emgrect(ss:es));
        noisethresh = (noise + (sdThreshold * noisesd));
        absmax = max(abs(emgnth(ss:es)));
        
        for j = ss : es
            if emgrect(j) > noisethresh
                mark = [mark j];
                RTdiff(i) = (j - event.trl(i,1));
                break
            end
        end
    end
    
    % make a new set of events
    for i = 1:length(mark)
        ev(i).type   = 'EMG_trig';
        ev(i).sample = mark(i);
        ev(i).value  = 1;
        ev(i).duration = [];
        ev(i).offset = [];
    end
    
    cfg                     = [];
    cfg.dataset             = dataset;    
    cfg.trialdef.prestim    = abs(bsln_toi(1));
    cfg.trialdef.poststim   = abs(actv_toi(2));
    
    event=ev;
    trl=[];
    for i=1:length(event)
                % add this to the trl definition
                begsample     = event(i).sample - cfg.trialdef.prestim*hdr.Fs;
                endsample     = event(i).sample + cfg.trialdef.poststim*hdr.Fs - 1;
                offset        = -cfg.trialdef.prestim*hdr.Fs;
                trigger       = event(i).value; % remember the trigger (=condition) for each trial
                if isempty(trl)
                    prevtrigger = nan;
                else
                    prevtrigger   = trl(end, 4); % the condition of the previous trial
                end
                trl(end+1, :) = [round([begsample endsample offset])  trigger prevtrigger];
    end
    
    cfg.trl = trl;
    cfg_data                = ft_definetrial(cfg);
    
else
    if ismember(mycfg.trigger,trig_names)
        fprintf('Found trigger: %s',mycfg.trigger);
    end
    cfg                     = [];
    cfg.dataset             = dataset;
    cfg.trialdef.eventtype  = mycfg.trigger;
    cfg.trialdef.prestim    = abs(bsln_toi(1));   ...  * ^
    cfg.trialdef.poststim   = abs(actv_toi(2));   ...  * ^
    cfg_data                = ft_definetrial(cfg);
end


% % Now we know the triggers exist, define trials
% %--------------------------------------------------------------------------
% cfg                     = [];
% cfg.dataset             = dataset;
% cfg.trialdef.eventtype  = mycfg.trigger;
% cfg.trialdef.prestim    = abs(bsln_toi(1));   ...  * ^
% cfg.trialdef.poststim   = abs(actv_toi(2));   ...  * ^
% cfg_data                = ft_definetrial(cfg);

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
    cfg.artfctdef.zvalue.cutoff      = 12;%5
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

nt = [size(cfg_data.trl,1) size(cfg_data.trlold,1)];

if isfield(mycfg,'reject_only') && mycfg.reject_only;
    fprintf('Reject only flagged: returning...\n');
    [tlck_actv,tlck_bsln]=deal([]);
    return
end

% Pre-process the data
%--------------------------------------------------------------------------
cfg_data.channel    = {'meg'};%{'-MRCNT*','-*STAT*','-MP*','-MM*','-MRSYN*'}';%{'MEG'};
cfg_data.precision  = 'single';
cfg_data.demean     = 'yes';
cfg_data.bpfilter   = 'yes';
cfg_data.bpfreq     = [1 100];
cfg_data.bpfilttype = 'but';
cfg_data.bpfiltdir  = 'twopass';
data_preproc        = ft_preprocessing(cfg_data);


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



