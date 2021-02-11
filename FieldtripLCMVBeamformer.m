function data = FieldtripLCMVBeamformer(mycfg)

CleanFieldtripFromPaths;

load(mycfg);

% Full path to the fieldtrip toolbox
%--------------------------------------------------------------------------
fieldtrip_path = '/cubric/software/MEG/fieldtrip-20190219/';
addpath(fieldtrip_path); cd(fieldtrip_path); ft_defaults

% Reconstruct full file names
%--------------------------------------------------------------------------
sub_folder    = mycfg.subject_dir; 
dataset       = [mycfg.subject_dir '/' mycfg.MEG];
mri_filename  = [mycfg.subject_dir '/' mycfg.MRI];

if ~isempty(mycfg.Headshape)
    head_filename = [mycfg.subject_dir '/' mycfg.Headshape];
else
    head_filename = [];
end

% Use default 4mm source grid if no .mat specified...
%--------------------------------------------------------------------------
if ~isfield(mycfg,'sourcemodel') || isempty(mycfg.sourcemodel)
      sourcemodel = 'standard_sourcemodel3d4mm.mat';
else; sourcemodel = mycfg.sourcemodel;
end

beam_method = 'lcmv';
cd(sub_folder);

% Where to save everything
SaveDir = mycfg.SaveSubDir;

% save outputs with a prepended name - i.e. so desynch anal doesnt
% overwrite rebound
if isfield(mycfg,'prepend')
    %mycfg.type = 'desync';
    prepend = mycfg.prepend;
else
    prepend = [];
end

% Are we running the beamformer or just compute and saving the leadfields,
% head & source models?
RunBeamformer = mycfg.ComputeBeamformer;


% Shouldn't need to change much passed this line...
%==========================================================================

% Load & Co-register the MRI
%--------------------------------------------------------------------------
[~,mri_name,~] = fileparts(mycfg.MRI);
mri = load([SaveDir mri_name '_CoReg']);

if exist([SaveDir, mri_name '_Segmented.mat'])
   fprintf('Loading existing segmented mri...\n');
   mri_segmented = load([SaveDir, mri_name '_Segmented']);
else
    % Segment the MRI
    %--------------------------------------------------------------------------
    cfg           = [];
    cfg.output    = {'brain'; 'skull'; 'scalp'};
    mri_segmented = ft_volumesegment(cfg, mri);
    mri_segmented.anatomy = mri.anatomy;

    % Save segmented mri
    save([SaveDir, mri_name '_Segmented'], '-struct', 'mri_segmented')
end

% Compute & Save Headmodel
%--------------------------------------------------------------------------
cfg        = [];
cfg.method = 'singleshell';
hdm = ft_prepare_headmodel(cfg, mri_segmented);
hdm = ft_convert_units(hdm, 'cm'); 

% Save headmodel
save([SaveDir, 'headmodel'], '-struct', 'hdm')

if exist([SaveDir, 'sourcemodel_template.mat']) 
    fprintf('Loading existing sourcemodel...\n');
    template_sourcemodel = load([SaveDir, 'sourcemodel_template']) ;
else
    % Compute & Save Sourcemodel
    %--------------------------------------------------------------------------
    if ~mycfg.UseKrish1mm

        %load(sourcemodel)
        %tmp = whos('-file',sourcemodel);
        %template_sourcemodel = eval(tmp.name);
        load(sourcemodel)
        template_sourcemodel = ft_convert_units(sourcemodel, 'cm');
        %clear(tmp.name,'tmp')

        % Save template sourcemodel
        save([SaveDir, 'sourcemodel_template'], '-struct', 'template_sourcemodel')

    elseif mycfg.UseKrish1mm

        addpath('/cubric/scratch/krish/NewAAL/');
        fprintf('Using Krish 1mm AAL90-Masked Source model\n');
        AtlasMNI='/cubric/software/MEG/fieldtrip-20161011/template/atlas/aal/ROI_MNI_V4.nii';
        AtlasLabels={'Calcarine_L' 'Calcarine_R' 'Cuneus_L' 'Cuneus_R' 'Lingual_L' 'Lingual_R' 'Occipital_Sup_L' 'Occipital_Sup_R' 'Occipital_Mid_L' 'Occipital_Mid_R' 'Occipital_Inf_L' 'Occipital_Inf_R' 'Fusiform_L' 'Fusiform_R'};
        %AtlasLabels = 'AllAtlasLabels';
        template =kMakeCustomTemplateSourceGrid([-6.4 6.4],[-11.5 -5.7],[-2.9 4.8],0.1,AtlasMNI,AtlasLabels,true);
        %template =kMakeCustomTemplateSourceGrid([-6.4 6.4],[-11.5 -5.7],[-2.9 4.8],0.4,AtlasMNI,AtlasLabels,true);
        template_sourcemodel = ft_convert_units(template, 'mm');

        save([SaveDir, 'sourcemodel_template'], '-struct', 'template_sourcemodel')        
    end
end

if exist([SaveDir, 'sourcemodel_indiv.mat'])
    fprintf('Loading existing sourcemodel (individ)...\n');
    sourcemodel = load([SaveDir, 'sourcemodel_indiv']);
else
    % Compute (warp) & Save the Sourcemodel
    %--------------------------------------------------------------------------
    cfg = [];
    cfg.grid.warpmni    = 'yes';
    cfg.grid.template   = template_sourcemodel;
    cfg.grid.nonlinear  = 'yes';
    cfg.mri             = mri;
    % cfg.inwardshift   = -1.5;
    sourcemodel         = ft_prepare_sourcemodel(cfg);

    % Save sourcemodel
    save([SaveDir, 'sourcemodel_indiv'], '-struct', 'sourcemodel')
end

if exist([SaveDir, 'grad.mat'])
    fprintf('Loading existing gradiometer definitions...\n');
    grad = load([SaveDir, 'grad']);
else
    % Sensor positions: Get grad structure from header file
    %--------------------------------------------------------------------------
    hdr  = ft_read_header(dataset);
    grad = hdr.grad;
    grad = ft_convert_units(grad, 'cm'); %make sure units match...

    % Save grad
    save([SaveDir, 'grad'], '-struct', 'grad')
end


if exist([SaveDir, 'leadfield.mat'])
    fprintf('Loading existing leadfields...\n');
    leadfield = load([SaveDir, 'leadfield']);
else
    % Compute & Save the Leadfields
    %--------------------------------------------------------------------------
    cfg                 = [];
    cfg.channel         = {'MEG'};
    cfg.grid            = sourcemodel;
    cfg.headmodel       = hdm; 
    cfg.grad            = grad;
    cfg.normalize       = 'yes'; 
    leadfield           = ft_prepare_leadfield(cfg);

    % Save leadfield
    save([SaveDir, 'leadfield'], '-struct', 'leadfield')
end


if ~RunBeamformer
    % we're done - everything is saved, return...
    return;
end
    

% Otherwise, carry on to loading the data, defining triggers ... 
%==========================================================================

hdr  = ft_read_header(dataset);

% Covariance frequencies and time-range
%--------------------------------------------------------------------------
cov_foi = mycfg.cov_fwin; %in Hz
cov_toi = mycfg.cov_toi;  %in seconds

% Stimulus and baseline time-range (for beamformer contrast)
%--------------------------------------------------------------------------
bsln_toi = mycfg.b_toi; %baseline
actv_toi = mycfg.a_toi; %stimulus

% Check if we've already epoched, cleaned and saved the data in fieldtrip
if exist([SaveDir, prepend, 'TheData.mat'])
    load([SaveDir, prepend, 'TheData'], 'cfg_data')
    fprintf('Loading existing cleaned/epoched data...\n');
else
    

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
        %cfg.trialdef.prestim    = abs(bsln_toi(1));
        %cfg.trialdef.poststim   = abs(actv_toi(2));
        
        cfg.trialdef.prestim  = abs(mycfg.cov_toi(1));
        cfg.trialdef.poststim = abs(mycfg.cov_toi(2));

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




    % Now we know the triggers exist, define trials
    %--------------------------------------------------------------------------
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

        if isfield(mycfg,'eogthreshold')
            eogcut = mycfg.eogthreshold;
        else
            eogcut = 5;
        end
        
        if isfield(mycfg,'jumpthreshold')
            jumpcut = mycfg.jumpthreshold;
        else
            jumpcut = 35;
        end
        
        fprintf('Using EOG thresh: %d\n',eogcut);
        fprintf('Using Peak2Peak thresh: %d\n',jumpcut);

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
        cfg.artfctdef.zvalue.cutoff  = jumpcut;
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
        cfg.artfctdef.zvalue.cutoff      = eogcut;
        cfg.artfctdef.zvalue.trlpadding  = 0;
        cfg.artfctdef.zvalue.artpadding  = 0.1;
        cfg.artfctdef.zvalue.fltpadding  = 0;

        % algorithmic parameters
        cfg.artfctdef.zvalue.bpfilter   = 'yes';
        cfg.artfctdef.zvalue.bpfilttype = 'but';
        cfg.artfctdef.zvalue.bpfreq     = [2 15];
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

    if isfield(mycfg,'ManualReject') && mycfg.ManualReject
        cfg_data        = ft_preprocessing(cfg_data);
        
        cfg          = [];
        cfg.method   = 'trial'; % 'trial';
        cfg.ylim     = [-2e-12 2e-12];
        cfg.channel = {'MEG'};
        
        cfg.metric  = 'zvalue';
        
        cfg.preproc.bpfilter    = 'yes';
        cfg.preproc.bpfreq      = [2 70];
        cfg.preproc.bpfiltord   =  4;
        cfg.preproc.bpfilttype  = 'but';
        %cfg.preproc.rectify     = 'yes';
        %cfg.preproc.boxcar      = 0.2;
        cfg_data     = ft_rejectvisual(cfg,cfg_data);
    end

    if any(prepend)
        save([SaveDir, prepend, 'TheData'], 'cfg_data')
    else
        save([SaveDir, 'TheData'], 'cfg_data');
    end

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









% Compute common weights
%--------------------------------------------------------------------------
cfg                             = [];
cfg.method                      = beam_method;% 'lcmv';
cfg.grid                        = leadfield;
cfg.headmodel                   = hdm;
cfg.grad                        = grad;
cfg.(beam_method).fixedori      = 'yes';
cfg.(beam_method).keepfilter    = 'yes';
cfg.(beam_method).projectnoise  = 'yes';
cfg.(beam_method).lambda        = '5%'; 
src                             = ft_sourceanalysis(cfg, data_tlck);
src                             = rmfield(src,'cfg'); 

if any(prepend)
    save([SaveDir, prepend, 'CommonWeights'], '-struct', 'src')
else
    save([SaveDir, 'CommonWeights'], '-struct', 'src')
end

% Get the (common) beamformer weights & save on their own for easy access
wts = src.avg.filter;
if any(prepend)
    save([SaveDir, prepend, 'CommonwightsMat'], 'wts')
else
    save([SaveDir, 'CommonwightsMat'], 'wts')
end

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


% Get source power estimates separately for bsln and actv - using the
% common weights
%--------------------------------------------------------------------------
cfg                             = [];
cfg.method                      = beam_method;
cfg.grid                        = leadfield;   % precomputed leadfield
cfg.grid.filter                 = wts;         % precomputed common weights
cfg.headmodel                   = hdm;         % headmodel
cfg.grad                        = grad;        % gradiomter positions
cfg.keeptrials                  = 'yes';
cfg.(beam_method).fixedori      = 'yes';
cfg.(beam_method).keepfilter    = 'no';
cfg.(beam_method).projectnoise  = 'yes';
cfg.(beam_method).lambda        = '5%'; 

% Pass the baseline data through 
src_bsln = ft_sourceanalysis(cfg, tlck_bsln);
src_bsln = rmfield(src_bsln,'cfg'); 
src_bsln.pos_template = template_sourcemodel.pos;
src_bsln.dim_template = template_sourcemodel.dim;

% Pass the stimulus data through 
src_actv = ft_sourceanalysis(cfg, tlck_actv);
src_actv = rmfield(src_actv,'cfg');
src_actv.pos_template = template_sourcemodel.pos;
src_actv.dim_template = template_sourcemodel.dim;

% Compute the actual CONTRAST:
% - Difference in power as a percentage change from baseline
%--------------------------------------------------------------------------
cfg           = [];
cfg.parameter = 'avg.pow';
cfg.operation = '((x1-x2)./x2)*100';
src_diff      = ft_math(cfg, src_actv, src_bsln);
src_diff.pos_native   = src_diff.pos;
src_diff.pos_template = template_sourcemodel.pos;

% Save source estimates...
if any(prepend)
    save([SaveDir, prepend, 'SourceBaseline'], '-struct', 'src_bsln')
    save([SaveDir, prepend, 'SourceActive']  , '-struct', 'src_actv')
    save([SaveDir, prepend, 'SourceContrast'], '-struct', 'src_diff')
else
    save([SaveDir, 'SourceBaseline'], '-struct', 'src_bsln')
    save([SaveDir, 'SourceActive']  , '-struct', 'src_actv')
    save([SaveDir, 'SourceContrast'], '-struct', 'src_diff')
end

% Interpolate (in native space) & Save
%--------------------------------------------------------------------------
cfg              = [];
cfg.parameter    = {'pow'};
cfg.interpmethod = 'nearest';
src_intrp_indiv  = ft_sourceinterpolate(cfg, src_diff, mri);

if any(prepend)
    save([SaveDir, prepend, 'SourceInterpIndiv'],'-struct','src_intrp_indiv');
else
    save([SaveDir, 'SourceInterpIndiv'],'-struct','src_intrp_indiv');    
end

% Interpolate (in template MNI space) & save
%--------------------------------------------------------------------------
templatemri_file = fullfile(fieldtrip_path, 'template/anatomy/single_subj_T1_1mm.nii');
templatemri      = ft_read_mri(templatemri_file);

cfg               = [];
cfg.parameter     = {'pow'};
cfg.interpmethod  = 'nearest';
src_intr          = ft_sourceinterpolate(cfg, src_diff, templatemri);
src_intr.coordsys = 'mni'; 
src_intr          = ft_convert_units(src_intr, 'mm'); 

if any(prepend)
    save([SaveDir, prepend, 'SourceInterpTemplate'],'src_intr');
else
    save([SaveDir, 'SourceInterpTemplate'],'src_intr');    
end

end

