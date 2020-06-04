function data = FieldtripLCMVBeamformer(mycfg)


% Full path to the fieldtrip toolbox
%--------------------------------------------------------------------------
fieldtrip_path = '/home/sapas10/fieldtrip-20170509/';
addpath(fieldtrip_path); cd(fieldtrip_path); ft_defaults

% Reconstruct full file names
%--------------------------------------------------------------------------
sub_folder    = mycfg.subject_dir; 
dataset       = [mycfg.subject_dir '/' mycfg.MEG];
mri_filename  = [mycfg.subject_dir '/' mycfg.MRI];
head_filename = [mycfg.subject_dir '/' mycfg.Headshape];

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

% Are we running the beamformer or just compute and saving the leadfields,
% head & source models?
RunBeamformer = mycfg.ComputeBeamformer;



% Shouldn't need to change much passed this line...
%==========================================================================

% Load & Co-register the MRI
%--------------------------------------------------------------------------

% Now assumes you passed a nifti in analyse format
mri = ft_read_mri(mri_filename,'dataformat','nifti');
mri = ft_convert_units(mri,'cm');

% Co-register the MRI manually first
cfg          = [];
cfg.method   = 'interactive';
cfg.coordsys = 'ctf';
mri = ft_volumerealign(cfg,mri);

% Read headshape data
head_data = ft_read_headshape(head_filename,'format','polhemus_pos');
head_data = ft_convert_units(head_data,'cm');

% ft_determine_coordsys(mri, 'interactive', 'no')
% ft_plot_headshape(head_data);

cfg            = [];
cfg.method                = 'headshape';
cfg.headshape.interactive = 'no'; 
cfg.headshape.icp         = 'yes'; 
cfg.headshape.headshape   = head_data;
cfg.coordsys              = 'ctf'; % 'ctf';
cfg.spmversion            = 'spm12';
cfg.viewresult            = 'yes';
mri = ft_volumerealign(cfg, mri); 

% note - to avoid a popup/interactive window, i had to set line 52 of 
% ft_determine_coordsys to interactive = no

% ft_plot_headshape(head_data);

cfg = [];
mri = ft_volumereslice(cfg, mri);

% Save the coregistered mri - in the MEG preproc for now?!
[~,mri_name,~] = fileparts(mycfg.MRI);
save([SaveDir mri_name '_CoReg'], '-struct', 'mri')


% %visualise mri?
% cfg = []; ft_sourceplot(cfg,mri)


% Segment the MRI
%--------------------------------------------------------------------------
cfg           = [];
cfg.output    = {'brain'; 'skull'; 'scalp'};
mri_segmented = ft_volumesegment(cfg, mri);
mri_segmented.anatomy = mri.anatomy;

% Save segmented mri
save([SaveDir, mri_name '_Segmented'], '-struct', 'mri_segmented')


% Compute & Save Headmodel
%--------------------------------------------------------------------------
cfg        = [];
cfg.method = 'singleshell';
hdm = ft_prepare_headmodel(cfg, mri_segmented);
hdm = ft_convert_units(hdm, 'cm'); 

% Save headmodel
save([SaveDir, 'headmodel'], '-struct', 'hdm')


% Compute & Save Sourcemodel
%--------------------------------------------------------------------------
load(sourcemodel)
tmp = whos('-file',sourcemodel);
template_sourcemodel = eval(tmp.name);
template_sourcemodel = ft_convert_units(template_sourcemodel, 'cm');
clear(tmp.name,'tmp')

% Save template sourcemodel
save([SaveDir, 'sourcemodel_template'], '-struct', 'template_sourcemodel')

% addpath('/cubric/scratch/krish/NewAAL/');
% fprintf('Using Krish 1mm AAL90-Masked Source model\n');
% AtlasMNI='/cubric/software/MEG/fieldtrip-20161011/template/atlas/aal/ROI_MNI_V4.nii';
% AtlasLabels={'Calcarine_L' 'Calcarine_R' 'Cuneus_L' 'Cuneus_R' 'Lingual_L' 'Lingual_R' 'Occipital_Sup_L' 'Occipital_Sup_R' 'Occipital_Mid_L' 'Occipital_Mid_R' 'Occipital_Inf_L' 'Occipital_Inf_R' 'Fusiform_L' 'Fusiform_R'};
% template =kMakeCustomTemplateSourceGrid([-6.4 6.4],[-11.5 -5.7],[-2.9 4.8],0.1,AtlasMNI,AtlasLabels,true);
% template_sourcemodel = ft_convert_units(template, 'mm');


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


% Sensor positions: Get grad structure from header file
%--------------------------------------------------------------------------
hdr  = ft_read_header(dataset);
grad = hdr.grad;
grad = ft_convert_units(grad, 'cm'); %make sure units match...

% Save grad
save([SaveDir, 'grad'], '-struct', 'grad')


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


if ~RunBeamformer
    % we're done - everything is saved, return...
    return;
end
    

% Otherwise, carry on to loading the data, defining triggers ... 
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

save([SaveDir, 'CommonWeights'], '-struct', 'src')

% Get the (common) beamformer weights & save on their own for easy access
wts = src.avg.filter;
save([SaveDir, 'CommonwightsMat'], 'wts')


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
save([SaveDir, 'SourceBaseline'], '-struct', 'src_bsln')
save([SaveDir, 'SourceActive']  , '-struct', 'src_actv')
save([SaveDir, 'SourceContrast'], '-struct', 'src_diff')


% Interpolate (in native space) & Save
%--------------------------------------------------------------------------
cfg              = [];
cfg.parameter    = {'pow'};
cfg.interpmethod = 'nearest';
src_intrp_indiv  = ft_sourceinterpolate(cfg, src_diff, mri);

save([SaveDir, 'SourceInterpIndiv'],'-struct','src_intrp_indiv');


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

save([SaveDir, 'SourceInterpTemplate'],'src_intr');


end

