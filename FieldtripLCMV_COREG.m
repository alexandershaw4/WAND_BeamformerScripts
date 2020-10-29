function FieldtripLCMV_COREG(mycfg)


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

% If pos (headfile) exist, just strep fid locs from it, otherwise manual
if ~isempty(head_filename)
    %fid_data = ft_read_headshape(head_filename,'format','polhemus_pos');
    %fid_data = ft_convert_units(fid_data,'cm');
    
    cfg = [];
    cfg.method = 'fiducial';
    cfg.coordsys = 'ctf';
    cfg.spmversion = 'spm12';
    
    %nasid = find(contains(lower(fid_data.fid.label),'nas'));
    %cfg.fiducial.nas = fid_data.fid.pos(nasid,:);
    
    %lpaid = find(contains(lower(fid_data.fid.label),'lpa'));
    %cfg.fiducial.lpa = fid_data.fid.pos(lpaid,:);

    %rpaid = find(contains(lower(fid_data.fid.label),'rpa'));
    %cfg.fiducial.rpa = fid_data.fid.pos(rpaid,:);
   
    cfg.method = 'interactive';
    mri = ft_volumerealign(cfg,mri);
    
else
    % Co-register the MRI manually first
    cfg          = [];
    cfg.method   = 'interactive';
    cfg.coordsys = 'ctf';
    mri = ft_volumerealign(cfg,mri);
end

if ~isempty(head_filename)
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
end

% note - to avoid a popup/interactive window, i had to set line 52 of 
% ft_determine_coordsys to interactive = no

% ft_plot_headshape(head_data);

%cfg = [];
%mri = ft_volumereslice(cfg, mri);

% Save the coregistered mri - in the MEG preproc for now?!
[~,mri_name,~] = fileparts(mycfg.MRI);
save([SaveDir mri_name '_CoReg'], '-struct', 'mri')
