function WAND_config(subject)

%CleanFieldtripFromPaths

wand    = '/cubric/collab/314_wand/bids/sourcedata/';
%subject = '314_05117';
task    = 'visual';    

% task = 'visual' 'auditorymotor' 'day3' 
%        'dotback' 'mmn' 'resting' 'simon' 'sternberg'


% Default data locations (in BIDS): MEG, MRI & Headshape
%--------------------------------------------------------------------------
mycfg.subject_dir = [wand subject];
MEG               = dir([wand subject '/meg/' task '/raw/*.ds']);
mycfg.MEG         = ['/meg/' task '/raw/' MEG(1).name]; 
niftis            = dir([wand subject '/mri/anat/*.nii']);
mycfg.MRI         = ['/mri/anat/' niftis(1).name];

try
    heads             = dir([wand subject '/meg/generic/*.pos']);
    mycfg.Headshape   = ['/meg/generic/' heads(1).name];
catch
    mycfg.Headshape   = [];
end

mycfg.SaveSubDir  = [wand subject '/meg/' task '/preproc/'];
mycfg.sourcemodel = []; % (empty = use 4mm, or spec fullfile to .mat)

% Alex temporary hack - can't save to the collab atm so make a home in
% scratch
mycfg.SaveSubDir = strrep(mycfg.SaveSubDir,'/cubric/collab/',...
    '/cubric/scratch/sapas10/');
unix(['mkdir -p ' mycfg.SaveSubDir]);


% You can use this code to just compute the leadfields and then stop 
% (saving everything into the preproc dir for later use), or you
% can set ComputeBeamformer = 1 and specify your data / contrast parameters
% below.

mycfg.ComputeBeamformer = 1;

% Optional: Only needed if mycfg.ComputeBeamformer = 1 ... 
%==========================================================================


% Define stimulus & frequency parameters
%--------------------------------------------------------------------------
mycfg.cov_fwin = [30 80];
mycfg.cov_toi  = [0 1.5];

mycfg.trigger  = 'stim_on';  % trig/stim 
mycfg.b_toi    = [-1.2 0];   % baseline times
mycfg.a_toi    = [0.3 1.5];  % active times

mycfg.UseKrish1mm = 1;
mycfg.Artifact = true;

% the manual local part - coregistration
if isempty(dir([mycfg.SaveSubDir '*CoReg.mat']))
    FieldtripLCMV_COREG(mycfg);
end

% save the input struct to load ont cluster
t = clock;
unique_name = [date '_' num2str(t(5)) num2str(round(t(6)))];
unique_name = [mycfg.SaveSubDir unique_name];

save(unique_name,'mycfg')

% Compute the lead fields and beamformer on the cluster
fprintf('Submitting job to cluster...\n');
cd(mycfg.SaveSubDir);
docluster_slurm_bfm('FieldtripLCMVBeamformer',unique_name);

% or run locally:
% FieldtripLCMVBeamformer(unique_name)




