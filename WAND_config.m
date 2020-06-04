CleanFieldtripFromPaths

wand    = '/cubric/collab/314_wand/bids/sourcedata/';
subject = '314_05117';
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
heads             = dir([wand subject '/meg/generic/*.pos']);
mycfg.Headshape   = ['/meg/generic/' heads(1).name];
mycfg.SaveSubDir  = [wand subject '/meg/' task '/preproc/'];
mycfg.sourcemodel = []; % (empty = use 4mm, or spec fullfile to .mat)


% You can use this code to just compute the leadfields and then stop 
% (saving everything into the preproc dir for later use), or you
% can set ComputeBeamformer = 1 and specify your data / contrast parameters
% below.

mycfg.ComputeBeamformer = 0;

% Optional: Only needed if mycfg.ComputeBeamformer = 1 ... 
%==========================================================================


% Define stimulus & frequency parameters
%--------------------------------------------------------------------------
mycfg.cov_fwin = [40 80];
mycfg.cov_toi  = [-1.5 1.5];

mycfg.trigger  = 'stim_on';  % trig/stim 
mycfg.b_toi    = [-1.2 0];   % baseline times
mycfg.a_toi    = [0.3 1.5];  % active times

mycfg.UseKrish1mm = 0;


% Other options: Epoch or Clean (uncomment as required)
%==========================================================================

DoEpoch    = false;
DoArtifact = false;

if DoEpoch
    
    % Epoch: Want to epoch using CTF tools (i.e. actually epoch the .ds)?
    %-------------------------------------------------------------------------
    trigname = 'stim_on';
    times    = [-2 2];
    dataset  = [mycfg.subject_dir mycfg.MEG];
    aMEG.epoch.epocher({dataset},trigname,times);
    [fp,fn,fe] = fileparts(mycfg.MEG); 
    mycfg.MEG  = [fp '/' fn 'Cut' fe];
end

if DoArtifact
    
    % Cleaning: Reject jump & EOG artifacts from continuous data 
    % note: this then finds the trials corresponding to bad segs and writes a
    % new .ds dataset
    %--------------------------------------------------------------------------
    trigname = 'stim_on';
    times    = [-2 2];
    dataset  = [mycfg.subject_dir mycfg.MEG];
    ftclean_ctf(dataset,trigname,times);

    % Fix segmentation violation fault
    bfix = [fileparts(mfilename('fullpath')) '/fixsve.sh'];
    system(['sh ' bfix ' "' dataset '"']);

    [fp,fn,fe] = fileparts(mycfg.MEG); 
    mycfg.MEG  = [fp '/' fn 'Fix' fe];
end

% Compute the leadsfields and/or beamformer
FieldtripLCMVBeamformer(mycfg);






