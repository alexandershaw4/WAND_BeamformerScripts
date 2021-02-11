function WAND_config_motor(subject)

%CleanFieldtripFromPaths

wand    = '/cubric/collab/314_wand/bids/sourcedata/';
%subject = '314_05117';
task    = 'auditorymotor';    

runoncluster = 1;
mycfg.ComputeBeamformer = 1;
mycfg.UseKrish1mm = 0;
mycfg.Artifact = true;
mycfg.ManualReject=0;

% New - rejection thresholds
mycfg.eogthreshold  = 5;  % 5
mycfg.jumpthreshold = 20;  % 35

% task = 'visual' 'auditorymotor' 'day3' 
%        'dotback' 'mmn' 'resting' 'simon' 'sternberg'


% Default data locations (in BIDS): MEG, MRI & Headshape
%--------------------------------------------------------------------------
if ischar(subject)
    mycfg.subject_dir = [wand subject];
    MEG               = dir([wand subject '/meg/' task '/raw/*.ds']);
    mycfg.MEG         = ['/meg/' task '/raw/' MEG(1).name]; 
    niftis            = dir([wand subject '/mri/anat/*.nii']);
    mycfg.MRI         = ['/mri/anat/' niftis(1).name];
elseif isstruct(subject)
    mycfg.subject_dir = [wand subject.id];
    MEG               = subject.meg;
    mycfg.MEG         = MEG;
    niftis            = subject.mri;
    mycfg.MRI         = subject.mri;
    subject = subject.id;
end
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
    '/cubric/newscratch/314_wand/');
unix(['mkdir -p ' mycfg.SaveSubDir]);


% since the visual has already been run and has the coreigstered MRIs, look
% in this subjects' visual folder for the MRI
%--------------------------------------------------------------------------
coreg = dir([strrep(mycfg.SaveSubDir,'auditorymotor','visual'), '/*CoReg.mat']);
if ~isempty(coreg)
    fprintf('Fetching CoReg MRI from Visual analysis directory\n');
    cpstr = ['cp ' coreg(1).folder '/' coreg(1).name ' ' mycfg.SaveSubDir];
    unix(cpstr); pause(1);
end


% You can use this code to just compute the leadfields and then stop 
% (saving everything into the preproc dir for later use), or you
% can set ComputeBeamformer = 1 and specify your data / contrast parameters
% below.


% Optional: Only needed if mycfg.ComputeBeamformer = 1 ... 
%==========================================================================


% Define stimulus & frequency parameters: DESYNC
%--------------------------------------------------------------------------
mycfg.cov_fwin = [15 30];
mycfg.cov_toi  = [-1.25 1.75];

mycfg.trigger  = 'emg';  % trig/stim 
mycfg.b_toi    = [-1.25 -0.5];   % baseline times
mycfg.a_toi    = [-0.25 0.5];  % active times


% the manual local part - coregistration
if isempty(dir([mycfg.SaveSubDir '*CoReg.mat']))
    FieldtripLCMV_COREG(mycfg);
end

% save the input struct to load ont cluster
t = clock;
unique_name = [date '_' num2str(t(5)) num2str(round(t(6)))];
unique_name = [mycfg.SaveSubDir unique_name];

mycfg.prepend = 'MRBD_';

save(unique_name,'mycfg')

% Compute the lead fields and beamformer on the cluster

cd(mycfg.SaveSubDir);
if runoncluster
    fprintf('Submitting job to cluster...\n');
    docluster_slurm_bfm('FieldtripLCMVBeamformer',unique_name);
else
    feval('FieldtripLCMVBeamformer',unique_name);
end

% Define stimulus & frequency parameters: REBUUND
%--------------------------------------------------------------------------
mycfg.cov_fwin = [15 30];
mycfg.cov_toi  = [-1.25 1.75];

mycfg.trigger  = 'emg';  % trig/stim 
mycfg.b_toi    = [-1.25 -0.5];   % baseline times
mycfg.a_toi    = [1 1.75];  % active times

% the manual local part - coregistration
%if isempty(dir([mycfg.SaveSubDir '*CoReg.mat']))
%    FieldtripLCMV_COREG(mycfg);
%end

% save the input struct to load ont cluster
t = clock;
unique_name = [date '_' num2str(t(5)) num2str(round(t(6)))];
unique_name = [mycfg.SaveSubDir unique_name];

mycfg.prepend = 'PMBR_';

save(unique_name,'mycfg')

% Compute the lead fields and beamformer on the cluster

cd(mycfg.SaveSubDir);
if runoncluster
    fprintf('Submitting job to cluster...\n');
    docluster_slurm_bfm('FieldtripLCMVBeamformer',unique_name);
else
    feval('FieldtripLCMVBeamformer',unique_name);
end
% or run locally:
% FieldtripLCMVBeamformer(unique_name)
