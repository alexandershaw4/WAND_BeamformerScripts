function data = AutoAnalysisMainVis(cfg)

% Parameters
%-------------------------------------------------------------------------
dataset  = [cfg.subjdir '/' cfg.MEGds];
MRI      = [ cfg.MRI  ];
Pol      = [ cfg.Pol  ];
trigname = cfg.trigname;
times    = cfg.times;
freq     = cfg.freq;
a_toi    = cfg.atoi;
b_toi    = cfg.btoi;

% WAND GAMMA WORKFLOW
%-------------------------------------------------------------------------

% (1.) epoch the dataset using CTF cmdline tools
% (2.) reject EOG and jump artifacts using fieldtrip (writes new ds)
% (3.) coreg, MR segment and beamformer

% check if epoched version already exists - 
[fp,fn,fe] = fileparts(dataset);
ndataset = [fp '/' fn 'Cut' fe];

if ~exist(ndataset)
    
    % Epoch
    %-------------------------------------------------------------------------
    aMEG.epoch.epocher({dataset},'stim_on',times);

    [fp,fn,fe] = fileparts(dataset); % updates dataset to cut version
    dataset = [fp '/' fn 'Cut' fe];
else
    fprintf('Epoched dataset already exists\n');
    dataset = ndataset;
end

[fp,fn,fe] = fileparts(dataset);
nMEGds      = [fp '/' fn 'Fix' fe];

if ~exist(nMEGds)

    % Artifact Reject
    %-------------------------------------------------------------------------
    aMEG.artifact.ftclean_ctf(dataset,trigname,times)

    % Fix segmentation fault
    system(['sh /home/sapas10/code/+aMEG/fixsve.sh "' dataset '"']);  

    [fp,fn,fe] = fileparts(dataset);     % updates dataset to cut version
    MEGds      = [fn 'Fix' fe];
else
    fprintf('Cleaned & SegFault-Fixed Dataset Exists\n');
    MEGds = [fn 'Fix' fe];
end

% Coregister (Polhemus), Segment, Leadfields, Beamformer
%-------------------------------------------------------------------------
%addpath('~/code/');aPaths.defaults

% mycfg.subject_dir = '/cubric/scratch/sapas10/WAND/Visual/314_87034/';
% mycfg.MEG         = '314_87034_Visual.ds';
% mycfg.MRI         = '314_87034.mri';
% mycfg.Headshape   = '314_87034.pos';

mycfg.subject_dir = cfg.subjdir;
mycfg.MEG         = MEGds;
mycfg.MRI         = MRI;
mycfg.Headshape   = Pol;

mycfg.cov_fwin = freq; %   [40 80];
mycfg.cov_toi  = times;%   [-1.5 1.5];

mycfg.trigger = trigname;% 'stim_on'; % trig/stim 
mycfg.b_toi = b_toi;%      [-1.2 0];    % baseline times
mycfg.a_toi = a_toi;%      [0.3 1.5];   % active times

mycfg.SaveSubDir = cfg.SaveSubDir; %'VisBeam';

mycfg.UseKrish1mm = 0;

% Compute the beamformer
data = aMEG.aLCMV.FieldtripLCMVBeamformer(mycfg);

% add local paths back in
addpath(genpath('~/spm12'))
addpath('~/code/');aPaths.defaults
addpath(genpath('~/code/SourceMesh'));
addpath(genpath('~/fieldtrip-20170509/'))

% write out the nifti vol
% vol=data.source_diff_vol.pow;
% ft_write_mri([cfg.subjdir '/VisGamPow.nii'],vol,'dataformat','nifti')

% % make source plot 
% load src_diff.mat
% 
% figure;
% p = pos(inside,:);
% c = pow(inside);
% p = [p(:,2), p(:,1) p(:,3)]; % swap x/y
% D = atemplate('mesh','def4','overlay',c,'sourcemodel',p)