function ftclean_ctf(dataset,trigname,times)
% - use fieldtrip to artifact clean a CTF MEG dataset
% - input epoched CTF dataset
% - writes the bad trials into the CTF classfile
%   (running newDs on the resulting dataset with remove the bad trials)
%
% AS

%dataset = '/cubric/scratch/sapas10/WAND/Visual/314_04037/314_04037_Visual.ds';
%epocher({dataset},'stim_on',[-2 2])
%[fp,fn,fe] = fileparts(dataset);
%dataset = [fp '/' fn 'Cut' fe];

% %define trials
cfg = [];
cfg.dataset             = dataset;
cfg.trialdef.eventtype  = trigname;    %'stim_on';
cfg.trialdef.prestim    = abs(times(1)) ;   %2;  % fix me
cfg.trialdef.poststim   = abs(times(2)) ;   %2;
cfg_data = ft_definetrial(cfg);


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




cfg.artfctdef.reject          = 'none', 'partial', 'complete', 'nan', or 'value' (default = 'complete')
cfg.artfctdef.minaccepttim    = when using partial rejection, minimum length
in seconds of remaining trial (default = 0.1)
cfg.artfctdef.crittoilim      = when using complete rejection, reject trial only when artifacts occur within
this time window (default = whole trial). This only works with in-memory data,
since trial time axes are unknown for data on disk.
cfg.artfctdef.feedback        = 'yes' or 'no' (default = 'no')
cfg.artfctdef.invert          = 'yes' or 'no' (default = 'no')
cfg.artfctdef.value           = scalar value to replace the data in the artifact segments (default = nan)
cfg.artfctdef.eog.artifact    = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_EOG
cfg.artfctdef.jump.artifact   = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_JUMP
cfg.artfctdef.muscle.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_MUSCLE
cfg.artfctdef.zvalue.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_ZVALUE
cfg.artfctdef.visual.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using FT_DATABROWSER
cfg.artfctdef.xxx.artifact    = Nx2 matrix with artifact segments, this could be added by your own artifact detection function

[cfg] = ft_rejectartifact(cfg, data)



% % compile list of segments either with jumps or with EOG contam
% %----------------------------------------------------------------
% trl = cfg.artfctdef.zvalue.trl;
% art = [artifact_EOG ; artifact_jump ];
% bad = [];
% 
% % find which trials corrspond to which bad segments
% 
% for i = 1:size(art,1) % each bad chunk
%     x0  = findthenearest( art(i,1),trl(:,1) );
%     x1  = findthenearest( art(i,2),trl(:,2) );
%     tbad = min([x0 x1]):max([x0 x1]);
%     %tbad
%     bad = unique([bad(:); tbad(:) ]);
% end
% 
% 
% [newexcludelength] = ctf_write_BadTrials(bad, dataset);






% % reject trials
% cfg = [];
% cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
% cfg.artfctdef.eog.artifact = artifact_EOG; %
% %cfg.artfctdef.jump.artifact = artifact_jump;
% %cfg.artfctdef.muscle.artifact = artifact_muscle;
% data_no_artifacts = ft_rejectartifact(cfg,data);
