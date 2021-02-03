function ftve = FindVisualPeakMakeVS(subdir)

cd(subdir)

load SourceContrast

[~,I]=max(pow);

loc_nat = pos(I,:);
loc_template = pos_template(I,:);

clearvars -except subdir I loc_nat loc_template

load CommonwightsMat.mat
% filter for voxel I: wts{i}*chan data

% reload the data with same options as when we created the beamformer

optfile = dir('*2021*.mat'); load(optfile.name,'mycfg')

% make active time the whole window of interest for the virtualelectrode
mycfg.a_toi = [-1 3];

% read, preprocess etc with same parameters as when generating beamformer
[tlck_actv,tlck_bsln] = ReloadDataOptionsForVS(mycfg);

cd(mycfg.SaveSubDir)

% generate virtual electrelectrode by passing data throguh weights for each
% trial
for i = 1:size(tlck_actv.trial,1)
    trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));
    time{i} = tlck_actv.time;
end

ftve.trial = trial;
ftve.time = time;

hdr = ft_read_header([mycfg.subject_dir mycfg.MEG]);
ftve.fsample = hdr.Fs;
ftve.label = 'V1';

% generate name with ID for VE
[fp,id]=fileparts(mycfg.subject_dir);

name = ['ftVE_WAND_VISUAL_' id];

save(name,'ftve','loc_nat','loc_template');

% export 
addpath(genpath('~/spm12'));
spm_eeg_VEft2spmLFP(ftve,['SPM_' name])

