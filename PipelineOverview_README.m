
% WAND visual gamma analysis pipeline:
%---------------------------------------------------------
%
% STEP 1: Coregistration & creation of source & head models, lead fields 
% and global filters:
%
% 1) TopLevelParseBatch.m: Script that finds data and checks what's already run
%    - calls checkWANDvisual.m: reads collab folder and sorts subjects
%    - calls (2) for each subject
% 2) WAND_config.m function: sets up sub dir, moves MRI and lunches:
%   3) FieldtripLCMV_COREG.m (if not already run): interactive fieldtrip
%      coregiatration, upon closing/completetion launches:
%   4) FieldtripLCMVBeamformer.m: runs segmentation, creates headmodel,
%   source model, lead fields and global beamformer weights (filters)
%   [note (3) runs entirely on the cluster]
%
%
% STEP 2: Identify peak-voxel from %-change contrats (active vs baseline),
% create virtual sensor, run timefrequency aalysis and find peaks
%
% 5) AutogenerateVirtualSensorsandTimeFreq.m - loops subjects, calls:
%    6) FindVisualPeakMakeVS.m - finds peak, generates virtual sensor
%       7) ReloadDataOptionsForVS.m - reloads data, computes VS
%    8) VisualVE_TimeFreq.m - runs timefreq analysis, created bert plot and
%    extracts peak gamma frequency amplitude 
%
% 9) Render the contrast on a full brain: 
% Visualisecontrast_ReembedPointsAtemplate.m