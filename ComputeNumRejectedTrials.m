function nt = ComputeNumRejectedTrials()

cd ~/WAND_VISUAL/ ;

% update list of subjects with completed beamformer weights
%!ls /cubric/newscratch/314_wand/314_wand/bids/sourcedata/*/meg/visual/preproc/CommonwightsMat.mat > VisList.txt
list = ReadverifyDatasets('VisList.txt');

for i = 1:length(list)
    
    fprintf('Reading & logging num bad...(%d/%d)\n',i,length(list));
    
    [fp,fn,fe] = fileparts(list{i});
    
    % extract sub ID
    subid = strrep(fp,'/cubric/newscratch/314_wand/314_wand/bids/sourcedata/','');
    subid = strrep(subid,'/meg/visual/preproc','');
    
    % move into subdir to make virtual electrode and run timefreq
    cd(fp);
    
    load(['ftVE_WAND_VISUAL_' subid])
    
    nt(i) = length(ftve.trial);
    id{i} = subid;
    
%     optfile = dir('*2021*.mat'); load(optfile.name,'mycfg')
% 
%     % make active time the whole window of interest for the virtualelectrode
%     mycfg.a_toi = [-1 3];
%     
%     mycfg.reject_only = 1; % just compute artifact rejection and return numbers
%     
%     % read, preprocess etc with same parameters as when generating beamformer
%      [~,~,nt(i,:)] = ReloadDataOptionsForVS(mycfg);
%      
%      id{i} = subid;
     
     
end
cd ~/WAND_VISUAL/
     save('VisNumReject','nt','id');