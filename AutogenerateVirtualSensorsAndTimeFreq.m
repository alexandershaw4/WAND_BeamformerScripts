
cd ~/WAND_VISUAL/

% update list of subjects with completed beamformer weights
!ls /cubric/newscratch/314_wand/314_wand/bids/sourcedata/*/meg/visual/preproc/CommonwightsMat.mat > VisList.txt
list = ReadverifyDatasets('VisList.txt');

for i = 1:length(list)
    
    [fp,fn,fe] = fileparts(list{i});
    
    % extract sub ID
    subid = strrep(fp,'/cubric/newscratch/314_wand/314_wand/bids/sourcedata/','');
    subid = strrep(subid,'/meg/visual/preproc','');
    
    % move into subdir to make virtual electrode and run timefreq
    cd(fp);
    
    % Create virtual sensor if necessary
    if ~exist(['ftVE_WAND_VISUAL_' subid]);
        ftve = FindVisualPeakMakeVS(fp);
    else
        ftve = load(['ftVE_WAND_VISUAL_' subid]);
    end
        
    % run timefrequency analysis & save bertogram
    pks = VisualVE_TimeFreq(ftve,[],[],1,['bertogram_' subid]);
    
    p = struct;
    p.id = subid;
    p.spike_ga = pks(1).PeakAmp;
    p.spike_gf = pks(1).PeakFreq;
    p.sust_ga = pks(2).PeakAmp;
    p.sust_gf = pks(2).PeakFreq;
    
    PeakData(i) = p;
    
end
    cd ~/WAND_VISUAL/
save('PeakData','PeakData');
    