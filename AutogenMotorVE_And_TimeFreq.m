
cd ~/WAND_MOTOR/

response = 'MRBD';

% update list of subjects with completed beamformer weights
str = ['ls /cubric/newscratch/314_wand/314_wand/bids/sourcedata/*/meg/auditorymotor/preproc/',...
    response '_CommonWeights.mat > ' response '_List.txt'];
unix(str);
list = ReadverifyDatasets([response '_List.txt']);

for i = 1:length(list)
    
    [fp,fn,fe] = fileparts(list{i});
    
    % extract sub ID
    subid = strrep(fp,'/cubric/newscratch/314_wand/314_wand/bids/sourcedata/','');
    subid = strrep(subid,'/meg/auditorymotor/preproc','');
    
    % move into subdir to make virtual electrode and run timefreq
    cd(fp);
    
    % Create virtual sensor if necessary
    if ~exist(['ftVE_WAND_' response '_L_' subid '.mat']) || ...
            ~exist(['ftVE_WAND_' response '_R_' subid '.mat'])
        ftve = FindMotorPeakMakeVS(fp);
        
        switch response
            case 'MRBD' ; ftve = ftve.mrbd;
            case 'PMBR' ; ftve = ftve.pmbr;
        end
        
    else
        f1 = load(['ftVE_WAND_' response '_L_' subid '.mat']);
        f2 = load(['ftVE_WAND_' response '_R_' subid '.mat']);
        
        % load the 2 fieldtrip structures and stack them as if they were
        % returned by FindMotorPeaksMakeVS
        ftve    = [];
        ftve{1} = f1.ftve;
        ftve{2} = f2.ftve;
    end
    
    %p.coord(i,:) = ftve.loc_template;
        
    switch response
        case 'MRBD' ; foi = [13 30]; woi = [-0.3 0.3];
        case 'PMBR' ; foi = [13 30]; woi = [1 1.75];
    end
    
    %if isfield(ftve,'ftve')
    %    ftve = ftve.ftve;
    %end
    
    % run timefrequency analysis & save bertogram
    name = [response '_L_' subid '_bert'];
    [pksL] = MotorVE_TimeFreq(ftve{1},foi,woi,1,name,response);
    name = [response '_R_' subid '_bert'];
    [pksR] = MotorVE_TimeFreq(ftve{2},foi,woi,1,name,response);
   
    
    close all; drawnow;
    
    p = struct;
    p.id = subid;
    p.MeanBetaL = pksL(1).MeanAmp;
    p.MeanBetaR = pksR(1).MeanAmp;
    
    PeakData(i) = p;
    
end
    cd ~/WAND_MOTOR/
save([response '_PeakData'],'PeakData');