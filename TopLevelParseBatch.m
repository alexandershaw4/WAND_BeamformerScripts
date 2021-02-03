
addpath('~/code/WAND_BeamformerScripts/');

ID = checkWANDvisual;

notrun = [];
missing_meg = [];
missing_mri = [];
missing_hs  = [];

for i = 2:2;%:length(ID)
    
    this = ID(i);
    
    fprintf('\nRunning Dataset: %d/%d\n', i, length(ID) );
    fprintf('ID: %s\n',this.id);
    
    if this.hasmeg && this.hasmri 
        fprintf('Has MRI, & MEG \n');
        
        local = strrep(this.path,'/cubric/collab/',...
        '/cubric/newscratch/314_wand/');
    
        file = 'CommonWeights';
        
    if exist(local) && exist([local '/' file])
        fprintf('ALREADY RUN!/n');
        continue;
    else
        
        WAND_config(this.id);
    end
        
    else
        
        % compile list to sort out afterward
        notrun = [notrun; this.id];
        
        if this.hasmeg == 0
            fprintf('Missing MEG\n');
            missing_meg = [missing_meg; this.id];
        end
        if this.hasmri == 0
            fprintf('Missing MRI\n');
            missing_mri = [missing_mri; this.id];
        end
        if this.hashead == 0
            fprintf('Missing headshape\n');
            missing_hs = [missing_hs; this.id];
        end
        
    end
end