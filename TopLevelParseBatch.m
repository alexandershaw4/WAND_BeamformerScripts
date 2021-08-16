
addpath('~/code/WAND_BeamformerScripts/');

force = 0;

task = 'visual'; % {'visual' or 'auditorymotor'}

switch task 
    case 'auditorymotor'; ID = checkWANDmotor;
    case 'visual'; ID = checkWANDvisual;
end

notrun = [];
missing_meg = [];
missing_mri = [];
missing_hs  = [];

for i = 1:length(ID)
    
    this = ID(i);
    
    fprintf('\nRunning Dataset: %d/%d\n', i, length(ID) );
    fprintf('ID: %s\n',this.id);
    
    if this.hasmeg && this.hasmri 
        fprintf('Has MRI, & MEG \n');
        
        local = strrep(this.path,'/cubric/collab/',...
        '/cubric/newscratch/314_wand/');
    
        switch task
            case 'visual'; file = 'SourceContrast';
            case 'auditorymotor'; file = 'PMBR_CommonWeights';
        end
        
        
    if exist(local) && exist([local '/meg/' lower(task) '/preproc/' file '.mat']) && ~force
        fprintf('ALREADY RUN!/n');
        hasrun{i} = this.id;
        continue;
    %elseif i == 33 || i == 71
    %    continue
    else
        if contains(this.mri,'newscratch')
            switch task
                case 'visual'; WAND_config(this);
                case 'auditorymotor' ; WAND_config_motor(this);
            end
        else
            switch task
                case 'visual'; WAND_config(this.id);
                case 'auditorymotor';  WAND_config_motor(this.id);
            end
        end
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