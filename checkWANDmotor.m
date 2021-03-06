function [ID]=checkWANDmotor()
% determine which WAND datasets have MEG visual, MRI and headshape files
%
% returns a structure with subject fields with full paths to each data part
% also returns a matrix of logical for what people have
%

basepth = '/cubric/collab/314_wand/bids/sourcedata';
bidsmeg = '/meg/auditorymotor/raw/';
bidsmri = '/mri/anat/';
bidshed = '/meg/generic/';

cd(basepth);

IDs = dir('314*'); IDs = {IDs.name}';

for i = 1:length(IDs)
    %this = ['id_' IDs{i};];
    
    ID(i).id = IDs{i};
    ID(i).path = [basepth '/' IDs{i}];
    
    % Check for MEG
    meg = ( dir([basepth '/' IDs{i} bidsmeg '/*.ds']) );
    
    if length(meg) > 0
        ID(i).meg = [meg(1).folder '/' meg(1).name];
        ID(i).hasmeg = 1;
    else
        ID(i).hasmeg = 0;
    end
    
    % Check for MRI
    mri = ( dir([basepth '/' IDs{i} bidsmri '/*.nii']) );
    
    if length(mri) > 0
        ID(i).mri = [mri(1).folder '/' mri(1).name];
        ID(i).hasmri = 1;
    else
        ID(i).hasmri = 0;
    end
    
    % double check  - sometimes i have copied the mri to the newscratch
    % for those who have meg and i found MRs
    if ID(i).hasmeg && ~ID(i).hasmri
        fprintf('checking newscratch for mri...\n')
        
        ns = ['/cubric/newscratch/314_wand/314_wand/bids/sourcedata/'];
        mri = ( dir([ns IDs{i} '/*.nii']) );
        
        if length(mri) > 0
            ID(i).mri = [mri(1).folder '/' mri(1).name];
            ID(i).hasmri = 1;
            fprintf('Found MRI for %s in newscratch!\n',IDs{i});
        end
        
    end
    
    
    % Check for headshape
    headshape = ( dir([basepth '/' IDs{i} bidshed '/*.pos']) );
    
    if length(headshape) > 0
        ID(i).head = [headshape(1).folder '/' headshape(1).name];
        ID(i).hashead = 1;
    else
        ID(i).hashead = 0;
    end    
    
end