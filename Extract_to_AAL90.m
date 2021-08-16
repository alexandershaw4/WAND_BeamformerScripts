function Extract_to_AAL90(pth)
% based on the data structure of the WAND data in newscratch (bids), find
% the correct files to reconstruct virtual electrodes/sensors for the 90
% AAL regions.
%
% rerquires the sourcemodel (coords), beamformer weights, data (saved in
% 'TheData.mat')
%
% give it the path to a subjects preproc bids dir, e.g.
% '/cubric/newscratch/314_wand/314_wand/bids/sourcedata/314_05117/meg/visual/preproc'

% extract sub id
id = strrep(pth,'/meg/visual/preproc','');
[~,id]=fileparts(id);

cd(pth);
sm  = load('sourcemodel_indiv.mat');
wt  = load('CommonwightsMat.mat');
aal = load('AAL_90_SOURCEMOD.mat');

dat = load('PreprocData.mat');

i = sm.inside;
wtm = wt.wts(i);

pos0 = aal.template_sourcemodel.pos;
pos1 = sm.pos(i,:);
    
% align the box boundaries
pos0 = pos0 - repmat(spherefit(pos0),[size(pos0,1) 1]);
pos1 = pos1 - repmat(spherefit(pos1),[size(pos1,1) 1]);

pos1  = pos1*10; % convert to mm
scale = mean([min(pos0); max(pos0)] ./ [min(pos1); max(pos1)]);
pos1 = pos1.*repmat(scale,[size(pos1,1),1]);

% make fake fieldtrip struct
virtsen = struct;
virtsen.fsample = dat.data_preproc.fsample;
virtsen.label = aal.AAL_Labels;
virtsen.pos = pos0;

for i = 1:90
    fprintf('Computing for atlas region %d/%d\n',i,90);
    [val,ind]=maxpoints(cdist(pos0(i,:),pos1),10,'min');
    theweights = cat(1,wtm{ind});
    
    for j = 1:length(dat.data_preproc.trial)
        virtsen.trial{j}(i,:) = mean( theweights*dat.data_preproc.trial{j} ,1);
        virtsen.time{j} = dat.data_preproc.time{j};
    end
end

name = ['AAL_Nodes_' id];
save(name,'virtsen');