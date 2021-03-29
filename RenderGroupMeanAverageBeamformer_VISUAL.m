cd ~/WAND_VISUAL/;

pth = '/cubric/newscratch/314_wand/314_wand/bids/sourcedata/*/meg/visual/preproc/';

vis = dir([pth 'SourceContrast.mat']); vis = strcat({vis.folder},'/', {vis.name});


for i = 1:length(vis)
        
    vp = load(vis{i},'pow');
    pow(i,:) = vp.pow(:);
    
    if i == 1
       load(vis{i},'inside');
       load(vis{i},'pos_template');
    end
end

i = inside;
p = pos_template;

% construct whole head 1mm grid to ovelray occipital regions on
addpath('/cubric/scratch/krish/NewAAL/');
load AAL_labels.mat;
AtlasLabels = labels;
template=kMakeCustomTemplateSourceGrid([-7.5 7.5],[-11.4 7.8],[-7.4 8.8],0.1,'/cubric/software/MEG/fieldtrip-20161011/template/atlas/aal/ROI_MNI_V4.nii',AtlasLabels);
sm = ft_convert_units(template, 'mm'); 

p0=sm.pos;
i0 = find(sm.inside);


% generate robust mean
pow(pow<0)=0;
addpath(genpath('~/spm12/'));
powmean = spm_robust_average(pow(:,i),1);


% stack into new matrices and render with atemplate
source = [p0(i0,:);p(i,:)];
overlay = [p0(i0,1)*0;powmean(:)];
addpath(genpath('~/code/SourceMesh/'))

%afigure;atemplate('sourcemodel',source,'overlay',overlay,'open');

% REDNER ON A BRAIN USING ATEMPLATE
figure('position',[343 119 1348 649]);

D0=atemplate('mesh','def2','sourcemodel',source,'overlay',overlay,'open');

savefig('GroupMeanVisualOverlay.fig');

export_fig('GroupMeanVisualOverlay','-dpng','-m4','-transparent')




