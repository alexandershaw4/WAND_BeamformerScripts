

% construct whole head 1mm grid to ovelray occipital regions on
addpath('/cubric/scratch/krish/NewAAL/');
load AAL_labels.mat;
AtlasLabels = labels;
template=kMakeCustomTemplateSourceGrid([-7.5 7.5],[-11.4 7.8],[-7.4 8.8],0.1,'/cubric/software/MEG/fieldtrip-20161011/template/atlas/aal/ROI_MNI_V4.nii',AtlasLabels);
sm = ft_convert_units(template, 'mm'); 

p0=sm.pos;
i0 = find(sm.inside);

%figure;
%scatter3(p0(i0,1),p0(i0,2),p0(i0,3),70,'MarkerFaceColor',[.4 .4 .4],'MarkerEdgeColor','none'); 
%alpha .4;hold on;

load('SourceContrast.mat')
i = inside;
p = pos_template;
%scatter3(p(i,1),p(i,2),p(i,3),100,pow(i),'filled')

% stack into new matrices and render with atemplate
source = [p0;p];
overlay = [p0(:,1)*0;pow];
addpath(genpath('~/code/SourceMesh/'))
afigure;atemplate('sourcemodel',source,'overlay',overlay,'open');