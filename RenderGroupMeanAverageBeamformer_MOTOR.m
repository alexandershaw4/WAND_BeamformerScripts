
pth = '/cubric/newscratch/314_wand/314_wand/bids/sourcedata/*/meg/auditorymotor/preproc/';

mrbd = dir([pth 'MRBD_SourceContrast.mat']); mrbd = strcat({mrbd.folder},'/', {mrbd.name});
pmbr = dir([pth 'PMBR_SourceContrast.mat']); pmbr = strcat({pmbr.folder},'/', {pmbr.name});

for i = 1:length(mrbd)
    
    load(mrbd{i},'pow');
    mrbdp(i,:) = pow(:);
    
    load(pmbr{i},'pow');
    pmbrp(i,:) = pow(:);
    
    if i == 1
       load(mrbd{i},'inside');
       load(mrbd{i},'pos_template');
    end
end
    
i = inside;
p = pos_template;

mrbdp(mrbdp>0)=0;
pmbrp(pmbrp<0)=0;

addpath(genpath('~/spm12/'));
mrbdpm = spm_robust_average(mrbdp(:,i),1);
pmbrpm = spm_robust_average(pmbrp(:,i),1);

lim1 = max(abs(mrbdpm));
lim2 = max(abs(pmbrpm));

figure('position',[566 410 1006 387]);
subplot(121); scatter3(p(i,1),p(i,2),p(i,3),100,mrbdpm,'filled');
caxis([-lim1 lim1]);colorbar;view(0,90);
title('Beta Desync');

subplot(122); scatter3(p(i,1),p(i,2),p(i,3),100,pmbrpm,'filled');
caxis([-lim2 lim2]);colorbar;view(0,90);
title('Beta Rebound');

% REDNER ON A BRAIN USING ATEMPLATE
addpath(genpath('~/code/SourceMesh/'))
figure('position',[343 119 1348 649]);

x = load('DenseAAL.mat');

nmesh         = gifti('BrainMesh_Ch2.gii');
mesh          = [];
mesh.vertices = nmesh.vertices;
mesh.faces    = nmesh.faces;
%mesh = spm_mesh_inflate(mesh,40);


s(1) = subplot(121);
D0=atemplate('mesh',mesh,'sourcemodel',p(i,:),'overlay',mrbdpm,'fighnd',s(1),'post_parcel',{x.v,x.vi});
D0.overlay.cb.CRange = D0.overlay.cb.CRange*1;
D0.overlay.cb.CLim = D0.overlay.cb.CRange;

s(2) = subplot(122);
D1=atemplate('mesh',mesh,'sourcemodel',p(i,:),'overlay',pmbrpm,'fighnd',s(2),'post_parcel',{x.v,x.vi});


% figure('position',[252 290 1434 490]);
% s(1) = subplot(121);
% atemplate('overlay',D0.post_parcel.ParVal,'sourcemodel',D0.post_parcel.pos,'fighnd',s(1));
% s(2) = subplot(122);
% atemplate('overlay',D1.post_parcel.ParVal,'sourcemodel',D1.post_parcel.pos,'fighnd',s(2));

