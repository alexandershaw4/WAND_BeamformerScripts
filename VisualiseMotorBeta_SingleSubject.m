

mrbd = load('MRBD_SourceContrast.mat');
pmbr = load('PMBR_SourceContrast.mat');

i = mrbd.inside;
p = mrbd.pos_template;

mrbdp = mrbd.pow; mrbdp(mrbdp>0)=0;
pmbrp = pmbr.pow; pmbrp(pmbrp<0)=0;

lim1 = max(abs(mrbdp));
lim2 = max(abs(pmbrp));

figure('position',[566 410 1006 387]);
subplot(121); scatter3(p(i,1),p(i,2),p(i,3),100,mrbdp(i),'filled');
caxis([-lim1 lim1]);colorbar;view(0,90);
title('Beta Desync');

subplot(122); scatter3(p(i,1),p(i,2),p(i,3),100,pmbrp(i),'filled');
caxis([-lim2 lim2]);colorbar;view(0,90);
title('Beta Rebound');

% REDNER ON A BRAIN USING ATEMPLATE
addpath(genpath('~/code/SourceMesh/'))
figure('position',[252 290 1434 490]);

x = load('DenseAAL.mat')

s(1) = subplot(121);
D0=atemplate('mesh','def2','sourcemodel',p(i,:),'overlay',mrbdp(i),'fighnd',s(1),'post_parcel',{x.v,x.vi});

s(2) = subplot(122);
D1=atemplate('mesh','def2','sourcemodel',p(i,:),'overlay',pmbrp(i),'fighnd',s(2),'post_parcel',{x.v,x.vi});


figure('position',[252 290 1434 490]);
s(1) = subplot(121);
atemplate('overlay',D0.post_parcel.ParVal,'sourcemodel',D0.post_parcel.pos,'fighnd',s(1));
s(2) = subplot(122);
atemplate('overlay',D1.post_parcel.ParVal,'sourcemodel',D1.post_parcel.pos,'fighnd',s(2));




