% Script to extract coordinates of VEs for the motor MRBD & PMBR

pth = '/cubric/newscratch/314_wand/314_wand/bids/sourcedata/*/meg/auditorymotor/preproc/';

mrbdl = dir([pth 'ftVE*MRBD_L*']);mrbdl = strcat({mrbdl.folder},'/', {mrbdl.name});
mrbdr = dir([pth 'ftVE*MRBD_R*']);mrbdr = strcat({mrbdr.folder},'/', {mrbdr.name});
pmbrl = dir([pth 'ftVE*PMBR_L*']);pmbrl = strcat({pmbrl.folder},'/', {pmbrl.name});
pmbrr = dir([pth 'ftVE*PMBR_R*']);pmbrr = strcat({pmbrr.folder},'/', {pmbrr.name});

for i = 1:length(mrbdl)
    x = load(mrbdl{i});    
    mrbdl_coord(i,:) = x.loc_templateL;
    x = load(mrbdr{i});    
    mrbdr_coord(i,:) = x.loc_templateR;    
    
    x = load(pmbrl{i});    
    pmbrl_coord(i,:) = x.loc_templateL;
    x = load(pmbrr{i});    
    pmbrr_coord(i,:) = x.loc_templateR;    
end


s = load([(fileparts(mrbdl{1})) '/sourcemodel_template.mat']);
p = s.pos(s.inside,:);

%x = load('AAL_SuperMesh_Sourcemodel.mat');
%pos = x.pos;
%mesh.vertices = x.pos;
%mesh.faces = x.face;

% nmesh         = read_nv;
% mesh          = [];
% mesh.vertices = nmesh.vertices;
% mesh.faces    = nmesh.faces;
                nmesh         = gifti('BrainMesh_Ch2.gii');
                mesh          = [];
                mesh.vertices = nmesh.vertices;
                mesh.faces    = nmesh.faces;
                mesh = spm_mesh_inflate(mesh,40);
                C = docurvature(mesh);
pos = mesh.vertices;


scale=[min(pos);max(pos)]./[min(p);max(p)];
aff = diag(mean(scale,1));

figure('position',[642   132   603   664]);
pt = patch(mesh);
pt.FaceVertexCData = C;
pt.EdgeColor = 'none';
pt.FaceColor='interp';
pt.FaceAlpha=0.5;
colormap(1-gray);
hold on;

cent = spherefit(mesh.vertices);

% draw normal lines
p1 = mrbdl_coord*aff;
p2 = mrbdr_coord*aff;
p3 = pmbrl_coord*aff;
p4 = pmbrr_coord*aff;

% Comute the normal line for each point and push it out to the boundary of
% t
% p11 = (p1*2)-(cent); % normal
% p13 = (p1*3)-(cent*2);
% dt = mean(diff([p1; p11; p13]),1);
% %b = [cent; p1; p11; p13];
% b = ((1:0.125:2)'*dt) + cent;
% l = line(b(:,1),b(:,2),b(:,3));
% 
% D=cdist(b,mesh.vertices);
% for i = 1:size(D,1)
%     [val(i),I(i)] = maxpoints(D(i,:),1,'min');
% end
% [~,cl]=min(val);
% newpnt = mesh.vertices(I(cl),:);


scatter3(p1(:,1),p1(:,2),p1(:,3),30,'b','filled');
scatter3(p2(:,1),p2(:,2),p2(:,3),30,'b','filled');
scatter3(p3(:,1),p3(:,2),p3(:,3),30,'r','filled');
scatter3(p4(:,1),p4(:,2),p4(:,3),30,'r','filled');



