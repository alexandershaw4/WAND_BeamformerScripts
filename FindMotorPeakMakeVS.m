function ftve = FindMotorPeakMakeVS(subdir)

cd(subdir)

% Beta Desync Frist
load MRBD_SourceContrast

% coordinates of subject sourcemodel and AAL atlas model
%-------------------------------------------------------------------
p = pos_template;
i = inside;
x = load('DenseAAL.mat');
v = x.v;
b = [min(p(i,:));max(p(i,:))];
for i = 1:3
    v(:,i) = rescale(v(:,i),b(1,i),b(2,i));
end

% Comute euclidean distance between every point in meshes
%-------------------------------------------------------------------
D = cdist(p,v);
Dm= D*0;
for i = 1:size(Dm,2)
    [~,I]=min(D(:,i));
    Dm(I,i) = 1;
end
clear D

% Generate a mask using AAL parcel indices
%-------------------------------------------------------------------
%roi = [1 2 57 58 19 20 69 70];% SupMot,Precentral,Postcentral,ParacentralLob
roiL = [1 57 19 69];
roiR = [2 58 20 70];


idL = x.vi*0;
idR = x.vi*0;
for i = 1:length(roiL)
    idL(find(x.vi==roiL(i)))=1;
    idR(find(x.vi==roiR(i)))=1;
end

% Proejct the mask into subject space & get peak
%-------------------------------------------------------------------
maskL = ~~(Dm*idL);
maskL(find(~inside))=0;

maskR = ~~(Dm*idR);
maskR(find(~inside))=0;

pow(isnan(pow))=0;
[~,IL] = max(abs(pow(maskL)));
theseL = find(maskL);
peakindexL = theseL;%(IL);
loc_natL = pos(peakindexL,:);
loc_templateL = pos_template(peakindexL,:);

[~,IR] = max(abs(pow(maskR)));
theseR = find(maskR);
peakindexR = theseR;%(IR);
loc_natR = pos(peakindexR,:);
loc_templateR = pos_template(peakindexR,:);

save('MeshDistanceMat','Dm');

% scatter3(p(inside,1),p(inside,2),p(inside,3),90,pp(inside),'filled');hold on
% colormap(cmocean('balance'));colorbar;
% caxis([-50 50])
% alpha .4;
% view(0,90)


% % proof
% figure;
% s=scatter3(p(i,1),p(i,2),p(i,3),40,'b','filled');hold on;
% s2=scatter3(p(peakindex,1),p(peakindex,2),p(peakindex,3),180,'r','filled');hold on


% now we have peak coordinate we can load data and generate virtual sensor
%-------------------------------------------------------------------
optfile = dir('*2021*.mat'); load(optfile(1).name,'mycfg')

if ~strcmp(mycfg.prepend,'MRBD_')
    % wrong one loaded... theres only 2
    load(optfile(2).name,'mycfg');
end

% make active time the whole window of interest for the virtualelectrode
mycfg.a_toi = [-1.25 2.5];

% read, preprocess etc with same parameters as when generating beamformer
[tlck_actv,~] = ReloadDataOptionsForVS(mycfg);

cd(mycfg.SaveSubDir)
load([mycfg.prepend 'CommonwightsMat.mat'])

% generate virtual electrelectrode by passing data throguh weights for each
% trial

method = 'peak'; % peak or pcomponent

I = peakindexL;
% LEFT
for i = 1:size(tlck_actv.trial,1)
    switch method
        case 'pcomponent'
            vertseries = cat(1,wts{I})*squeeze(tlck_actv.trial(i,:,:));
            [u,s] = spm_svd(vertseries);
            s = diag(s);
            npc = find(cumsum(s)./sum(s) >= 0.9);
            npc = npc(1);

            if npc > 1
                trial{i} = sum(u(:,npc)'*vertseries,1);
            else
                trial{i} = u(:,1)'*vertseries;
            end
        case 'peak'
            peakindexL = theseL(IL);
            I = peakindexL;
            trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));         
            loc_natL = pos(peakindexL,:);
            loc_templateL = pos_template(peakindexL,:);
    end
    time{i} = tlck_actv.time;    
end

ftveL.trial = trial;
ftveL.time = time;

clear trial time

I = peakindexR;
% RIGHT
for i = 1:size(tlck_actv.trial,1)
        switch method
        case 'pcomponent'
            vertseries = cat(1,wts{I})*squeeze(tlck_actv.trial(i,:,:));
            [u,s] = spm_svd(vertseries);
            s = diag(s);
            npc = find(cumsum(s)./sum(s) >= 0.9);
            npc = npc(1);

            if npc > 1
                trial{i} = sum(u(:,npc)'*vertseries,1);
            else
                trial{i} = u(:,1)'*vertseries;
            end
        case 'peak'
            peakindexR = theseR(IR);
            I = peakindexR;
            trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));
            loc_natR = pos(peakindexR,:);
            loc_templateR = pos_template(peakindexR,:);
        end
    time{i} = tlck_actv.time;    
end

ftveR.trial = trial;
ftveR.time = time;

hdr = ft_read_header([mycfg.subject_dir mycfg.MEG]);

ftveL.fsample = hdr.Fs;
ftveL.label = 'MotorL_MRBD';
ftveR.fsample = hdr.Fs;
ftveR.label = 'MotorR_MRBD';



% generate name with ID for VE
[fp,id]=fileparts(mycfg.subject_dir);

nameL = ['ftVE_WAND_MRBD_L_' id];
nameR = ['ftVE_WAND_MRBD_R_' id];

ftve = ftveL;
save(nameL,'ftve','loc_natL','loc_templateL');
ftve = ftveR;
save(nameR,'ftve','loc_natR','loc_templateR');

ftve_mrbd = {ftveL ftveR};

clearvars -except subdir ftve_mrbd loc_natL loc_templateL loc_natR loc_templateR method

% Now on to PMBR response - same again........
%==============================================================
cd(subdir)

% Beta Desync Frist
load PMBR_SourceContrast

% coordinates of subject sourcemodel and AAL atlas model
%-------------------------------------------------------------------
p = pos_template;
i = inside;
x = load('DenseAAL.mat');
v = x.v;
b = [min(p(i,:));max(p(i,:))];
for i = 1:3
    v(:,i) = rescale(v(:,i),b(1,i),b(2,i));
end

% Comute euclidean distance between every point in meshes
%-------------------------------------------------------------------
D = cdist(p,v);
Dm= D*0;
for i = 1:size(Dm,2)
    [~,I]=min(D(:,i));
    Dm(I,i) = 1;
end
clear D % free up large matrix memory

% Generate a mask using AAL parcel indices
%-------------------------------------------------------------------
%roi = [1 2 60 61 57 58]; % L/R sup motor, precentral and postcentral
roiL = [1 57 19 69];
roiR = [2 58 20 70];


idL = x.vi*0;
idR = x.vi*0;
for i = 1:length(roiL)
    idL(find(x.vi==roiL(i)))=1;
    idR(find(x.vi==roiR(i)))=1;
end

% Proejct the mask into subject space & get peak
%-------------------------------------------------------------------
maskL = ~~(Dm*idL);
maskL(find(~inside))=0;

maskR = ~~(Dm*idR);
maskR(find(~inside))=0;


[~,IL] = max(abs(pow(maskL)));
theseL = find(maskL);
peakindexL = theseL(IL);
loc_natL = pos(peakindexL,:);
loc_templateL = pos_template(peakindexL,:);

[~,IR] = max(abs(pow(maskR)));
theseR = find(maskR);
peakindexR = theseR(IR);
loc_natR = pos(peakindexR,:);
loc_templateR = pos_template(peakindexR,:);


% now we have peak coordinate we can load data and generate virtual sensor
%-------------------------------------------------------------------
optfile = dir('*2021*.mat'); load(optfile(2).name,'mycfg')

if ~strcmp(mycfg.prepend,'PMBR_')
    % wrong one loaded... theres only 2
    load(optfile(1).name,'mycfg');
end

% make active time the whole window of interest for the virtualelectrode
mycfg.a_toi = [-1.25 2.5];

% read, preprocess etc with same parameters as when generating beamformer
[tlck_actv,~] = ReloadDataOptionsForVS(mycfg);

cd(mycfg.SaveSubDir)
load([mycfg.prepend 'CommonwightsMat.mat'])

% generate virtual electrelectrode by passing data throguh weights for each
% trial

I = peakindexL;
% LEFT
for i = 1:size(tlck_actv.trial,1)
        switch method
        case 'pcomponent'
            vertseries = cat(1,wts{I})*squeeze(tlck_actv.trial(i,:,:));
            [u,s] = spm_svd(vertseries);
            s = diag(s);
            npc = find(cumsum(s)./sum(s) >= 0.9);
            npc = npc(1);

            if npc > 1
                trial{i} = sum(u(:,npc)'*vertseries,1);
            else
                trial{i} = u(:,1)'*vertseries;
            end
        case 'peak'
            peakindexL = theseL(IL);
            I = peakindexL;
            trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));
            loc_natL = pos(peakindexL,:);
            loc_templateL = pos_template(peakindexL,:);
        end
    time{i} = tlck_actv.time;    
end

ftveL.trial = trial;
ftveL.time = time;

clear trial time

I = peakindexR;
% RIGHT
for i = 1:size(tlck_actv.trial,1)
        switch method
        case 'pcomponent'
            vertseries = cat(1,wts{I})*squeeze(tlck_actv.trial(i,:,:));
            [u,s] = spm_svd(vertseries);
            s = diag(s);
            npc = find(cumsum(s)./sum(s) >= 0.9);
            npc = npc(1);

            if npc > 1
                trial{i} = sum(u(:,npc)'*vertseries,1);
            else
                trial{i} = u(:,1)'*vertseries;
            end
        case 'peak'
            peakindexR = theseR(IR);
            I = peakindexR;
            trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));
            loc_natR = pos(peakindexR,:);
            loc_templateR = pos_template(peakindexR,:);
        end
    time{i} = tlck_actv.time;    
end

ftveR.trial = trial;
ftveR.time = time;

hdr = ft_read_header([mycfg.subject_dir mycfg.MEG]);

ftveL.fsample = hdr.Fs;
ftveL.label = 'MotorL_PMBR';
ftveR.fsample = hdr.Fs;
ftveR.label = 'MotorR_PMBR';

% generate name with ID for VE
[fp,id]=fileparts(mycfg.subject_dir);

nameL = ['ftVE_WAND_PMBR_L_' id];
nameR = ['ftVE_WAND_PMBR_R_' id];

ftve = ftveL;
save(nameL,'ftve','loc_natL','loc_templateL');
ftve = ftveR;
save(nameR,'ftve','loc_natR','loc_templateR');

ftve_pmbr = {ftveL ftveR};

%ftve = {ftve_mrbd ftve_pmbr};

clear ftve
ftve.mrbd = ftve_mrbd;
ftve.pmbr = ftve_pmbr;

end










