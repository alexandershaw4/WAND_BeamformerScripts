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
roi = [1 2 60 61 57 58]; % L/R sup motor, precentral and postcentral
id = x.vi*0;
for i = 1:length(roi)
    id(find(x.vi==roi(i)))=1;
end

% Proejct the mask into subject space & get peak
%-------------------------------------------------------------------
mask = ~~(Dm*id);
[~,I] = max(abs(pow(mask)));
these = find(mask);
peakindex = these(I);
loc_nat = pos(peakindex,:);
loc_template = pos_template(peakindex,:);

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
I = peakindex;
for i = 1:size(tlck_actv.trial,1)
    trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));
    time{i} = tlck_actv.time;
end

ftve.trial = trial;
ftve.time = time;

hdr = ft_read_header([mycfg.subject_dir mycfg.MEG]);
ftve.fsample = hdr.Fs;
ftve.label = 'Motor';

% generate name with ID for VE
[fp,id]=fileparts(mycfg.subject_dir);

name = ['ftVE_WAND_MRBD_' id];

save(name,'ftve','loc_nat','loc_template');

ftve_mrbd = ftve;

clearvars -except subdir ftve_mrbd loc_nat loc_template

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
roi = [1 2 60 61 57 58]; % L/R sup motor, precentral and postcentral
id = x.vi*0;
for i = 1:length(roi)
    id(find(x.vi==roi(i)))=1;
end

% Proejct the mask into subject space & get peak
%-------------------------------------------------------------------
mask = ~~(Dm*id);
[~,I] = max(abs(pow(mask)));
these = find(mask);
peakindex = these(I);
loc_nat = pos(peakindex,:);
loc_template = pos_template(peakindex,:);

% % proof
% figure;
% s=scatter3(p(i,1),p(i,2),p(i,3),40,'b','filled');hold on;
% s2=scatter3(p(peakindex,1),p(peakindex,2),p(peakindex,3),180,'r','filled');hold on


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
I = peakindex;
for i = 1:size(tlck_actv.trial,1)
    trial{i} = wts{I}*squeeze(tlck_actv.trial(i,:,:));
    time{i} = tlck_actv.time;
end

ftve.trial = trial;
ftve.time = time;

hdr = ft_read_header([mycfg.subject_dir mycfg.MEG]);
ftve.fsample = hdr.Fs;
ftve.label = 'Motor';

% generate name with ID for VE
[fp,id]=fileparts(mycfg.subject_dir);

name = ['ftVE_WAND_PMBR_' id];

save(name,'ftve','loc_nat','loc_template');

ftve_pmbr = ftve;

ftve = {ftve_mrbd ftve_pmbr};



end










