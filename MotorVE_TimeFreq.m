function [pks,p] = MotorVE_TimeFreq(ftve,foi,woi,saveim,name,response)
% run time freuency analysis on the virtual sensor and produce a bertogram,
% return peak freauency and amplitude in frequency winfoe foi (e.g. [30
% 80]) and time window woi (e.g. [.300 .800] ms);

if nargin < 5 || isempty(response)
    response = 'MRBD';
end

if nargin < 4 || isempty(saveim)
    saveim=0;
end

if nargin < 3 || isempty(woi)
    woi = [-0.25 0.5; 1 1.75];
end

if nargin < 2 || isempty(foi)
    foi = [15 30];
end

% load the virtual sensor and run a timefrequency analysis
cfg.baseline = 'relchange';
cfg.sampletimes = ftve.time{1};
cfg.fsample = ftve.fsample;
cfg.filterorder = 3;
FoI = 1:.125:50;

MatDat = cat(1,ftve.trial{:});

tf = bert_singlechannel(MatDat,cfg,FoI,[-1.25 -0.5]);

agram = tf.agram;

% suppress effects in the opposite direction of interest
switch upper(response)
    case 'MRBD'; agram(agram>0)=0;
    case 'PMBR'; agram(agram<0)=0;
end

figure;p=plotbert(ftve.time{1},tf.freqs,agram);

if saveim
    export_fig(name,'-dpng','-m2','-nocrop');
    savefig(name);
end

m = agram;

for i = 1:size(woi,1)
    pks(i).woi = woi(i,:);
    pks(i).foi = foi;
    
    it = [findthenearest(woi(i,1),ftve.time{1}):findthenearest(woi(i,2),ftve.time{1})];
    
    iw = [findthenearest(foi(1),FoI):findthenearest(foi(2),FoI)];
    
    window = mean( m(iw,it),2 );
    
    %[~,I] = mean(window);
    %I=I(1);
    
    %pks(i).PeakFreq = FoI(iw(I));
    pks(i).MeanAmp = mean(window);
    
end
    