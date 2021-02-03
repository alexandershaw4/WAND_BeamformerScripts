function pks = VisualVE_TimeFreq(ftve,foi,woi,saveim,name)
% run time freuency analysis on the virtual sensor and produce a bertogram,
% return peak freauency and amplitude in frequency winfoe foi (e.g. [30
% 80]) and time window woi (e.g. [300 800] ms);

if nargin < 4 || isempty(saveim)
    saveim=0;
end

if nargin < 3 || isempty(woi)
    woi = [0 300; 300 1000];
end

if nargin < 2 || isempty(foi)
    foi = [30 80];
end

% load the virtual sensor and run a timefrequency analysis
cfg.baseline = 'relchange';

FoI = 3:.25:100;

tf = bert_singlechannel(ftve,cfg,FoI,[-1 0]);

figure;plotbert(ftve.time{1},tf.freqs,tf.agram,100);

if saveim
    export_fig(name,'-dpng','-m2','-nocrop');
end

m = tf.agram;

for i = 1:size(woi,1)
    pks(i).woi = woi(i,:);
    pks(i).foi = foi;
    
    it = [findthenearest(woi(i,1),ftve.time{1}):findthenearest(woi(i,2),ftve.time{1})];
    
    iw = [findthenearest(foi(1),FoI):findthenearest(foi(2),FoI)];
    
    window = mean( m(iw,it),2 );
    
    [~,I] = max(window);
    I=I(1);
    
    pks(i).PeakFreq = FoI(iw(I));
    pks(i).PeakAmp = window(I);
    
end
    

