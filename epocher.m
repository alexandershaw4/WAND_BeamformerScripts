function epocher(List,Marker,Times,outpath)
% Use after AddMarker_CTF.m
% to cut around the markers
%
% AS

if ~iscell(List)
    Datasets = ReadverifyDatasets(List);
else
    Datasets = List;
end

if nargin < 4 || isempty(outpath)
    outpath = [];
end

for i = 1:length(Datasets)
    fprintf('\nEpoching dataset %d (%s)\n',i,Datasets{i});
    [p,n,e] = fileparts(Datasets{i});
    cd(p);
    epoch([p,'/',n,e],Marker,Times,outpath);
end

end

function epoch(Dataset,Marker,Times,outpath)

[PATHSTR,NAME,EXT] = fileparts(Dataset);

if ~isempty(outpath)
    DatasetOut = [outpath '/' NAME 'Cut.ds'];
else
    DatasetOut = [PATHSTR '/' NAME 'Cut.ds'] ;
end

cmdstring = ['newDs -f -marker ' Marker ' -time ' num2str(Times(1)) ' ' num2str(Times(2)) ' ' Dataset ' ' DatasetOut];

unix(cmdstring)

%Keep a reference for the user of where the original dataset came from just in case!
% fid = fopen(['UNCUT_DATASET'], 'w');
% fprintf(fid, '%s\n', Dataset);
% fclose(fid);
% unix(cmdstring);

end
