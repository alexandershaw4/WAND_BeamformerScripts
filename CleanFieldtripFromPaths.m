% Clean out all fieldtrips from paths

phrases = {'ft_defaults' 'ft_read_mri' 'spm' 'ft_apply_montage'};

for j = 1:length(phrases)
    
    these = which(phrases{j},'-all');
    for i = 1:length(these)
        p = fileparts(these{i});
        p
        rmpath(genpath(p));
    end
    
end

CUB = dir('/cubric/software/MEG/fieldtrip*');

for i = 1:length(CUB)
    try
        pth = [CUB(i).folder '/' CUB(i).name];
        rmpath(genpath(pth));
    end
end

clear CUB j i p phrases pth these