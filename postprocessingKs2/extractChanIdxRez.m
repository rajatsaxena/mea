% basepath
basepath = 'X:\SWIL-Exp-Rajat\Spikesorted-SWIL';
% files
fnames = {'SWIL105PPC', 'SWIL11PPC', 'SWIL12PPC', 'SWIL13PPC', 'SWIL15PPC',  ...
    'SWIL18PPC', 'SWIL19PPC', 'SWIL20PPC', 'SWIL22PPC', 'SWIL23PPC', ...
    'SWIL24PPC', 'SWIL25PPC', 'SWIL26PPC', 'SWIL105VC', 'SWIL11VC', ... 
    'SWIL12VC', 'SWIL13VC', 'SWIL15VC', 'SWIL18VC', 'SWIL19VC', ... 
    'SWIL20VC', 'SWIL22VC', 'SWIL23VC', 'SWIL24VC', 'SWIL25VC', 'SWIL26VC'};

% iterate throug filesname
for i=1:length(fnames)
    fname = fnames{i};
    disp('Processing....');
    disp(fname);
    dat = load(fullfile(basepath, fname, 'rez.mat'));
    xcoords = dat.rez.xcoords;
    ycoords = dat.rez.ycoords;
    clear dat
    opfname = fullfile(basepath, fname, 'proc-channelmap.mat');
    save(opfname, 'xcoords', 'ycoords');
end