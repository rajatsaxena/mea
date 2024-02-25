%% Load the rez file from Kilosort to extract the correct channel numbers

dir='X:\SWIL-Exp-Rajat\Spikesorted-SWIL'; 
[file,path]=uigetfile('*.mat','Select a rez mat file',dir);
load(strcat(path, file)); % this will take a while bc the file is huge

chnum=rez.ops.chanMap-1; % need to -1 so that channel numbers index starts from 0 (matches Sotiris')
shanknum=rez.ops.kcoords; % gives you the shank numbers I think for each channel

save(fullfile(path,'channel_numbers'),'chnum'); % these are the correct channel numbers that match to channel_positions.npy (DO NOT USE channel_map.npy)
save(fullfile(path,'shank_num'),'shanknum');
