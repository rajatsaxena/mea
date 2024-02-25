% Example pipeline for processing histology
% taken from Andy Peters, modified by Winny Ning
% further modified by Rajat Saxena

% init variables
aname='SWIL22'; % animal name
region='PPC'; % brain region (PPC/ VC)
apfs=30000.0; %spike data sampling rate 
ksparentdir='T:\SWIL-Rajat\Spikesorted-SWIL';
goodchans = false; %set true to load good units only, otherwise load good/mua
lfpfs=1000; % lfp sampling rate
lfpparentdir = 'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL'; % lfp directory
lfpts = [1800 3600]; % start and end times (seconds) for lfp data
probe = 'intan_256F';
probemap = load(fullfile('.\ProbeFiles','chanMap_intan256F_bottom_SM.mat'));

%% 1) Load CCF and set paths for slide and slice images
addpath(genpath('npy-matlab-master'));
addpath(genpath('Allen_Histology')); 
addpath(genpath('brewermap'));
addpath(genpath('allenCCF_repo_functions'));

% Load CCF atlas
allen_atlas_path = '.\Allen_Histology'; % path to CCF atlas
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area
st = AP_loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % table of all labels

% Set paths for histology images and directory to save slice/alignment
dir=fullfile('T:\SWIL-Rajat\Histology-SWIL',strcat(aname,region)); 
[file,im_path]=uigetfile('*.tif','Select a histology image',dir); % im_path = path to folder with images goes here
addpath(im_path); % make image is in matlab path
slice_path = [im_path(1:end-1) filesep 'slices']; % later on will make slices folder in path

%% 2) Preprocess slide images to produce slice images

% Set resize factor
% resize_factor = []; % (slides ome.tiff: auto-resize ~CCF size 10um/px)
resize_factor = 1; % (slides tiff: resize factor)

% Set slide or slice images
slice_images = false; % (images are slides - extract individual slices)
slice_images = true; % (images are already individual slices)

% Preprocess images
AP_process_histology(im_path,resize_factor,slice_images);

% (optional) Rotate, center, pad, flip slice images 
% Draw reference line at midline to rotate/center image a bit
AP_rotate_histology(slice_path);

%% 3) Align CCF to slices

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path); % manually match AP coronal slice, can also match angles if slicing unevenly

% Align CCF slices and histology slices
AP_auto_align_histology_ccf(slice_path);

% (second: curate manually)
% manually add additional points (>3) on histology and CCF slices to make alignment better
AP_manual_align_histology_ccf(tv,av,st,slice_path);


%% 4) Utilize aligned CCF

% Display aligned CCF over histology slices
AP_view_aligned_histology(st,slice_path);

% Display histology within 3D CCF
AP_view_aligned_histology_volume(tv,av,st,slice_path,1);

% Get probe trajectory from histology, convert to CCF coordinates
% Enter num of probes + draw vertical lines on histology images
% Probe 1 will be medial, Probe 2 is lateral S1
AP_get_probe_histology(tv,av,st,slice_path);

%% 5) Align histology to electrophysiology data and save brain regions

% This will match trajectory of probe in CCF with ephys signatures.
% Requires kilosort output of spike_times, spike_templates and template_depths.
[spike_times, spike_templates, templates, channel_positions, template_depths, template_depths_flipped] = AP_loadKs2(ksparentdir, strcat(aname,region), apfs, goodchans);

% find missing channel numbers that were out of brain or dead
% rather than loading rez.mat, optimize this: RS
channelmap = AP_loadChannelNum(probemap, channel_positions);

% load lfp data
lfpallchans=double(readNPY(fullfile(lfpparentdir, 'SWILRound5', aname, strcat(aname,'-lfp.npy')))); % load in lfp data
lfpallchans=lfpallchans(:,lfpts(1)*lfpfs:lfpts(2)*lfpfs); % subset lfp
if strcmp(region,'PPC')
    lfp=lfpallchans(256:512,:);
else
    lfp=lfpallchans(1:256,:);
end

% repeat for each shank
use_probe = 1; % shank num
shmap = channelmap(channelmap(:,2)==use_probe,:);
lfp_sh = lfp(shmap(:,1),:); % lfp for each shank

% For each shank, sort channels according to DEPTH
% depth is not flipped yet
[~,ind]=sort(shmap(:,4),'descend'); % col 3 (Y) is DV depth, sort from superficial (2125, 1st row) to deep (0), same as kilosort (**DV will flip later)
shmap_sorted = shmap(ind,:);
lfp_sh=lfp_sh(ind,:); %%% need to sort this by depth index because lfp_sh was sorted by chnum

% flipped ycoords positions for each shanks
Ychannel_positions_sh = 2125-(shmap_sorted(:,4)); % flipped ycoords from either shank 1 or 2, this is already sorted Y position from superficial (0) to deep (2125)

% align probe histology with depth, lfp, mua
figname = fullfile(slice_path,strcat('slice',string(use_probe),'.fig'));
AP_align_probe_histology(st,slice_path, ...
    spike_times,spike_templates,template_depths_flipped, ...
    lfp_sh,Ychannel_positions_sh, ...
    use_probe, figname);


% Analyze the output mat files and get the laminar boundaries
% probe_ccf.trajectory_coords = 3D CCF coords of the probe trajectory
% probe_ccf.trajectory_areas = annotated CCF areas for each point
% probe_ccf.probe_depths = relative depth of the probe to that point (0 is
% superficial (headstage) and 2125 is tip of probe (deepest layers). Note
% this is reversed from Kilosort output
addpath(genpath(slice_path));
load('probe_ccf.mat');

% get regions for all channels
brainregions=st(probe_ccf(use_probe).trajectory_areas(:),:).safe_name; 
% find nearest region for each channel
idx = interp1(probe_ccf(use_probe).probe_depths',1:length(probe_ccf(use_probe).probe_depths),shmap(:,4),'nearest');
areaID{use_probe}=brainregions(idx);
% save output areas for each channel in npy file
if length(areaID)==2
    writeNPY(fullfile(slice_path,strcat(aname,'-brainRegions.npy')));
end

%% 6) Plot figure of CCF brain structures in 2D and 3D with probe tracts 
% (code from Muhang Li, Riken; modified by Winny Ning; Rajat)
%%% Get brain areas in search, then draw

% Prompt for which structures to show (only structures which are
% labelled in the slice-spacing downsampled annotated volume)
% structure_search = lower(inputdlg('Search structures'));
% structure_match = find(contains(lower(st.safe_name),structure_search)); % find structure name in st.safe_name
% 
% selected_structure = listdlg('PromptString','Select a structure to plot:', ...
%     'ListString',st.safe_name(structure_match),'ListSize',[520,500], ...
%     'SelectionMode','single');
% 
% plot_structure = structure_match(selected_structure);

% Posterior parietal association (
% Field CA1 (458); DG (474); PPC/anterior area (340/347); RSCag (299); AM
% (172); AL (165); V1 (186); Vis ctx (158); Visrl (354); lateral visual
% area (179); Laterointermediate (207), Auditory primary (137), TeA (361)

if strcmp(region,'PPC')
    structures = [458 474 347 299 172]; 
else
    structures = [458 474 165 172 179 186 354 347 207 137 361];
end
[ccf_3d_axes, ccf_axes] = AP_plot_CCF(tv,av,st,structures);
hold on;

%%% Draw probes (output from AP_histology)
probe_color = 'r';
AP_draw_probe(ccf_3d_axes, ccf_axes, probe_ccf, probe_color);
