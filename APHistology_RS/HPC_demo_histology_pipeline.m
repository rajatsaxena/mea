% Example pipeline for processing histology
% taken from Andy Peters, modified by Winny Ning

%% 1) Load CCF and set paths for slide and slice images
cd 'W:\Shared drives\Bruce McNaughton Lab\MATLAB Code';
addpath(genpath('npy-matlab-master'));
addpath(genpath('Allen_Histology')); 
addpath(genpath('brewermap'));

% Load CCF atlas
allen_atlas_path = 'W:\Shared drives\Bruce McNaughton Lab\MATLAB Code\Allen_Histology'; % path to CCF atlas
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area
st = AP_loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % table of all labels

% Set paths for histology images and directory to save slice/alignment
dir='Y:\Winny\Histology'; 
[file,im_path]=uigetfile('*.tif','Select a histology image',dir); % im_path = path to folder with images goes here
addpath(im_path); % make image is in matlab path
slice_path = [im_path(1:end-1) filesep 'slices']; % later on will make slices folder in path

%% 2) Preprocess slide images to produce slice images

% Set white balance and resize slide images, extract slice images
% (Note: this resizes the images purely for file size reasons - the CCF can
% be aligned to histology no matter what the scaling. If pixel size is
% available in metadata then automatically scales to CCF resolution,
% otherwise user can specify the resize factor as a second argument)

% Set resize factor
% resize_factor = []; % (slides ome.tiff: auto-resize ~CCF size 10um/px)
resize_factor = 1; % (slides tiff: resize factor)

% Set slide or slice images
% slice_images = false; % (images are slides - extract individual slices)
slice_images = true; % (images are already individual slices)

% Preprocess images
AP_process_histology(im_path,resize_factor,slice_images);

% (optional) Rotate, center, pad, flip slice images 
% Draw reference line at midline to rotate/center image a bit
AP_rotate_histology(slice_path);

%% 3) Align CCF to slices

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path); % manually match AP coronal slice, can also match angles if slicing unevenly
% 5.34 S1,  7.46 HPC for TR8
% TR10 HPC 8.10, TR11 HPC 7.25

% Align CCF slices and histology slices
% (first: automatically, by outline)
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

%% 5) Align histology to electrophysiology data

% This will match trajectory of probe in CCF with ephys signatures.
% Requires kilosort output of spike_times, spike_templates and template_depths.
ksdir=uigetdir('C:\Users\Justin\Desktop'); % kilosort output folder
% ksdir='Y:\Winny\Data\TR8_2-22-2023\TR8_HPC_SpikeSorted_N=2_Batch192_RS';
addpath(genpath(ksdir));
cd(ksdir);

spike_times_samples=double(readNPY('spike_times.npy')); % spiketimes are in integers here
% NOTE: this is all the raw spike_times from Kilosort, regardless of whether it's a 'good' cluster
% This contains spiking data from BOTH shanks 
Fs=30000.0; %30 kHz sampling rate 
spike_times=spike_times_samples./Fs; % convert spike timestamps from samples to seconds
spike_templates=double(readNPY('spike_templates.npy')); % vector specifying the identity of the template that was used to extract each spike
% later can split single spiking data into each shank

%%% spike_templates are assumed that its 1-indexed!!!
% check out: https://github.com/cortex-lab/spikes/blob/master/analysis/templatePositionsAmplitudes.m


%%%%%%%%%%% if you want just subset the good clusters, need to also subset spike_templates
clusterinfo = tdfread('cluster_info.tsv'); % all info for each cluster 
clustergroup = string(clusterinfo.group); 
% Get only good units from Phy curation
goodclusteridx = find(clustergroup=="good    "); % index 
goodclusterID=clusterinfo.cluster_id(goodclusteridx); % 0 based, cluster ID

spikecluster=readNPY('spike_clusters.npy'); % cluster ID for each spike

% Get spike timestamps for all good clusters
Fs=30000.0; %30 kHz sampling rate 
goodclus_spikets=[]; goodclus_templates=[];

for c=1:length(goodclusterID) % 122 HPC units for TR8
    spikes=spike_times_samples(spikecluster==goodclusterID(c));
    goodclus_spikets=[goodclus_spikets; spikes./Fs]; % convert timestamps from samples to seconds

    goodtemplates=spike_templates(spikecluster==goodclusterID(c));
    goodclus_templates=[goodclus_templates; goodtemplates];
end

spike_times=goodclus_spikets;
spike_templates=goodclus_templates;

%%%%%%%%%%%%%%%%

load('channel_numbers.mat'); % chnum, correct channel numbers (index starts from 0), this will be less than 256 if you disable channels in Kilosort
channel_positions=readNPY('channel_positions.npy'); % col 1 is x coords, col 2 is y coords, the DV will be flipped later on!! 
% this will be flipped in template_depths later on
channelmap=[];
channelmap(:,1)=chnum; % channel number
channelmap(:,2)=channel_positions(:,1); % x coords
channelmap(:,3)=channel_positions(:,2); % y coords

% channelmap(:,2)= 2125 - channelpositions(:,2); % reverse the ycoords, so that 0 is superficial, 2125 is deepest (tip of shank)
% figure; scatter(channelmap(:,2),channelmap(:,3)); %set(gca,'Ydir','reverse'); 

lfpdir=uigetdir('Y:\Winny\Data'); % kilosort output folder
cd(lfpdir);
alllfp=double(readNPY('TR8_lfp.npy')); % load in lfp data
lfp=alllfp(257:end,:); % only take PPC/HPC which is LAST 256ch (257:512)
lfpts=double(readNPY('TR8_lfpts.npy')); % lfp timestamps (s)

load('trialts.mat'); 
timeidx=[trialts{1}(1) trialts{end}(end)]; % take a time snippet, entire behavior epoch
indices=find(lfpts>=timeidx(1) & lfpts<=timeidx(end));
lfpts_beh=lfpts(indices);
lfp_beh=lfp(:,indices); % for all 256 channels in HPC, for a time snippet


% Next, need to sort LFP channels by depth and shank
% Split by shank first
% for 256F HPC contacts are facing you in RH, so shank 1 (1 to 128) is most lateral  
% SHANK 1
s1ind=find(channel_positions(:,1)<=20); % INDICES, use xcoords to find most lateral shank, only 101 ch, rely on x position (NOT Channel map...
% if you disable ch in kilosort, chmap will be very wrong, just consecutive numbers
sh1=[];
sh1(:,1)=chnum(s1ind,:); % channel number on shank 1
sh1(:,2:3)=channel_positions(s1ind,:); % XY coords on shank 1
sh1_lfp=lfp_beh(sh1(:,1)+1,:); % lateral HPC shank, this index needs to be CHAN NUM from 0 to 127 *** but need to +1 for index, can't index from 0
% in order from 0 to 127

% SHANK 2, most medial HPC 
s2ind=find(channel_positions(:,1)>20); % indices for MEDIAL HPC shank
sh2=[];
sh2(:,1)=chnum(s2ind,:); % 1st col is channel num
sh2(:,2:3)=channel_positions(s2ind,:); % 2nd col is X and 3rd col is Y (depth)
sh2_lfp=lfp_beh(sh2(:,1)+1,:); % medial HPC shank, this index is CHAN NUM, need to +1 for index because it starts at 0


% For each shank, sort channels according to DEPTH, depth is not flipped
% yet
[~,ind1]=sort(sh1(:,3),'descend'); % col 3 (Y) is DV depth, sort from superficial (2125, 1st row) to deep (0), same as kilosort (**DV will flip later)
sh1_sorted=sh1(ind1,:); % Sorted shank 1
sh1_lfp_sorted=sh1_lfp(ind1,:); %%% need to sort this by depth index because sh1_lfp was sorted from 0 to 127

[~,ind2]=sort(sh2(:,3),'descend');
sh2_sorted=sh2(ind2,:);
sh2_lfp_sorted=sh2_lfp(ind2,:);

% To plot sorted LFP, from superficial on top to deep 
tidx=[19 20];
figure; 
for i=1:size(sh1_lfp_sorted,1)
    plot(lfpts_beh(1000*tidx(1):1000*tidx(2)),sh1_lfp_sorted(i,[1000*tidx(1):1000*tidx(2)])-(i*500)); hold on; % sup lfp on top
end
xlim([lfpts_beh(1000*tidx(1)) lfpts_beh(1000*tidx(2))]);


templates=double(readNPY('templates.npy')); % 3d dimension = ch across both shanks, note this is 205 ch
% templateidx=readNPY('templates_ind.npy'); % [nTemplates, nTempChannels], this just goes from 0 to nchannels-1 (useless)

%%% TO GET TEMPLATE DEPTHS 
% Get depth of each template
% (get min-max range for each channel)
template_chan_amp = squeeze(range(templates,2));
% (zero-out low amplitude channels)
template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
% (get center-of-mass on thresholded channel amplitudes)
template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2); % note these channel Y pos is UNSORTED
template_depths_flipped = 2125 - template_depths; % flip depths, this uses ch_pos

use_probe = 1; % choose shank 1 (lateral), need to run separately for each shank
lfp_sh = sh1_lfp_sorted; % choose LFP from one shank, LATERAL SHANK
Ychannel_positions_sh = 2125-(sh1_sorted(:,3)); % flipped ycoords from either shank 1 or 2, this is already sorted Y position from superficial (0) to deep (2125)
% might start around 400 microns if you have disabled channels

% use_probe = 2; % choose shank 2 (medial HPC)
% lfp_sh = sh2_lfp; % choose LFP from one shank
% Ychannel_positions_sh = 2125-(sh2(:,3)); % flipped ycoords from either shank 1 or 2


AP_align_probe_histology(st,slice_path, ...
    spike_times,spike_templates,template_depths_flipped, ...
    lfp_sh,Ychannel_positions_sh, ...
    use_probe);

% can move up and down to adjust position of histology with ephys data and
% then hit esc to save and quit. this will only update the probe_ccf.mat
% output file, which is the only one we care about

%% NOT NEEDED, DON'T RUN THIS: Extract slices from full-resolution images
% (not worth it at the moment, each slice is 200 MB)
%AP_grab_fullsize_histology_slices(im_path)

% this is just a general function to convert histology coords to CCF (ex.
% mark injection site on histology and want to translate into CCF coords)
histology_points  = {[probe_ccf(1).trajectory_coords(:,1),probe_ccf(1).trajectory_coords(:,2)]} ; % dunno what the inputs are????
% cell contains n points x 2 ([x,y]). E.g., for two slices: 
% {[100,200],[]} will convert the point x = 100, y = 200 on the first
% slide into CCF coordinates

% this transforms image coordinates to CCF, takes a cell array of points
% corresponding to each slice and transforms it to CCF space. Affine
% transform, same linear transform across the whole slice

% Convert points in histology images to CCF coordinates
ccf_points = AP_histology2ccf(histology_points,slice_path);
% Concatenate points and round to nearest integer coordinate
ccf_points_cat = round(cell2mat(ccf_points));
% Get indicies from subscripts
ccf_points_idx = sub2ind(size(av),ccf_points_cat(:,1),ccf_points_cat(:,2),ccf_points_cat(:,3));
% Find annotated volume (AV) values at points
ccf_points_av = av(ccf_points_idx);
% Get areas from the structure tree (ST) at given AV values
ccf_points_areas = st(ccf_points_av,:).safe_name;


%%%% To convert probe_ccf.mat into .npy for IBL ephys alignment
% Code from Andy Peters
% Pick filename
probe_ccf = load('D:\repos\AP_histology\sample slides\CA3-MEC-RE-Fibre-1122\Save Slices\probe_ccf.mat')
%filename = uiputfile('probe_ccf.npy');
filename = 'D:\repos\AP_histology\sample slides\CA3-MEC-RE-Fibre-1122\Save Slices\probe_ccf.npy';
% Write probe_ccf coordinates as NPY file
writeNPY(probe_ccf.trajectory_coords,filename)

% Convert probe_ccf.mat into .npy for IBL ephys alignment for each probe
% Pick filename
[filename, path] = uiputfile;
% Write probe_ccf coordinates as NPY file (separately for each probe)
for curr_probe = 1:length(probe_ccf)
    curr_filename = sprintf('%s%s_probe%d.npy',path,filename,curr_probe);
    writeNPY(probe_ccf(curr_probe).trajectory_coords,curr_filename);
end

%% 6) Analyze the output mat files and get the laminar boundaries
% returns 3 mat files in slices folder: probe_ccf, histology_ccf,
% atlas2histology_tform

addpath(genpath(slice_path));
load('probe_ccf.mat');

% probe_ccf.trajectory_coords = 3D CCF coords of the probe trajectory
% probe_ccf.trajectory_areas = annotated CCF areas for each point
% probe_ccf.probe_depths = relative depth of the probe to that point (0 is
% superficial (headstage) and 2125 is tip of probe (deepest layers). Note
% this is reversed from Kilosort output

CCFregions={probe_ccf.trajectory_areas}; % each shank is in a cell, need to decipher the brain region with the tree
probedepth={probe_ccf.probe_depths};

regions=st(probe_ccf(1).trajectory_areas(:),:).safe_name; % this gets the regions for all 625 values 

%%%% run this part only 
boundaries_depth={probe_ccf.trajectory_area_boundaries_DV}; % in microns this is across the whole column chunk, need to subset
figure;
shanknum=1;
for i=1:length(boundaries_depth{shanknum})
    yline(boundaries_depth{shanknum}(i),'r','LineWidth',2); hold on;
end
set(gca,'Ydir','reverse');
ylim([0  2125]);

% Only pick out regional boundaries that are from 0 to 2125 microns (span of probe)
endbound=boundaries_depth{shanknum}(find(boundaries_depth{shanknum}>2125,1)); % go 1 more than 2125 to span the range of electrodes
HPC_boundaries=cellfun(@(x) x(x<=endbound),boundaries_depth,'UniformOutput',false);  

%% 7) Load unitdata and loop through cells to categorize units that are superfical vs deep
% load('TR8_HPCunitdata');
counter1=0; counter2=0; counter3=0; counter4=0; counter5=0; counter6=0; counter7=0; counter8=0;

for c=1:length(unitdata)
    unitdepth=2125-unitdata(c).depth; % flip the depths so that 0 is superficial, 2125 is deep

    if unitdepth<=HPC_boundaries{1}(2)
        unitdata(c).region='L1'; % PPC layer 1
        counter1=counter1+1;

    elseif unitdepth<=HPC_boundaries{1}(3)
        unitdata(c).region='L2/3'; % PPC layer 2/3
        counter2=counter2+1;

    elseif unitdepth<=HPC_boundaries{1}(4)
        unitdata(c).region='L4'; % PPC layer 4
        counter3=counter3+1;

    elseif unitdepth<=HPC_boundaries{1}(5)
        unitdata(c).region='L5'; % PPC limb layer 5
        counter4=counter4+1;

    elseif unitdepth<=HPC_boundaries{1}(6) || unitdepth<=HPC_boundaries{1}(7)
        unitdata(c).region='L6'; % PPC layer 6a and 6b
        counter5=counter5+1;

    elseif unitdepth<=HPC_boundaries{1}(8) ||  unitdepth<=HPC_boundaries{1}(9) || unitdepth<=HPC_boundaries{1}(10)
        unitdata(c).region='CC'; % Corpus callosum
        counter6=counter6+1;

    elseif unitdepth<=HPC_boundaries{1}(11)
        unitdata(c).region='CA1'; % HPC CA1
        counter7=counter7+1;

    elseif unitdepth<=HPC_boundaries{1}(12)
        unitdata(c).region='DG'; % DG
        counter8=counter8+1;
    end
end

sum([counter1,counter2,counter3,counter4,counter5,counter6,counter7,counter8])

% 0 units in layer 1
% 10 units in layer 2/3
% 6 units in layer 4
% 64 units in layer 5
% 33 units in layer 6
% 0 units in corpus callosum
% 8 units in CA1 HPC
% 1 unit in DG

%% Another (better) method of assigning region and layer for each unit
for c=1:length(unitdata)
    unitdepth=2125-unitdata(c).depth;
    idx=interp1(probe_ccf(1).probe_depths',1:length(probe_ccf(1).probe_depths),unitdepth,'nearest'); % find closest probe depth
    areaID=probe_ccf(1).trajectory_areas(idx); % num for the region
    unitdata(c).regname=cell2mat(st.acronym(areaID)); % get region name from st
end

%% Plot pie chart of how many cells in each layer
cell_counts=[counter1; counter2; counter3; counter4; counter5; counter6; counter7; counter8]; 
sum(cell_counts); % 122 HPC cells
idx=find(cell_counts~=0);

pie_labels = {};
figure; p=pie(cell_counts(idx),'%.1f%%');  % p=pie(cell_counts(idx),'%.1f%%'); 
title('PPC and HPC cells'); 
pText = findobj(p,'Type','text'); 
percentValues = get(pText,'String'); 
pie_labels = {'L1: '; 'L2/3: '; 'L4: '; 'L5: '; 'L6: '; 'CC: '; 'CA1: '; 'DG: '};
    
combinedtxt = strcat(pie_labels(idx),percentValues); 
for x=1:size(combinedtxt,1)
    pText(x).String = combinedtxt(x);
    pText(x).FontSize=13;
end

set(gca,'FontSize',14);


X = categorical({'L1'; 'L2/3'; 'L4'; 'L5'; 'L6'; 'CC'; 'CA1'; 'DG'});
X = reordercats(X,{'L1'; 'L2/3'; 'L4'; 'L5'; 'L6'; 'CC'; 'CA1'; 'DG'});
figure; bar(X,cell_counts); box off; xlabel('Brain Region');
ylabel('Number of cells'); title('TR8'); set(gca,'FontSize',14);

%% Save out unitdata struct with region labels
dir='Y:\Winny\Data\TR8_2-22-2023';
myname='TR8_HPCunitdata_8-13-2023';
save(fullfile(dir,myname),'unitdata');

%% 8) Plot figure of CCF brain structures in 2D and 3D with probe tracts 
% (code from Muhang Li, Riken; modified by Winny Ning)

% Load CCF atlas
allen_atlas_path = 'W:\Shared drives\Bruce McNaughton Lab\MATLAB Code\Allen_Histology'; % path to CCF atlas
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area
st = AP_loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % table of all labels

% Draw CCF structures, 3-view outline of the brain and selected structures

% Set up figure axes, plot brain outline
figure('Color','w');
% Set up 2D axes
ccf_axes = gobjects(3,1);
ccf_axes(1) = subplot(1,4,1,'YDir','reverse'); hold on; axis image off;
ccf_axes(2) = subplot(1,4,2,'YDir','reverse'); hold on; axis image off;
ccf_axes(3) = subplot(1,4,3,'YDir','reverse'); hold on; axis image off;

% Set up 3D axes
ccf_3d_axes = subplot(1,4,4);
set(ccf_3d_axes,'ZDir','reverse'); hold(ccf_3d_axes,'on');
axis vis3d equal off manual
view([-30,25]); axis tight;
h = rotate3d(ccf_3d_axes); h.Enable = 'on';

% Plot 3D/2D brain outlines
% (2D Plot)
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(av,[],curr_view)) > 1));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2),x(:,1),'k','linewidth',2),curr_outline)
    % (draw 1mm scalebar)
    %     line(ccf_axes(curr_view),[0,0],[0,100],'color','k','linewidth',2);
end
linkaxes(ccf_axes);

% (3D Plot)
% wireframe
% [~, brain_outline] = plotBrainGrid([],ccf_3d_axes);

% mesh
slice_spacing = 5;
brain_volume = ...
    bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
    1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
brain_outline_patchdata = isosurface(permute(brain_volume,[3,1,2]),0.5);
brain_outline = patch( ...
    ccf_3d_axes, ...
    'Vertices',brain_outline_patchdata.vertices*slice_spacing, ...
    'Faces',brain_outline_patchdata.faces, ...
    'FaceColor',[0.8,0.5,0.5],'EdgeColor','none','FaceAlpha',0.1);


%%% Get brain areas in search, then draw
% Primary somatosensory area upper limb
% Posterior parietal association
% Field CA1

% Prompt for which structures to show (only structures which are
% labelled in the slice-spacing downsampled annotated volume)
structure_search = lower(inputdlg('Search structures'));
structure_match = find(contains(lower(st.safe_name),structure_search)); % find structure name in st.safe_name

selected_structure = listdlg('PromptString','Select a structure to plot:', ...
    'ListString',st.safe_name(structure_match),'ListSize',[520,500], ...
    'SelectionMode','single');

plot_structure = structure_match(selected_structure);

% Get all areas within and below the selected hierarchy level
plot_structure_id = st.structure_id_path{plot_structure};
plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
    st.structure_id_path));

% Get structure color and volume
slice_spacing = 5;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{plot_structure},2,[])')./255;
plot_ccf_volume = ismember(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);

% Plot 2D structure
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2)*slice_spacing, ...
        x(:,1)*slice_spacing,'color',plot_structure_color,'linewidth',2),curr_outline)
end

% Plot 3D structure
structure_3d = isosurface(permute(plot_ccf_volume,[3,1,2]),0);
structure_alpha = 0.2;
patch(ccf_3d_axes, ...
    'Vertices',structure_3d.vertices*slice_spacing, ...
    'Faces',structure_3d.faces, ...
    'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);

% sgtitle(string(structure_search{1})); % label title with region name


%%% Draw probes (output from AP_histology)

% load probe_ccf mat file
[probe_ccf_fn,probeccf_path] = uigetfile('*.mat','Pick probe histology file','Y:\Winny\Histology');
addpath(genpath(probeccf_path));
load(probe_ccf_fn);
probe_color = 'r';

% Loop through probes and draw
for curr_probe = 1:length(probe_ccf)

    % Get line of best fit through mean of marked points
    probe_coords_mean = mean(probe_ccf(curr_probe).points,1);
    xyz = bsxfun(@minus,probe_ccf(curr_probe).points,probe_coords_mean);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);

    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end

    % Evaluate line of best fit (length of probe to deepest point)
    [~,deepest_probe_idx] = max(probe_ccf(curr_probe).points(:,2));
    probe_deepest_point = probe_ccf(curr_probe).points(deepest_probe_idx,:);
    probe_deepest_point_com_dist = pdist2(probe_coords_mean,probe_deepest_point);
    probe_length_ccf = 2125/10; % mm / ccf voxel size, orig was 3840/10 %%%% CCF scaling is 10 microns per pixel
    % changed this part

    probe_line_eval = probe_deepest_point_com_dist - [probe_length_ccf,0];
    probe_line = (probe_line_eval'.*histology_probe_direction') + probe_coords_mean;

    % Draw probe in 3D view
    line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
        'linewidth',2,'color',probe_color)

    % Draw probes on coronal + saggital
    line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',probe_color);
    line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',probe_color);

    % Draw probe mean on horizontal
    plot(ccf_axes(2), probe_coords_mean(:,3),probe_coords_mean(:,1), ...
        '.','MarkerSize',20,'color',probe_color);
end
 
hold on; % repeat for the next region HPC



%% If you want to use all 128 channels INCLUDING disabled LFP ch for the figure
LFP=lfp_beh(1:128,:); % pick just shank 1
load('Y:\Winny\Data\chmapfull.mat'); % full channel map
chmapfullshank1=chmapfull(1:128,:,:); % pick shank 1 info
[~,ind1]=sort(chmapfullshank1(:,3),'descend'); % col 3 (Y) is DV depth, sort from superficial (2125) to deep (0), same as kilosort (**DV will flip later)
chmapsorted=chmapfullshank1(ind1,:); % Sorted shank 1
shank1lfp=LFP(ind1,:);
lfp_sh=shank1lfp; % this is all lfp 
Ychannel_positions_sh = 2125-(chmapsorted(:,3)); % flipped ycoords from either shank 1 or 2, this is already sorted Y position from superficial (0) to deep (2125)

%% Ignore below: Method to get table with cluster ID, depth, and corresponding area name (NOT working yet, requires lots of data formatting)
histinfo={probe_ccf.trajectory_areas}; % probe track coordinates?
AllenCCFPath=allen_atlas_path;
sp=loadKSdir(ksdir); % issue no RecSes in this
clusinfo=tdfread('cluster_info.tsv'); % all info for each cluster 
removenoise=1; % check this??
surfacefirst=1; % superficial surface = lowest index (default zero). Deeper in brain = highest index
LFPDir=[]; % can assign lfp folder, if empty don't use
treeversion=2;
% treeversion: which Allen Brain Tree structure to use? (default 2, = 2017; 1 = older)
trackcoordinates={probe_ccf.trajectory_coords};
% trackcoordinates: for as many datapoints as in histinfo, the X/Y/Z
% coordinates of the probe track in CCF

% function from Enny van Beest
Depth2AreaPerUnit = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,removenoise,surfacefirst,LFPDir,treeversion,trackcoordinates);








