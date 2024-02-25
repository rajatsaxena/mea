%% Used for LFP alignment if you have some channels out of the brain in cortex
% requires snippet of all LFP channels loaded in order of most superficial
% to most deep (10s of data)

%What this does: channels outside the brain are usually super correlated
%with a steep drop-off once a channel enters the brain, so this gets the
%correlation across channels, finds where the correlation with superficial
%channels drops off, and sets that as the 'start of the brain'. That gives
%an anchor point where you can match the probe to the CCF, it then finds
%where the probe is along the probe_ccf trajectory based on the length of
%the probe, and gets the CCF areas and borders across the probe.

% code taken from GitHub Andy Peters
lfp_channel_positions=sh1(:,3);

% Cortex alignment: assume top of probe out of brain and very
% correlated, look for big drop in correlation from top-down
lfp_corr = corrcoef(double(transpose(sh1_lfp-nanmedian(sh1_lfp,1)))); % compute LFP median-subtracted correlation
lfp_corr_diag = lfp_corr;
lfp_corr_diag(triu(true(size(lfp_corr)),0)) = NaN;
lfp_corr_from_top = nanmean(lfp_corr_diag,2)';

n_lfp_medfilt = 10;
lfp_corr_from_top_medfilt = medfilt1(lfp_corr_from_top,n_lfp_medfilt);
lfp_corr_from_top_medfilt_thresh = ...
    (max(lfp_corr_from_top_medfilt) - ...
    range(lfp_corr_from_top_medfilt)*0.2);
ctx_start = lfp_channel_positions( ...
    find(lfp_corr_from_top_medfilt > lfp_corr_from_top_medfilt_thresh,...
    1,'last'));

%if verbose
figure;
imagesc(lfp_channel_positions,lfp_channel_positions,lfp_corr_diag);
axis image
colormap(brewermap([],'*RdBu'));
caxis([-1,1])
xlabel('Depth (\mum)');
ylabel('Depth (\mum)');
c = colorbar;
ylabel(c,'Med. sub. correlation');
line(xlim,[ctx_start,ctx_start],'color','k','linewidth',2);
line([ctx_start,ctx_start],ylim,'color','k','linewidth',2);
title('LFP correlation and cortex start');
%end

% (give leeway for cortex to start relative to first unit - if within
% this leeway then back up cortex start, if outside then keep but throw
% a warning)
ctx_lfp_spike_leeway = 100; % um leeway for lfp/unit match
ctx_lfp_spike_diff = ctx_start-min(template_depths);
if ctx_lfp_spike_diff > ctx_lfp_spike_leeway
    warning('LFP-estimated cortex start is after first unit %.0f um', ...
        ctx_lfp_spike_diff);
elseif ctx_lfp_spike_diff > 0 && ctx_lfp_spike_diff < ctx_lfp_spike_leeway
    ctx_start = min(template_depths)-1; % (back up to 1um before first unit)
end

%%% If histology is aligned (from AP histology), get areas by depth

% (load the CCF structure tree)
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Load probe CCF alignment
load('probe_ccf.mat');
probe_ccf=probe_ccf(1); % just get one shank first

% Get area names and borders across probe
[~,dv_sort_idx] = sort(probe_ccf.trajectory_coords(:,2)); % probe trajectory Ycoords 
dv_voxel2um = 10*0.945; % CCF DV estimated scaling
probe_trajectory_depths = ...
    pdist2(probe_ccf.trajectory_coords, ...
    probe_ccf.trajectory_coords((dv_sort_idx == 1),:))*dv_voxel2um;

probe_depths = probe_trajectory_depths + ctx_start;

% Get recorded areas and boundaries
% original code ignores layer distinctions)
probe_areas = unique(regexprep( ...
    st(probe_ccf.trajectory_areas(probe_depths > 0 & ...
    probe_depths < max(channel_positions(:,2))),:).safe_name, ...
    ' layer .*',''));

% Want to get region and layer boundaries
% though this should skip layer 1
probe_areas = unique(regexprep( ...
    st(probe_ccf.trajectory_areas(probe_depths > 0 & ...
    probe_depths < max(channel_positions(:,2))),:).safe_name, ...
    ' layer ',''));

%%%%% FIX HERE
probe_area_boundaries = cellfun(@(area) prctile(probe_depths(contains(st(probe_ccf.trajectory_areas,:).safe_name,area)),[0,100]), ...
    probe_areas,'UniformOutput',false);

% orig
probe_area_boundaries = cellfun(@(area) ...
    prctile(probe_depths(contains( ...
    st(probe_ccf.trajectory_areas,:).safe_name,area)),[0,100]), ...
    probe_areas,'uni',false);
% gives boundaries from 400 to 1852 microns is S1 upper limb
