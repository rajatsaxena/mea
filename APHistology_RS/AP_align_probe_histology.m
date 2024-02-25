function AP_align_probe_histology(st,slice_path, ...
    spike_times,spike_templates,template_depths, ...
    lfp,lfp_channel_positions,use_probe, figname)
% AP_align_probe_histology(st,slice_path,spike_times,spike_templates,template_depths,lfp,lfp_channel_positions,use_probe)

% If no probe specified, use probe 1
if ~exist('use_probe','var') || isempty(use_probe)
    use_probe = 1;
end

% Load probe CCF
probe_ccf_fn = [slice_path filesep 'probe_ccf.mat']; % from previous step (drawn probe)
load(probe_ccf_fn);

% Get normalized log spike n
[~,~,spike_templates_reidx] = unique(spike_templates); % num of unique templates
norm_template_spike_n = mat2gray(log10(accumarray(spike_templates_reidx,1)+1)); % orig

% Get multiunit correlation
n_corr_groups = 35; % 35 works well; % adapted for 256F probe. Orig was 40 for neuropixels, how many bins for depth
max_depths = 2125; % adapted for 256F probe. Orig was 3840 for neuropixels (hardcode, sometimes kilosort2 drops channels)
depth_group_edges = linspace(0,max_depths,n_corr_groups+1);
depth_group = discretize(template_depths,depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
unique_depths = 1:length(depth_group_edges)-1; % 35 unique depths

spike_binning = 0.01; % in seconds, binning 10 milliseconds increments from first to last spike
corr_edges = nanmin(spike_times):spike_binning:nanmax(spike_times); 
corr_centers = corr_edges(1:end-1) + diff(corr_edges);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths)
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times( ...
        ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
        corr_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

gui_fig = figure('color','w','KeyPressFcn',@keypress);

% Plot spike depth vs rate
unit_ax = subplot('Position',[0.1,0.1,0.1,0.8]);
% scatter(norm_template_spike_n,template_depths,15,'k','filled'); % dimensions not matching PROBLEM
%%% THIS IS A MAJOR PROBLEM ** especially if you subset only the good
%%% clusters

% really shitty fix
scatter(norm_template_spike_n,template_depths(1:length(norm_template_spike_n)),15,'k','filled'); % dimensions not matching, 8 values off
set(unit_ax,'YDir','reverse'); % original
ylim([0,max_depths]);
xlabel('N spikes')
title('Template depth & rate')
set(unit_ax,'FontSize',12)
ylabel('Depth (\mum)');

% Plot multiunit correlation
multiunit_ax = subplot('Position',[0.2,0.1,0.3,0.8],'YDir','reverse');
imagesc(depth_group_centers,depth_group_centers,mua_corr); % 35x35
caxis([0,max(mua_corr(mua_corr ~= 1))]); colormap(hot);
ylim([0,max_depths]);
set(multiunit_ax,'YTick',[]);
title('MUA correlation');
set(multiunit_ax,'FontSize',12)
xlabel(multiunit_ax,'Multiunit depth'); % xaxis 0 means (superficial layers), shank tip is 2125 (deep)

% Plot LFP median-subtracted correlation
lfp_moving_median = 10; % original was 10, channels to take sliding median
lfp_ax = subplot('Position',[0.5,0.1,0.3,0.8],'YDir','reverse');
imagesc(lfp_channel_positions,lfp_channel_positions, ...
    corrcoef((movmedian(zscore(double(lfp),[],2),lfp_moving_median,1) - ...
    nanmedian(zscore(double(lfp),[],2),1))'));
% values=corrcoef((movmedian(zscore(double(lfp),[],2),lfp_moving_median,1) - nanmedian(zscore(double(lfp),[],2),1))');
% the lfp corr matrix will definitely be diff sizes 101*101 dimensions from channelpos

xlim([0,max_depths]);
ylim([0,max_depths]);
set(lfp_ax,'YTick',[]);

title('LFP correlation'); % renamed from lfp power to correlation
set(lfp_ax,'FontSize',12);
caxis([-1,1]);
xlabel(lfp_ax,'Depth (\mum)'); 
colormap(lfp_ax,brewermap([],'*RdBu'));

% Link all y-axes
linkaxes([unit_ax,multiunit_ax,lfp_ax],'y');

% Plot probe areas (interactive)
% (load the colormap - located in the repository, find by associated fcn)
allenCCF_path = fileparts(which('allenCCFbregma'));
cmap_filename = [allenCCF_path filesep 'allen_ccf_colormap_2017.mat'];
load(cmap_filename);

probe_areas_ax = subplot('Position',[0.8,0.1,0.05,0.8]);

% Convert probe CCF coordinates to linear depth (*10 to convert to um)
% (use the dorsal-most coordinate as the reference point)
[~,dv_sort_idx] = sort(probe_ccf(use_probe).trajectory_coords(:,2));

probe_trajectory_depths = ...
    pdist2(probe_ccf(use_probe).trajectory_coords, ...
    probe_ccf(use_probe).trajectory_coords((dv_sort_idx == 1),:))*10;

trajectory_area_boundary_idx = ...
    [1;find(diff(double(probe_ccf(use_probe).trajectory_areas)) ~= 0)+1];
trajectory_area_boundaries = probe_trajectory_depths(trajectory_area_boundary_idx);
trajectory_area_centers = (trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries)/2);
trajectory_area_labels = st.safe_name(probe_ccf(use_probe).trajectory_areas(trajectory_area_boundary_idx));

[~,area_dv_sort_idx] = sort(trajectory_area_centers);

image([],probe_trajectory_depths,probe_ccf(use_probe).trajectory_areas);
colormap(probe_areas_ax,cmap);
caxis([1,size(cmap,1)]);
set(probe_areas_ax,'YTick',trajectory_area_centers(area_dv_sort_idx), ...
    'YTickLabels',trajectory_area_labels(area_dv_sort_idx));
set(probe_areas_ax,'XTick',[]);
set(probe_areas_ax,'YAxisLocation','right')

ylim([0,max_depths]);
ylabel({'Probe areas','(Arrow/shift keys to move)','(Escape: save & quit)'});
set(probe_areas_ax,'FontSize',10)

% Draw boundary lines at borders (and undo clipping to extend across all)
boundary_lines = gobjects;
for curr_boundary = 1:length(trajectory_area_boundaries)
    boundary_lines(curr_boundary,1) = line(probe_areas_ax,[-13.5,1.5], ...
        repmat(trajectory_area_boundaries(curr_boundary),1,2),'color','b','linewidth',2);
end
set(probe_areas_ax,'Clipping','off');

% Package into gui
gui_data = struct;
gui_data.probe_ccf_fn = probe_ccf_fn; 

gui_data.probe_ccf = probe_ccf;
gui_data.use_probe = use_probe;

gui_data.unit_ax = unit_ax;
gui_data.multiunit_ax = multiunit_ax;

gui_data.probe_areas_ax = probe_areas_ax;
gui_data.probe_areas_ax_ylim = ylim(probe_areas_ax);
gui_data.probe_trajectory_depths = probe_trajectory_depths;

gui_data.trajectory_area_boundaries = trajectory_area_boundaries; %% ADDED THIS LINE to save out blue lines

gui_data.figname = figname;

% Upload gui data
guidata(gui_fig,gui_data);

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

% Set amounts to move by with/without shift
if any(strcmp(eventdata.Modifier,'shift'))
    y_change = 100;
else
    y_change = 1;
end

switch eventdata.Key
    
    % up/down: move probe areas
    case 'uparrow'
        new_ylim = gui_data.probe_areas_ax_ylim - y_change;
        ylim(gui_data.probe_areas_ax,new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig,gui_data);
    case 'downarrow'
        new_ylim = gui_data.probe_areas_ax_ylim + y_change;
        ylim(gui_data.probe_areas_ax,new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig,gui_data);
        
    % escape: save and quit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')
            
            probe_ccf = gui_data.probe_ccf;
            
            % Get the probe depths corresponding to the trajectory areas
            probe_depths = gui_data.probe_trajectory_depths - ...
                gui_data.probe_areas_ax_ylim(1); % how much it shifted in Y         
            
            probe_ccf(gui_data.use_probe).probe_depths = probe_depths;

            boundary_depths=gui_data.trajectory_area_boundaries - gui_data.probe_areas_ax_ylim(1); %% ADDED THIS LINE to save out blue lines

            probe_ccf(gui_data.use_probe).trajectory_area_boundaries_DV=boundary_depths; %% ADDED THIS LINE to save out blue lines

            % Save the appended probe_ccf structure
            save(gui_data.probe_ccf_fn,'probe_ccf');

            % Close the figure
            savefig(gui_data.figname);
            close(gui_fig);
            
        end
end

end


