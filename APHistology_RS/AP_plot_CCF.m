function [ccf_3d_axes, ccf_axes] = AP_plot_CCF(tv,av,st,structures)
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
        'FaceColor',[0.83,0.5,0.5],'EdgeColor','none','FaceAlpha',0.1);
    % [0.8,0.5,0.5] original
    
    for s=1:length(structures)
        plot_structure = structures(s);

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

    end

end