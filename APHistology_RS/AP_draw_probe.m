function AP_draw_probe(ccf_3d_axes, ccf_axes, probe_ccf, probe_color)
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
end