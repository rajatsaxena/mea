function [spike_times, spike_templates, templates, channel_positions, template_depths, template_depths_flipped] = AP_loadKs2(ksparent, aname, apfs, goodchans)
    ksdir=fullfile(ksparent,aname); % kilosort output folder
    spike_times=double(readNPY(fullfile(ksdir,'spike_times.npy'))); % spiketimes are in integers here
    spike_templates=double(readNPY(fullfile(ksdir,'spike_templates.npy'))); % vector specifying the identity of the template that was used to extract each spike
    spike_clusters=readNPY(fullfile(ksdir,'spike_clusters.npy')); % cluster ID for each spike
    % later can split single spiking data into each shank

    %%%%%%%%%%% if you want just subset the good clusters, need to also subset spike_templates
    clusterinfo = tdfread(fullfile(ksdir,'cluster_info.tsv')); % all info for each cluster 
    clustergroup = string(clusterinfo.group); 
    % Get all units except noise from Phy curation
    if goodchans
        goodclusterID=clusterinfo.cluster_id(find(clustergroup=="good "));
    else
        goodclusterID=clusterinfo.cluster_id(find(clustergroup~="noise ")); % 0 based, cluster ID
    end
    % Get spike timestamps for all good clusters
    goodclus_spikets=[]; goodclus_templates=[];
    for c=1:length(goodclusterID) 
        spikes=spike_times(spike_clusters==goodclusterID(c));
        goodclus_spikets=[goodclus_spikets; spikes./apfs]; % convert timestamps from samples to seconds

        goodtemplates=spike_templates(spike_clusters==goodclusterID(c));
        goodclus_templates=[goodclus_templates; goodtemplates];
    end
    spike_times=goodclus_spikets;
    spike_templates=goodclus_templates;
    clear goodclus_spikets goodclus_templates
    
    % load templates
    templates=double(readNPY(fullfile(ksdir,'templates.npy'))); % 3d dimension = ch across both shanks, note this is 205 ch
    
    % load channel_positions
    channel_positions = readNPY(fullfile(ksdir,'channel_positions.npy'));
    
    % TO GET TEMPLATE DEPTHS
    % template shape: [nTemplates, nTimePoints, nTempChannels]
    % Get depth of each template
    % (get min-max range for each channel)
    template_chan_amp = squeeze(range(templates,2));
    % (zero-out low amplitude channels)
    template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
    template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
    % (get center-of-mass on thresholded channel amplitudes)
    template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2); % note these channel Y pos is UNSORTED
    template_depths_flipped = 2125 - template_depths; % flip depths, this uses ch_pos

end