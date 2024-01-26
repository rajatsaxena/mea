% basepath
basepath = 'X:\SWIL-Exp-Rajat\Spikesorted-SWIL';
% files
fnames = {'SWIL105PPC', 'SWIL11PPC', 'SWIL12PPC', 'SWIL13PPC', 'SWIL15PPC',  ...
    'SWIL18PPC', 'SWIL19PPC', 'SWIL20PPC', 'SWIL22PPC', 'SWIL23PPC', ...
    'SWIL24PPC', 'SWIL25PPC', 'SWIL26PPC', 'SWIL105VC', 'SWIL11VC', ... 
    'SWIL12VC', 'SWIL13VC', 'SWIL15VC', 'SWIL18VC', 'SWIL19VC', ... 
    'SWIL20VC', 'SWIL22VC', 'SWIL23VC', 'SWIL24VC', 'SWIL25VC', 'SWIL26VC'};
 
% iterate through files
for i=1:length(fnames)
    fname = fnames{i};
    analysisdirpath = fullfile(basepath,fname);
    cd(analysisdirpath);
    disp('Processing....');
    disp(fname);

    % create session structure and validate the session files
    session = sessionTemplate(analysisdirpath,'showGUI',false);
    % general info
    session.general.name = fname;
    session.general.basePath = analysisdirpath;
    session.general.experimenters = 'RS';
    session.general.investigator = 'BLM';
    session.general.projects = 'SWIL';
    session.general.sessionType = 'Acute';
    % animal info
    session.animal.name = fname;
    session.animal.sex = 'Male';
    session.animal.species = 'Mouse';
    session.animal.strain = 'C57B1/6';
    % extracellular info
    if exist(fullfile(analysisdirpath, 'temp_wh.dat'), 'file')
        session.extracellular.fileName = 'temp_wh.dat';
        session.extracellular.nChannels = length(session.extracellular.spikeGroups.channels{1}) + length(session.extracellular.spikeGroups.channels{2});
    else
        session.extracellular.fileName = strcat(fname,'merged.bin');
    end
    % validate session data
    validateSessionStruct(session);

    % load spikes and waveforms for good cluster
    spikes = loadSpikes('session',session);
    
    % read cluster info
    cluinfo = tdfread(fullfile(analysisdirpath,'cluster_info.tsv'));

    % calculate waveform metrics
    sr = session.extracellular.sr;
    waveform_metrics = calc_waveform_metrics(spikes,sr);

    % ACG and CCG metrics
    sr = session.extracellular.sr;
    acg_metrics = calc_ACG_metrics(spikes,sr);

    % fit ACGs
    fit_params = fit_ACG(acg_metrics.acg_narrow);
    close all;

    % cell-type classification
    % All cells are initially assigned as Pyramidal cells
    putativeCellType = repmat({'Pyramidal Cell'},1,length(spikes.cluID));
    % Cells are reassigned as interneurons by below criteria 
    % Narrow interneuron assigned if troughToPeak <= 0.425 ms (preferences.putativeCellType.troughToPeak_boundary)
    putativeCellType(waveform_metrics.troughtoPeak <= 0.425) = repmat({'Narrow Interneuron'},sum(waveform_metrics.troughtoPeak <= 0.425),1);
    % acg_tau_rise > 6 ms (preferences.putativeCellType.acg_tau_rise_boundary) and troughToPeak > 0.425 ms
    putativeCellType(fit_params.acg_tau_rise > 6 & waveform_metrics.troughtoPeak > 0.425) = repmat({'Wide Interneuron'},sum(fit_params.acg_tau_rise > 6 & waveform_metrics.troughtoPeak > 0.425),1);

    % save output data in a csv file
    dat.cluID = spikes.cluID';
    dat.cellType = putativeCellType';
    dat.acgRiseTime = fit_params.acg_tau_rise';
    dat.wfampCE = spikes.peakVoltage';
    dat.wfpeakToTrough = waveform_metrics.peaktoTrough';
    dat.wftroughToPeak = waveform_metrics.troughtoPeak';
    dat.wfabRatio = waveform_metrics.ab_ratio';
    dat.wfpolarity = waveform_metrics.polarity';
    dat.wfderivative = waveform_metrics.derivative_TroughtoPeak';
    dat.wfpeakA = waveform_metrics.peakA';
    dat.wfpeakB = waveform_metrics.peakB';
    dat.wftrough = waveform_metrics.trough';
    dat.tmi = acg_metrics.thetaModulationIndex';
    dat.burstRoyer = acg_metrics.burstIndex_Royer2012';
    dat.burstDoublets = acg_metrics.burstIndex_Doublets';
    dat = struct2table(dat);
    opfname = fullfile(basepath, 'analyzedMetrics', strcat(fname,'-CellExplorerUnitMetrics.csv'));
    writetable(dat,opfname);
    clear dat;
end