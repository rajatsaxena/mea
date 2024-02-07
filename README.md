# mea
**McNaughton lab ephys data analysis.** 

**ORDER** \
extractWfAcgMetrics.m -> postProcessKs2.py -> plotPooledMetrics.py -> analyzeBehavFinal.py -> plotBehaviorAll.py -> mainRatemaps.py

Note that many of these scripts are modified from different sources such as Allen Institute, IBL, cortex-lab, Intan, NWB, etc. I've tried to give credit whenever possible.


* **preprocessingKs2**:
  * IntanRhdToNWB.py: use rhdtonwb-SWIL.py script as a template for converting Intan RHD files to NWB data format
  * ProcessRhd: used to convert Intan RHD files to merged binary files (for Kilosort) after shift alignment, and extract merged 30x sub-sampled LFP, digital, and analog signals
  * rhd_to_bin: optimized version of ProcessRhd by Rob bain to create merged binary files for raw data, sub-sampled lfp file, analog, and digital signals
    
* **postprocessingKs2**: scripts for post-processing after spike-sorting and manual curation using Kilsort and Phy
  * utilsClusterQual.py: functions for cluster quality metrics (ISI violations, presence ratio, amplitude cutoff, firing rate threshold)
  * utilsWaveformMetrics.py: functions for calculating waveform metrics (amplitude, peak-to-trough, repolarization slope, spread)
  * postProcessKs2.py: functions for running postprocess kilosort quality checks and saving final processed metrics
  * extractWfAcgMetrics.m: matlab script based on CellExplorer (Petersen et al., 2021) to extract waveform and autocorrelogram metrics
  * plotPooledMetrics.py: script to save pooled unit metrics for all animals and save waveform metrics figure for putative excitatory & inhibitory neuron
    
* **Behavior**: scripts for SWIL behavior data analysis
  * behavutils.py: supplement utils for processing behavior .mat file output from smoothwalk
  * analyzeBehavFinal.py: file to process SWIL behavior data and save npy with behavior for each environments and aligned intan timestamps
  * plotBehaviorAll.py: plot speed across all VR environments pooled across all animals

* **LFP**: scripts for SWIL lfp data analysis
  * decimateLFP_fs24.py.py: script to decimate LFPs by 24x for all recorded channels
  * mea.py: script with utils functions
  * detect_peaks.npy: script to detect peaks (similar to matlab)

* **Ratemaps** : scripts for generating rate for SWIL project
  * rmaputils.py: utils for generating ratemaps
  * mainRatemaps.py: file to generate and store ratemaps and save processed data
