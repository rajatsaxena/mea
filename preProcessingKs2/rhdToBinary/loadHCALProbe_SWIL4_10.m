dirnames = {'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound2\SWIL4-TD1\RawData';
            'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound2\SWIL5-TD2\RawData';
            'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound2\SWIL6-TD3\RawData';
            'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound2\SWIL7-TD4\RawData';
            };
for j=1:length(dirnames)
    dirname = dirnames{j};
    x = dir(fullfile(dirname,'*rhd'));
    
    for i=1:length(x)
        filename = x(i).name;
        [vis, hpc] = read_Intan_RHD2000_file_int_dualProbe(fullfile(dirname,filename));
        
        filename(end-3:end) = '.bin';
        fname = strcat(filename(1:end-4), '_VC', filename(end-3:end));
        fid_al = fopen(fullfile(dirname,fname),'w');
        fwrite(fid_al,vis,'int16');
        fclose(fid_al);
        
        fname = strcat(filename(1:end-4), '_PPC', filename(end-3:end));
        fid_hpc = fopen(fullfile(dirname,fname),'w');
        fwrite(fid_hpc,hpc,'int16');
        fclose(fid_hpc);
    end
end