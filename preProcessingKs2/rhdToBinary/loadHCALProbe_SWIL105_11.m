dirnames = {'./Data-SWIL/SWILRound3/SWIL11/RawData';
                    './Data-SWIL/SWILRound3/SWIL11/SWIL11_220413_175914';
                    };
                
for j=1:length(dirnames)
    dirname = dirnames{j};
    x = dir(fullfile(dirname,'*rhd'));
    
    for i=1:length(x)
        filename = x(i).name;
        [arr1, arr2] = read_Intan_RHD2000_file_int_dualProbe_512Ch(fullfile(dirname,filename));
        
        vis = horzcat([arr2(129:end,:);arr1(129:end,:)]);
        hpc = horzcat([arr2(1:128,:);arr1(1:128,:)]);
        clear arr1 arr2
        
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