function [fid] = ...
    create_test_file(base_filename, test_filename, pml_type, RO, ...
    number_of_pml_layers, tolerance,eigenvalue)
    fid = fopen(base_filename, 'r');
    if (fid ==-1) 
        return;
    end
    fod = fopen(test_filename, 'w');
    
    tline = fgetl(fid);
    while ischar(tline)
        fprintf(fod, [tline '\n']);
        if strfind(tline, 'pml_type')
            fprintf(fod, [num2str(pml_type) '\n']);
            tline = fgetl(fid);
        end
        if strfind(tline, 'RO')
            fprintf(fod, [num2str(RO) '\n']);
            tline = fgetl(fid);
        end
        if strfind(tline, 'number_of_pml_layers')
            fprintf(fod, [num2str(number_of_pml_layers) '\n']);
            tline = fgetl(fid);
        end
        if strfind(tline, 'tolerance')
            fprintf(fod, [num2str(tolerance) '\n']);
            tline = fgetl(fid);
        end
          if strfind(tline, 'eigenvalue')
            fprintf(fod, [num2str(eigenvalue) '\n']);
            tline = fgetl(fid);
        end
       
       tline = fgetl(fid);
    end
    fclose(fid);
    fclose(fod);
end