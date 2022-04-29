clc; close all; clear all;

base_filename = 'sphere.txt';
%base_filename = 'cube.txt';
%base_filename = 'cylinder.txt';

% parameter values for new test file

pml_type = 2;
RO = 1e-7;
number_of_pml_layers = 8;
tolerance = 0.001;
eigenvalue = 1;

number_of_tests = 20;

% Commands that goes into the batch file
batch_file_list = [];

% Here we are sweeping the test over tolerance value
% tolerance_list = logspace(-1, -5, number_of_tests);
eigenvalue_list = linspace(0.9,0.92,number_of_tests);

for ind = 1:number_of_tests
    
    %tolerance = tolerance_list(ind);
    eigenvalue = eigenvalue_list(ind);
    
    [filepath,filename,ext] = fileparts(base_filename);
    test_filename = [filename '_' num2str(ind) '.txt'];

    fid = create_test_file(base_filename, test_filename, pml_type, RO, ...
        number_of_pml_layers, tolerance, eigenvalue);
    
    if (fid == -1)
        disp('Base file cannot be opened!');
        return;
    end
    
    batch_file_list = strvcat(batch_file_list, ['fdfd.exe ' test_filename]);
     com = ['fdfd.exe ' test_filename]; 
     [status, result] = system(com);

end

fod = fopen('run.bat', 'w');
for ind = 1:number_of_tests
    fprintf(fod, '%s\n', batch_file_list(ind,:));
end
fclose(fod);

    

