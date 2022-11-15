clc; close all; clear all;

base_filename = 'sphere.txt';
analytical_solution = 'sphere_analytical';

% parameter values for new test file
pml_type = 2;
RO = 1e-7;
number_of_pml_layers = 8;

number_of_tests = 20;

% Commands that goes into the batch file
batch_file_list = [];

% Here we are sweeping the test over tolerance value
tolerance_list = logspace(-1, -5, number_of_tests);
error_list = tolerance_list*0;
simulation_time_list = tolerance_list*0;
is_sweeping_tolerance = true;

for ind = 1:number_of_tests
    
    tolerance = tolerance_list(ind);
    
    [filepath,filename,ext] = fileparts(base_filename);
    test_filename = [filename '_' num2str(ind) '.txt'];
    rcs_theta_filename = ['farfield_theta_' filename '_' num2str(ind)];
    simulation_time_filename = ['sim_time_' filename '_' num2str(ind)];

    fid = create_test_file(base_filename, test_filename, pml_type, RO, ...
        number_of_pml_layers, tolerance);
    
    if (fid == -1)
        disp('Base file cannot be opened!');
        return;
    end
    
%     batch_file_list = strvcat(batch_file_list, ['fdfd.exe ' test_filename]);
     [status, result] = system(['fdfd.exe ' test_filename],'-echo');

     eval(analytical_solution);
     hold on;
     eval(rcs_theta_filename);
     hold off;
     
     error_list(ind) = get_rms_error(ydata, stt);
     eval(simulation_time_filename);
     simulation_time_list(ind) = simulation_time;
     
end

if (is_sweeping_tolerance)
    figure(1);
    semilogx(tolerance_list, error_list,'b-', 'linewidth',2);
    xlabel('tolerance');
    ylabel('error');
    figure(2);
    semilogx(tolerance_list, simulation_time_list,'b-', 'linewidth',2);
    xlabel('tolerance');
    ylabel('simulation time (s)');
end
