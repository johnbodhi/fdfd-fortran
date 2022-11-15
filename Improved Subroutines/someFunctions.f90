module someFunctions
	contains
		integer function assignParameters (string1, string2)
		use global
		implicit none
		character(80) string1, string2;
        integer tmp;

		if ( string1 .eq. 'frequency' ) then
			read(string2,*) frequency; 
		end if
		if ( string1 .eq. 'computation_domain_size_in_x_direction' ) then
			read(string2,*) Sx; 
		end if
		if ( string1 .eq. 'computation_domain_size_in_y_direction' ) then
			read(string2,*) Sy; 
		end if
		if ( string1 .eq. 'computation_domain_size_in_z_direction' ) then
			read(string2,*) Sz; 
		end if
		if ( string1 .eq. 'cell_size_in_x_direction_dx' ) then
			read(string2,*) dx; 
		end if
		if ( string1 .eq. 'cell_size_in_y_direction_dy' ) then
			read(string2,*) dy; 
		end if
		if ( string1 .eq. 'cell_size_in_z_direction_dz' ) then
			read(string2,*) dz; 
        end if
        
        if ( string1 .eq. 'waveform_amplitude' ) then
			read(string2,*) A; 
        end if
        
         if ( string1 .eq. 'surface_ref' ) then
			read(string2,*) ref; 
         end if
         
        if ( string1 .eq. 'surface_tra' ) then
			read(string2,*) tra; 
        end if
          
        if ( string1 .eq. 'reflection_plane' ) then
			read(string2,*) rp; 
        end if
           
        if ( string1 .eq. 'transmission_plane' ) then
			read(string2,*) tp; 
        end if
        
             
        if ( string1 .eq. 'source_plane' ) then
			read(string2,*) sp; 
        end if
        
        
		if ( string1 .eq. 'number_of_pml_layers' ) then
			read(string2,*) Npml; 
		end if

		if ( string1 .eq. 'initial_value_filename' ) then
			read(string2,*) initial_value_filename; 
		end if

		if ( string1 .eq. 'result_filename' ) then
			read(string2,*) result_filename; 
		end if

		if ( string1 .eq. 'geometry_filename' ) then
			read(string2,*) geometry_filename; 
		end if

		if ( string1 .eq. 'materials_filename' ) then
			read(string2,*) materials_filename; 
		end if
		if ( string1 .eq. 'run_on_gpu' ) then
			read(string2,*) tmp;
            if (tmp.eq.0) then
				run_on_gpu = .false.;
            else
				run_on_gpu = .true.;
            end if
		end if
		if ( string1 .eq. 'print_residual' ) then
			read(string2,*) print_residual; 
		end if

		if ( string1 .eq. 'geometry' ) then
			if ( string2 .eq. 'sphere' ) then
				numberOfSpheres = numberOfSpheres + 1;
			end if 
		end if
		if ( string1 .eq. 'radius' ) then
			read(string2,*) radius(numberOfSpheres); 
		end if
		if ( string1 .eq. 'center_x' ) then
			read(string2,*) centerX(numberOfSpheres); 
		end if
		if ( string1 .eq. 'center_y' ) then
			read(string2,*) centerY(numberOfSpheres); 
		end if
		if ( string1 .eq. 'center_z' ) then
			read(string2,*) centerZ(numberOfSpheres); 
        end if
    

        if ( string1 .eq. 'electric_conductivity' ) then
			read(string2,*) electric_conductivity(numberOfSpheres); 
		end if
		if ( string1 .eq. 'magnetic_conductivity' ) then
			read(string2,*) magnetic_conductivity(numberOfSpheres); 
		end if
		if ( string1 .eq. 'sphere_eps_r_real' ) then
			read(string2,*) sphere_eps_r_real(numberOfSpheres); 
		end if
		if ( string1 .eq. 'sphere_eps_r_imag' ) then
			read(string2,*) sphere_eps_r_imag(numberOfSpheres); 
		end if
		if ( string1 .eq. 'sphere_mu_r_real' ) then
			read(string2,*) sphere_mu_r_real(numberOfSpheres); 
		end if
		if ( string1 .eq. 'sphere_mu_r_imag' ) then
			read(string2,*) sphere_mu_r_imag(numberOfSpheres); 
		end if
		if ( string1 .eq. 'sphere_kappa_real' ) then
			read(string2,*) sphere_kappa_real(numberOfSpheres); 
		end if
		if ( string1 .eq. 'sphere_kappa_imag' ) then
			read(string2,*) sphere_kappa_imag(numberOfSpheres); 
		end if


		if ( string1 .eq. 'geometry' ) then
			if ( string2 .eq. 'cube' ) then
				numberOfCubes = numberOfCubes + 1;
			end if 
		end if

		if ( string1 .eq. 'cube_min_x' ) then
			read(string2,*) cube_min_x(numberOfCubes); 
		end if

		if ( string1 .eq. 'cube_min_y' ) then
			read(string2,*) cube_min_y(numberOfCubes); 
		end if

		if ( string1 .eq. 'cube_min_z' ) then
			read(string2,*) cube_min_z(numberOfCubes); 
		end if

		if ( string1 .eq. 'cube_max_x' ) then
			read(string2,*) cube_max_x(numberOfCubes); 
		end if

		if ( string1 .eq. 'cube_max_y' ) then
			read(string2,*) cube_max_y(numberOfCubes); 
		end if

		if ( string1 .eq. 'cube_max_z' ) then
			read(string2,*) cube_max_z(numberOfCubes); 
		end if

		if ( string1 .eq. 'cube_eps_r_real' ) then
			read(string2,*) cube_eps_r_real(numberOfCubes);
		end if
		if ( string1 .eq. 'cube_eps_r_imag' ) then
			read(string2,*) cube_eps_r_imag(numberOfCubes); 
		end if
		if ( string1 .eq. 'cube_mu_r_real' ) then
			read(string2,*) cube_mu_r_real(numberOfCubes); 
		end if
		if ( string1 .eq. 'cube_mu_r_imag' ) then
			read(string2,*) cube_mu_r_imag(numberOfCubes); 
		end if
		if ( string1 .eq. 'cube_kappa_real' ) then
			read(string2,*) cube_kappa_real(numberOfCubes); 
		end if
		if ( string1 .eq. 'cube_kappa_imag' ) then
			read(string2,*) cube_kappa_imag(numberOfCubes); 
		end if



		if ( string1 .eq. 'geometry' ) then
			if ( string2 .eq. 'cylinder' ) then
				numberOfcylinders = numberOfcylinders + 1;
			end if 
		end if

		if ( string1 .eq. 'cylinder_center_x' ) then
			read(string2,*) cylinder_center_x(numberOfcylinders); 
		end if

		if ( string1 .eq. 'cylinder_center_y' ) then
			read(string2,*) cylinder_center_y(numberOfcylinders); 
		end if

		if ( string1 .eq. 'cylinder_radius' ) then
			read(string2,*) cylinder_radius(numberOfcylinders); 
		end if

		if ( string1 .eq. 'cylinder_min_z' ) then
			read(string2,*) cylinder_min_z(numberOfcylinders); 
		end if

		if ( string1 .eq. 'cylinder_max_z' ) then
			read(string2,*) cylinder_max_z(numberOfcylinders); 
		end if

		if ( string1 .eq. 'cylinder_eps_r_real' ) then
			read(string2,*) cylinder_eps_r_real(numberOfcylinders);
		end if
		if ( string1 .eq. 'cylinder_eps_r_imag' ) then
			read(string2,*) cylinder_eps_r_imag(numberOfcylinders); 
		end if
		if ( string1 .eq. 'cylinder_mu_r_real' ) then
			read(string2,*) cylinder_mu_r_real(numberOfcylinders); 
		end if
		if ( string1 .eq. 'cylinder_mu_r_imag' ) then
			read(string2,*) cylinder_mu_r_imag(numberOfcylinders); 
		end if
		if ( string1 .eq. 'cylinder_kappa_real' ) then
			read(string2,*) cylinder_kappa_real(numberOfcylinders); 
		end if
		if ( string1 .eq. 'cylinder_kappa_imag' ) then
			read(string2,*) cylinder_kappa_imag(numberOfCylinders); 
		end if

		if ( string1 .eq. 'printstep' ) then
			read(string2,*) printstep; 
		end if

		if ( string1 .eq. 'number_of_iterations' ) then
			read(string2,*) number_of_iterations; 
		end if
		if ( string1 .eq. 'tolerance' ) then
			read(string2,*) tolerance; 
		end if

		assignParameters = 1;
		end function assignParameters

end module someFunctions