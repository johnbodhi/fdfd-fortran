
!==============================================================================
subroutine set_geometry
	use global
	implicit none
	integer mi,mj,mk, tmi,tmj,tmk, gNx, gNy, gNz, ierr, ind, med_ind;
	character(80) string1, string2;
	integer offset_x, offset_y, offset_z;

	open (1,FILE = geometry_filename,IOSTAT=ierr, ACTION ='READ');
	if (ierr .gt. 0) then
			print 403, geometry_filename;
403   FORMAT ('Project ', A, ' can not be opened');
	stop;
	end if

	read (1, *, iostat = ierr ) string1;
	if ((string1.eq."Nx") .and. (ierr .eq. 0 ))then
		read (1, *, iostat = ierr ) string1;
		read(string1,*) gNx; 
	else
		print 404, initial_value_filename;
		stop;
	endif

	read (1, *, iostat = ierr ) string1;
	if ((string1.eq."Ny") .and. (ierr .eq. 0 )) then
		read (1, *, iostat = ierr ) string1;
		read(string1,*) gNy; 
	else
		print 404, initial_value_filename;
		stop;
	endif

	read (1, *, iostat = ierr ) string1;
	if ((string1.eq."Nz") .and. (ierr .eq. 0 )) then
		read (1, *, iostat = ierr ) string1;
		read(string1,*) gNz; 
	else
		print 404, initial_value_filename;
		stop;
	endif

	offset_x = 8;
	offset_y = 8;
	offset_z = 8;
	read (1, *, iostat = ierr ) string1;
	if ((string1.eq."Medium") .and. (ierr .eq. 0 )) then
		do mi = 1,gNx,1
			do mj = 1,gNy,1
				do mk = 1,gNz,1
					read(1,212) med_ind;
					Medium(mi+offset_x,mj+offset_y,mk+offset_z) = med_ind + 1;
					!Add 1 to  med_ind because in the model index starts from zero
				end do
			end do
		end do 
	else
		print 404, initial_value_filename;
		stop;
	endif

	close(1);

212   FORMAT (I5);
404   FORMAT ('Invalid file : ', A);

end subroutine set_geometry
!==============================================================================

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

		if ( string1 .eq. 'RO' ) then
			read(string2,*) RO; 
        end if

        if ( string1 .eq. 'pml_type' ) then
			read(string2,*) pml_type; 
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

!==============================================================================

!==============================================================================
subroutine setPMLboundaries
	use global
	implicit none
	integer mi,mj,mk;
	real sigma, dpml;

! To set the boundaries PML
    !pml_type=2; 
    ! pml_type = 0 constant cond., pml_type = 1 linear cond., pml_type = 2 parabolic cond.
	
    ! RO=1e-7;
    
    dpml=Npml*dx;  sigma_max=-eps_o*c*(pml_type+1)*log(RO)/(2*dpml);

! set the x=0 PML
do mi=0,Npml,1
            sigma=sigma_max * ((mi*dx)/dpml)**(pml_type+1);
            eps_yx (Npml-mi+1,:,:) = eps_o - j * sigma/ w;     
            eps_zx (Npml-mi+1,:,:) = eps_o - j * sigma/ w;     
end do
do mi=0,Npml-1,1
            sigma=sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
            mu_yx (Npml-mi,:,:) =  (eps_o - j * sigma/ w) * mu_o / eps_o;
            mu_zx (Npml-mi,:,:) =  (eps_o - j * sigma/ w) * mu_o / eps_o;
end do

!  set the x=N boundary PML
do mi=0,Npml,1
            sigma=sigma_max * ((mi*dx)/dpml)**(pml_type+1);
            eps_yx (mi+Nx-Npml,:,:) = eps_o - j * sigma / w;     
            eps_zx (mi+Nx-Npml,:,:) = eps_o - j * sigma / w;     
end do   
do mi=0,Npml-1,1
            sigma=sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
            mu_yx (mi+Nx-Npml,:,:) = (eps_o - j * sigma/ w) * mu_o / eps_o;
            mu_zx (mi+Nx-Npml,:,:) = (eps_o - j * sigma/ w) * mu_o / eps_o;
end do   
! set the y=0 PML
do mi=0,Npml,1
            sigma=sigma_max * ((mi*dx)/dpml)**(pml_type+1);
            eps_xy (:,Npml-mi+1,:) = eps_o - j * sigma/ w;     
            eps_zy (:,Npml-mi+1,:) = eps_o - j * sigma/ w;     
end do
do mi=0,Npml-1,1
            sigma=sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
            mu_xy (:,Npml-mi,:) = (eps_o - j * sigma/ w) * mu_o / eps_o;
            mu_zy (:,Npml-mi,:) = (eps_o - j * sigma/ w) * mu_o / eps_o;
end do

!  set the y=N boundary PML
do mi=0,Npml,1
            sigma=sigma_max * ((mi*dx)/dpml)**(pml_type+1);
            eps_xy (:,mi+Ny-Npml,:) = eps_o - j * sigma / w;     
            eps_zy (:,mi+Ny-Npml,:) = eps_o - j * sigma / w;     
end do   
do mi=0,Npml-1,1
            sigma=sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
            mu_xy (:,mi+Ny-Npml,:) = (eps_o - j * sigma/ w) * mu_o / eps_o;
            mu_zy (:,mi+Ny-Npml,:) = (eps_o - j * sigma/ w) * mu_o / eps_o;
end do   
! set the z=0 PML
do mi=0,Npml,1
            sigma=sigma_max * ((mi*dx)/dpml)**(pml_type+1);
            eps_xz (:,:,Npml-mi+1) = eps_o - j * sigma/ w;     
            eps_yz (:,:,Npml-mi+1) = eps_o - j * sigma/ w;     
end do
do mi=0,Npml-1,1
            sigma=sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
            mu_xz (:,:,Npml-mi) = (eps_o - j * sigma/ w) * mu_o / eps_o;
            mu_yz (:,:,Npml-mi) = (eps_o - j * sigma/ w) * mu_o / eps_o;
end do
!  set the z=N boundary PML
do mi=0,Npml,1
            sigma=sigma_max * ((mi*dx)/dpml)**(pml_type+1);
            eps_xz (:,:,mi+Nz-Npml) = eps_o - j * sigma / w;     
            eps_yz (:,:,mi+Nz-Npml) = eps_o - j * sigma / w;     
end do   
do mi=0,Npml-1,1
            sigma=sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
            mu_xz (:,:,mi+Nz-Npml) = (eps_o - j * sigma/ w) * mu_o / eps_o;
            mu_yz (:,:,mi+Nz-Npml) = (eps_o - j * sigma/ w) * mu_o / eps_o;
end do   
!  ===========================================

end subroutine setPMLboundaries

!==============================================================================
subroutine setCoefficients
	use global 
	implicit none
    integer mi, mj, mk
    
	allocate(Cexhy(Nx,Ny,Nz), Cexhz(Nx,Ny,Nz), Cexex(Nx,Ny,Nz));
	allocate(Ceyhz(Nx,Ny,Nz), Ceyhx(Nx,Ny,Nz), Ceyey(Nx,Ny,Nz));
	allocate(Cezhx(Nx,Ny,Nz), Cezhy(Nx,Ny,Nz), Cezez(Nx,Ny,Nz));

	allocate(Chxey(Nx,Ny,Nz), Chxez(Nx,Ny,Nz), Chxhx(Nx,Ny,Nz));
	allocate(Chyez(Nx,Ny,Nz), Chyex(Nx,Ny,Nz), Chyhy(Nx,Ny,Nz));
	allocate(Chzex(Nx,Ny,Nz), Chzey(Nx,Ny,Nz), Chzhz(Nx,Ny,Nz));


	allocate(tmpx(Nx,Ny,Nz), tmpy(Nx,Ny,Nz), tmpz(Nx,Ny,Nz));
	allocate(Bex(Nx,Ny,Nz), Bey(Nx,Ny,Nz), Bez(Nx,Ny,Nz));
	allocate(Bhx(Nx,Ny,Nz), Bhy(Nx,Ny,Nz), Bhz(Nx,Ny,Nz));
	
	tmpx = (0.0,0.0); tmpy = (0.0,0.0); tmpz = (0.0,0.0);
	
    Cexhy = 1/(j*w*dz*eps_xz); Cexhz = 1/(j*w*dy*eps_xy); Cexex = -(eps_xi - eps_o) / eps_xi;  
    Ceyhz = 1/(j*w*dx*eps_yx); Ceyhx = 1/(j*w*dz*eps_yz); Ceyey = -(eps_yi - eps_o) / eps_yi;  
    Cezhx = 1/(j*w*dy*eps_zy); Cezhy = 1/(j*w*dx*eps_zx); Cezez = -(eps_zi - eps_o) / eps_zi;  

    Chxey = 1/(j*w*dz*mu_xz); Chxez = 1/(j*w*dy*mu_xy); Chxhx = (mu_xi - mu_o) / mu_xi;  
    Chyez = 1/(j*w*dx*mu_yx); Chyex = 1/(j*w*dz*mu_yz); Chyhy = (mu_yi - mu_o) / mu_yi;  
    Chzex = 1/(j*w*dy*mu_zy); Chzey = 1/(j*w*dx*mu_zx); Chzhz = (mu_zi - mu_o) / mu_zi;  
    
    ! adjustments
    Cexhy(:,1,:) = (0.0,0.0); Cexhy(:,ny,:) = (0.0,0.0); Cexhy(:,:,1) = (0.0,0.0); Cexhy(:,:,nz) = (0.0,0.0); Cexhy(nx,:,:) = (0.0,0.0);
    Cexhz(:,1,:) = (0.0,0.0); Cexhz(:,ny,:) = (0.0,0.0); Cexhz(:,:,1) = (0.0,0.0); Cexhz(:,:,nz) = (0.0,0.0); Cexhz(nx,:,:) = (0.0,0.0);
    Cexex(:,1,:) = (0.0,0.0); Cexex(:,ny,:) = (0.0,0.0); Cexex(:,:,1) = (0.0,0.0); Cexex(:,:,nz) = (0.0,0.0); Cexex(nx,:,:) = (0.0,0.0);

    Ceyhz(1,:,:) = (0.0,0.0); Ceyhz(nx,:,:) = (0.0,0.0); Ceyhz(:,:,1) = (0.0,0.0); Ceyhz(:,:,nz) = (0.0,0.0); Ceyhz(:,ny,:) = (0.0,0.0);
    Ceyhx(1,:,:) = (0.0,0.0); Ceyhx(nx,:,:) = (0.0,0.0); Ceyhx(:,:,1) = (0.0,0.0); Ceyhx(:,:,nz) = (0.0,0.0); Ceyhx(:,ny,:) = (0.0,0.0);
    Ceyey(1,:,:) = (0.0,0.0); Ceyey(nx,:,:) = (0.0,0.0); Ceyey(:,:,1) = (0.0,0.0); Ceyey(:,:,nz) = (0.0,0.0); Ceyey(:,ny,:) = (0.0,0.0);

    Cezhx(1,:,:) = (0.0,0.0); Cezhx(nx,:,:) = (0.0,0.0); Cezhx(:,1,:) = (0.0,0.0); Cezhx(:,ny,:) = (0.0,0.0); Cezhx(:,:,nz) = (0.0,0.0);
    Cezhy(1,:,:) = (0.0,0.0); Cezhy(nx,:,:) = (0.0,0.0); Cezhy(:,1,:) = (0.0,0.0); Cezhy(:,ny,:) = (0.0,0.0); Cezhy(:,:,nz) = (0.0,0.0); 
    Cezez(1,:,:) = (0.0,0.0); Cezez(nx,:,:) = (0.0,0.0); Cezez(:,1,:) = (0.0,0.0); Cezez(:,ny,:) = (0.0,0.0); Cezez(:,:,nz) = (0.0,0.0);

    Chxey(:,ny,:) = (0.0,0.0); Chxey(:,:,nz) = (0.0,0.0); 
    Chxez(:,ny,:) = (0.0,0.0); Chxez(:,:,nz) = (0.0,0.0); 
    Chxhx(:,ny,:) = (0.0,0.0); Chxhx(:,:,nz) = (0.0,0.0); 

    Chyex(nx,:,:) = (0.0,0.0); Chyex(:,:,nz) = (0.0,0.0); 
    Chyez(nx,:,:) = (0.0,0.0); Chyez(:,:,nz) = (0.0,0.0); 
    Chyhy(nx,:,:) = (0.0,0.0); Chyhy(:,:,nz) = (0.0,0.0); 

    Chzex(nx,:,:) = (0.0,0.0); Chzex(:,ny,:) = (0.0,0.0); 
    Chzey(nx,:,:) = (0.0,0.0); Chzey(:,ny,:) = (0.0,0.0); 
    Chzhz(nx,:,:) = (0.0,0.0); Chzhz(:,ny,:) = (0.0,0.0); 
    
    Bex = Cexex * Einc_x; Bey = Ceyey * Einc_y; Bez = Cezez * Einc_z;
    Bhx = Chxhx * Hinc_x; Bhy = Chyhy * Hinc_y; Bhz = Chzhz * Hinc_z;
    
    Bex(1:nx,2:ny,2:nz) = Bex(1:nx,2:ny,2:nz) - &
                    (-Cexhz(1:nx,2:ny,2:nz)*Bhz(1:nx,2:ny,2:nz)+Cexhz(1:nx,2:ny,2:nz)*Bhz(1:nx,1:ny-1,2:nz) &
                     +Cexhy(1:nx,2:ny,2:nz)*Bhy(1:nx,2:ny,2:nz)-Cexhy(1:nx,2:ny,2:nz)*Bhy(1:nx,2:ny,1:nz-1));  

    Bey(2:nx,1:ny,2:nz) = Bey(2:nx,1:ny,2:nz) - &
                    (-Ceyhx(2:nx,1:ny,2:nz)*Bhx(2:nx,1:ny,2:nz)+Ceyhx(2:nx,1:ny,2:nz)*Bhx(2:nx,1:ny,1:nz-1) &
                     +Ceyhz(2:nx,1:ny,2:nz)*Bhz(2:nx,1:ny,2:nz)-Ceyhz(2:nx,1:ny,2:nz)*Bhz(1:nx-1,1:ny,2:nz));  
                     
    Bez(2:nx,2:ny,1:nz) = Bez(2:nx,2:ny,1:nz) - &
                    (-Cezhy(2:nx,2:ny,1:nz)*Bhy(2:nx,2:ny,1:nz)+Cezhy(2:nx,2:ny,1:nz)*Bhy(1:nx-1,2:ny,1:nz) &
                     +Cezhx(2:nx,2:ny,1:nz)*Bhx(2:nx,2:ny,1:nz)-Cezhx(2:nx,2:ny,1:nz)*Bhx(2:nx,1:ny-1,1:nz));  


    Bvec(0*Nxyz+1:1*Nxyz) = reshape(Bex, (/Nxyz/));
    Bvec(1*Nxyz+1:2*Nxyz) = reshape(Bey, (/Nxyz/));
    Bvec(2*Nxyz+1:3*Nxyz) = reshape(Bez, (/Nxyz/));


end subroutine setCoefficients

!===========================================================================================
!===========================================================================================
subroutine readProjectData 
	use global
	use ifport
!	use iflib
	use someFunctions
	implicit none

	character(80) string1, string2;
    integer narguments, ierr, parameterFound;
	
	narguments = nargs( ); 
	if (narguments .le. 1) then
		print *, 'Please enter the project file name';
		stop;
	endif
	

	geometry_filename = "no_geometry_file";
	initial_value_filename = "no_initial_value_file";
	result_filename = "no_results_file";
	materials_filename = "no_materials_file";
    run_on_gpu = .false.;
   	print_residual=.true.;

	call getarg(1,string1);	read(string1,*) projectFileName; 
	print *, 'Project name is ', projectFileName;
	open (1,FILE = projectFileName,IOSTAT=ierr, ACTION ='READ');
	if (ierr .gt. 0) then
			print 400, projectFileName;
400   FORMAT ('Project ', A, ' can not be opened');
	stop;
	end if
	read (1, *, iostat = ierr ) string1;
	do while (ierr .eq. 0 ) 
		read (1, *, iostat = ierr ) string2;
		parameterFound = assignParameters(string1, string2);
		read (1, *, iostat = ierr ) string1;
	end do

	close(1);
end subroutine readProjectData 
!===========================================================================================
! define and initialize the constants and arrays used in the program
subroutine initialize
	use global
	implicit none
	integer i;

	j=(0.0,1.0)
	pi = 3.1415926536;
	c = 2.99792458e+8; 
	eps_o = 8.854187817e-12;
	mu_o  = pi*4e-7;
	nu_o = sqrt(mu_o/eps_o);
	w = 2*pi*frequency;
	ko = w/c;
	Nx = floor(Sx/dx) + 1; Ny = floor(Sy/dy) + 1; Nz = floor(Sz/dz) + 1;
	Nx = (((Nx-1)/16)+1)*16;
	Ny = (((Ny-1)/16)+1)*16;
	
	Nxyz = Nx*Ny*Nz;
	Nu = 3*Nxyz;
	ell = 2;

    Nxm1=Nx-1; Nym1=Ny-1; Nzm1=Nz-1;
    Nxm2=Nx-2; Nym2=Ny-2; Nzm2=Nz-2;
    
	allocate(eps_r(Nxm1,Nym1,Nzm1), mu_r(Nxm1,Nym1,Nzm1));
	eps_r = eps_o;
	mu_r = mu_o;

	allocate(Escat_x(Nx,Ny,Nz), Escat_y(Nx,Ny,Nz), Escat_z(Nx,Ny,Nz));
	allocate(Hscat_x(Nx,Ny,Nz), Hscat_y(Nx,Ny,Nz), Hscat_z(Nx,Ny,Nz));
	allocate(Einc_x(Nx,Ny,Nz), Einc_y(Nx,Ny,Nz), Einc_z(Nx,Ny,Nz));
	allocate(Hinc_x(Nx,Ny,Nz), Hinc_y(Nx,Ny,Nz), Hinc_z(Nx,Ny,Nz));
	allocate(Medium(Nx,Ny,Nz));
	Escat_x  = (0.0,0.0); Escat_y  = (0.0,0.0); Escat_z  = (0.0,0.0);  
	Einc_x  = (0.0,0.0); Einc_y  = (0.0,0.0); Einc_z  = (0.0,0.0);  
	Hinc_x  = (0.0,0.0); Hinc_y  = (0.0,0.0); Hinc_z  = (0.0,0.0);  
	Medium = 1;
	
	allocate(Bvec(Nu), Xvec(Nu));
	allocate(work(Nu,3+2*(ell+1)));
	Xvec = (0.0,0.0);
	JMdistanceFromPml = 3;
	li = Npml + JMdistanceFromPml; lj=Npml + JMdistanceFromPml; lk=Npml + JMdistanceFromPml;
	ui = nx-(Npml + JMdistanceFromPml + 1); uj=ny-(Npml + JMdistanceFromPml + 1); uk=nz-(Npml + JMdistanceFromPml + 1);


	allocate(cjxyp (ui-li,1,uk-lk), cjxzp (ui-li,uj-lj,1), cjyxp (1,uj-lj,uk-lk));
	allocate(cjyzp (ui-li,uj-lj,1), cjzxp (1,uj-lj,uk-lk), cjzyp (ui-li,1,uk-lk));
	allocate(cjxym (ui-li,1,uk-lk), cjxzm (ui-li,uj-lj,1), cjyxm (1,uj-lj,uk-lk));
	allocate(cjyzm (ui-li,uj-lj,1), cjzxm (1,uj-lj,uk-lk), cjzym (ui-li,1,uk-lk));
	allocate(cmxyp (ui-li,1,uk-lk), cmxzp (ui-li,uj-lj,1), cmyxp (1,uj-lj,uk-lk));
	allocate(cmyzp (ui-li,uj-lj,1), cmzxp (1,uj-lj,uk-lk), cmzyp (ui-li,1,uk-lk));
	allocate(cmxym (ui-li,1,uk-lk), cmxzm (ui-li,uj-lj,1), cmyxm (1,uj-lj,uk-lk));
	allocate(cmyzm (ui-li,uj-lj,1), cmzxm (1,uj-lj,uk-lk), cmzym (ui-li,1,uk-lk));
	allocate(theta(361), farEtheta(361), farEphi(361),sigmaThetaTheta(361),sigmaPhiTheta(361)); 
	phi = 0.0; forall (i=1:361) theta(i) = (i-1)*pi/180;


	allocate(eps_xy(Nx,Ny,Nz), eps_xz(Nx,Ny,Nz), eps_xi(Nx,Ny,Nz));
	allocate(eps_yz(Nx,Ny,Nz), eps_yx(Nx,Ny,Nz), eps_yi(Nx,Ny,Nz));
	allocate(eps_zx(Nx,Ny,Nz), eps_zy(Nx,Ny,Nz), eps_zi(Nx,Ny,Nz));

	allocate(mu_xy(Nx,Ny,Nz), mu_xz(Nx,Ny,Nz), mu_xi(Nx,Ny,Nz));
	allocate(mu_yz(Nx,Ny,Nz), mu_yx(Nx,Ny,Nz), mu_yi(Nx,Ny,Nz));
	allocate(mu_zx(Nx,Ny,Nz), mu_zy(Nx,Ny,Nz), mu_zi(Nx,Ny,Nz));

	eps_xy = eps_o; eps_xz = eps_o; eps_xi = eps_o;
	eps_yz = eps_o; eps_yx = eps_o; eps_yi = eps_o;
	eps_zx = eps_o; eps_zy = eps_o; eps_zi = eps_o;

	mu_xy = mu_o; mu_xz = mu_o; mu_xi = mu_o;
	mu_yz = mu_o; mu_yx = mu_o; mu_yi = mu_o;
	mu_zx = mu_o; mu_zy = mu_o; mu_zi = mu_o;

end subroutine initialize

!==============================================================================
subroutine setIncidentField
	use global
	implicit none
	integer indx,indy,indz;

! traveling in the + z direction
	forall (indx=1:Nx,indy=1:Ny,indz=1:Nz) Einc_x(indx,indy,indz) = exp(-j*ko*indz*dz); 
	forall (indx=1:Nx,indy=1:Ny,indz=1:Nz) Hinc_y(indx,indy,indz) = exp(-j*ko*(indz+0.5)*dz); 
!	forall (indx=1:Nxm1,indy=2:Nym1,indz=2:Nzm1) Einc_x(indx,indy,indz) = exp(-j*ko*indz*dz); 
!	forall (indx=1:Nxm1,indy=1:Ny,indz=1:Nzm1) Hinc_y(indx,indy,indz) = exp(-j*ko*(indz+0.5)*dz); 
    Hinc_y = Hinc_y / nu_o;

end subroutine setIncidentField

!==============================================================================
subroutine createObjects
	use global
	implicit none
	integer r;
	integer mi,mj,mk,ind;
	real coord_x, coord_y, coord_z, dist;


! create cubes
	do r=1,numberOfCubes,1
		do mi=1,Nxm1,1
			coord_x = (mi-0.5) * dx;
			do mj=1,Nym1,1
				coord_y = (mj-0.5) * dy;
				do mk=1,Nzm1,1
					coord_z = (mk-0.5) * dy;
					if ((coord_x >= cube_min_x(r)).and.(coord_x <= cube_max_x(r)).and. &
					(coord_y >= cube_min_y(r)).and.(coord_y <= cube_max_y(r)).and.  &
					(coord_z >= cube_min_z(r)).and.(coord_z <= cube_max_z(r))) then
						eps_r(mi,mj,mk) = eps_o * (cube_eps_r_real(r) + j*cube_eps_r_imag(r));
						mu_r(mi,mj,mk)  = mu_o * (cube_mu_r_real(r) + j*cube_mu_r_imag(r));
					end if
				end do
			end do
		end do
	end do

! create spheres
	do r=1,numberOfSpheres,1
		do mi=1,Nxm1,1
			coord_x = (mi-0.5) * dx;
			do mj=1,Nym1,1
				coord_y = (mj-0.5) * dy;
				do mk=1,Nzm1,1
					coord_z = (mk-0.5) * dy;
					dist = sqrt( ((coord_x - centerX(r) )**2 ) + ( ( coord_y - centerY(r))**2 ) &
						 & + (( coord_z - centerZ(r))**2 ) ); 
					if (dist <= radius(r)) then
						eps_r(mi,mj,mk) = eps_o * (sphere_eps_r_real(r) + j*sphere_eps_r_imag(r));
						mu_r(mi,mj,mk)  = mu_o * (sphere_mu_r_real(r) + j*sphere_mu_r_imag(r));
					end if
				end do
			end do
		end do
	end do

! create cylinders
	do r=1,numberOfcylinders,1
		do mi=1,Nxm1,1
			coord_x = (mi-0.5) * dx;
			do mj=1,Nym1,1
				coord_y = (mj-0.5) * dy;
				do mk=1,Nzm1,1
					coord_z = (mk-0.5) * dy;
					dist = sqrt( ((coord_x - cylinder_center_x(r) )**2 ) + ( ( coord_y - cylinder_center_y(r))**2 ) ); 
					if ((dist <= cylinder_radius(r)).and.(coord_z >= cylinder_min_z(r)).and.(coord_z <= cylinder_max_z(r))) then
						eps_r(mi,mj,mk) = eps_o * (cylinder_eps_r_real(r) + j*cylinder_eps_r_imag(r));
						mu_r(mi,mj,mk)  = mu_o * (cylinder_mu_r_real(r) + j*cylinder_mu_r_imag(r));
					end if
				end do
			end do
		end do
	end do

	eps_xi(2:Nxm1, 2:Nym1, 2:Nzm1) = (eps_r(2:Nxm1, 2:Nym1, 2:Nzm1) + eps_r(2:Nxm1, 1:Nym2, 2:Nzm1) &
								+ eps_r(2:Nxm1, 1:Nym2, 1:Nzm2) + eps_r(2:Nxm1, 2:Nym1, 1:Nzm2))/4;
	eps_xy = eps_xi;
	eps_xz = eps_xi;


	eps_yi(2:Nxm1, 2:Nym1, 2:Nzm1) = (eps_r(2:Nxm1, 2:Nym1, 2:Nzm1) + eps_r(1:Nxm2, 2:Nym1, 2:Nzm1) &
								+ eps_r(2:Nxm1, 2:Nym1, 1:Nzm2) + eps_r(1:Nxm2, 2:Nym1, 1:Nzm2))/4;
	eps_yz = eps_yi;
	eps_yx = eps_yi;

	eps_zi(2:Nxm1, 2:Nym1, 2:Nzm1) = (eps_r(2:Nxm1, 2:Nym1, 2:Nzm1) + eps_r(1:Nxm2, 2:Nym1, 2:Nzm1) &
								+ eps_r(2:Nxm1, 1:Nym2, 2:Nzm1) + eps_r(1:Nxm2, 1:Nym2, 2:Nzm1))/4;
	eps_zx = eps_zi;
	eps_zy = eps_zi;

	mu_xi(2:Nxm1, 2:Nym1, 2:Nzm1) = (mu_r(2:Nxm1, 2:Nym1, 2:Nzm1) + mu_r(1:Nxm2, 2:Nym1, 2:Nzm1))/2;
	mu_xy = mu_xi; 
	mu_xz = mu_xi; 

	mu_yi(2:Nxm1, 2:Nym1, 2:Nzm1) = (mu_r(2:Nxm1, 2:Nym1, 2:Nzm1) + mu_r(2:Nxm1, 1:Nym2, 2:Nzm1))/2;
	mu_yz = mu_yi; 
	mu_yx = mu_yi; 

	mu_zi(2:Nxm1, 2:Nym1, 2:Nzm1) = (mu_r(2:Nxm1, 2:Nym1, 2:Nzm1) + mu_r(2:Nxm1, 2:Nym1, 1:Nzm2))/2;
	mu_zx = mu_zi; 
	mu_zy = mu_zi; 

	deallocate(eps_r, mu_r);

end subroutine createObjects


!==============================================================================
subroutine solveSparseMatrices
	use global
	use ifport
	implicit none
	integer ind, mi, mj, mk, Nzero

	integer       ldw;
	integer       ldrw, ierr;
	complex*16, allocatable :: rwork(:,:);
	integer,	allocatable :: iwork(:);

	real*8 droptol

    REAL(8) elapsed_time
	external      matvec, precond;
    elapsed_time = TIMEF( )

	nonzero_x=.false.;
	Nxyz = Nx * Ny * Nz;

	ldw = 3*Nx*Ny*Nz*(3+2*(ell+1));
	!ldrw = ell+1;
	ldrw = (ell+1)*(3+2*(ell+1));

	if (printstep.eq.0) printstep=4;

    if (run_on_gpu) then
		print *, 'This version of FDFD can run only on CPU!';
        stop;
	else
    	call zbcg2 (print_residual,ell,Nu,Xvec,nonzero_x,Bvec, matvec,precond,tolerance, &
    	                  & number_of_iterations,work,info)   	
   	end if

	print *, 'iterative solve result ', info;
	if (info>0) stop;

    elapsed_time = TIMEF( )
    PRINT *, 'Program has used', elapsed_time, ' seconds'

    time_filename = 'sim_time_' // trim(projectFileName(1:scan(projectFileName,'.')-1)) // '.m';

    open (1, FILE = time_filename);
     write (1,*) '%nx =', nx, ', ny=', ny, ', nz=', nz
     write (1,*) '%Total simulation time is =', elapsed_time, ' seconds'
     write (1,*) '%Total simulation time is =', elapsed_time/60, ' minutes'
     write (1,*) 'simulation_time =', elapsed_time, ';'
!     write (1,*) 'Number of iterations is =', elapsed_time
    close (1)
    
    Hscat_x = Bhx - tmpx;
    Hscat_y = Bhy - tmpy;
    Hscat_z = Bhz - tmpz;


    Escat_x = reshape(Xvec(0*Nxyz+1:1*Nxyz), (/Nx,Ny,Nz/));
    Escat_y = reshape(Xvec(1*Nxyz+1:2*Nxyz), (/Nx,Ny,Nz/));
    Escat_z = reshape(Xvec(2*Nxyz+1:3*Nxyz), (/Nx,Ny,Nz/));


end subroutine solveSparseMatrices
!==============================================================================
subroutine createResultFile
	use global
	implicit none
	integer mi,mj,mk;

	if (result_filename.eq."no_results_file") then
		return;
	endif


	open (1, FILE = result_filename);
	write(1,*) 'Nx'
	write(1,*)  Nx
	write(1,*) 'Ny'
	write(1,*)  Ny
	write(1,*) 'Nz'
	write(1,*)  Nz
	write(1,*) 'Escat_x Data'

	do mi = 1,Nx,1
		do mj = 1,Ny,1
			do mk = 1,Nz,1
				write(1,210) real(Escat_x(mi,mj,mk)), aimag(Escat_x(mi,mj,mk));
			end do
		end do
	end do 

	write(1,*) 'Escat_y Data'
	do mi = 1,Nx,1
		do mj = 1,Ny,1
			do mk = 1,Nz,1
				write(1,210) real(Escat_y(mi,mj,mk)), aimag(Escat_y(mi,mj,mk));
			end do
		end do
	end do 
	write(1,*) 'Escat_z Data'
	do mi = 1,Nx,1
		do mj = 1,Ny,1
			do mk = 1,Nz,1
				write(1,210) real(Escat_z(mi,mj,mk)), aimag(Escat_z(mi,mj,mk));
			end do
		end do
	end do 

	close(1);
	210   FORMAT (E15.8,' ',E15.8);

end subroutine createResultFile
!==============================================================================

!==============================================================================
subroutine matvec(n,x,y)
	use global
	implicit none
    integer n; complex*16 x(nx,ny,nz,3), y(nx,ny,nz,3);

    tmpx(:,1:Nym1,1:Nzm1) = Chxez(:,1:Nym1,1:Nzm1)*(x(:,2:Ny,1:Nzm1,3)-x(:,1:Nym1,1:Nzm1,3)) &
            + Chxey(:,1:Nym1,1:Nzm1)*(-x(:,1:Nym1,2:Nz,2)+x(:,1:Nym1,1:Nzm1,2));

    tmpy(1:Nxm1,:,1:Nzm1) = Chyex(1:Nxm1,:,1:Nzm1)*(x(1:Nxm1,:,2:Nz,1)-x(1:Nxm1,:,1:Nzm1,1)) &
            + Chyez(1:Nxm1,:,1:Nzm1)*(-x(2:Nx,:,1:Nzm1,3)+x(1:Nxm1,:,1:Nzm1,3));

    tmpz(1:Nxm1,1:Nym1,:) = Chzey(1:Nxm1,1:Nym1,:)*(x(2:Nx,1:Nym1,:,2)-x(1:Nxm1,1:Nym1,:,2)) &
            + Chzex(1:Nxm1,1:Nym1,:)*(-x(1:Nxm1,2:Ny,:,1)+x(1:Nxm1,1:Nym1,:,1));

    y(1:Nx,2:Ny,2:Nz,1) = x(1:Nx,2:Ny,2:Nz,1) &
            - (Cexhz(1:Nx,2:Ny,2:Nz)*(-tmpz(1:Nx,2:Ny,2:Nz )+tmpz(1:Nx,1:Nym1,2:Nz)) &
            +  Cexhy(1:Nx,2:Ny,2:Nz)*(tmpy(1:Nx,2:Ny,2:Nz)-tmpy(1:Nx,2:Ny,1:Nzm1)));

    y(2:Nx,1:Ny,2:Nz,2) = x(2:Nx,1:Ny,2:Nz,2) &
            - (Ceyhx(2:Nx,1:Ny,2:Nz)*(-tmpx(2:Nx,1:Ny,2:Nz)+tmpx(2:Nx,1:Ny,1:Nzm1)) &
            +  Ceyhz(2:Nx,1:Ny,2:Nz)*(tmpz(2:Nx,1:Ny,2:Nz)-tmpz(1:Nxm1,1:Ny,2:Nz)));

    y(2:Nx,2:Ny,1:Nz,3) = x(2:Nx,2:Ny,1:Nz,3) &
            - (Cezhy(2:Nx,2:Ny,1:Nz)*(-tmpy(2:Nx,2:Ny,1:Nz)+tmpy(1:Nxm1,2:Ny,1:Nz)) &
            +  Cezhx(2:Nx,2:Ny,1:Nz)*(tmpx(2:Nx,2:Ny,1:Nz)-tmpx(2:Nx,1:Nym1,1:Nz)));

end subroutine matvec
!==============================================================================
subroutine calculateJandM

		use global;
		implicit none
		integer Nx1, Ny1, Nz1, mi, mj, mk;

		Nx1 = Nx-1; Ny1 = Ny-1; Nz1 = Nz-1;


		cmyxp(1,:,:) =  0.5 * (Escat_z (ui,lj:uj-1,lk:uk-1) + Escat_z (ui,lj+1:uj,lk:uk-1));
		cmzxp(1,:,:) = -0.5 * (Escat_y (ui,lj:uj-1,lk:uk-1) + Escat_y (ui,lj:uj-1,lk+1:uk));
		cmxyp(:,1,:) = -0.5 * (Escat_z (li:ui-1,uj,lk:uk-1) + Escat_z (li+1:ui,uj,lk:uk-1));
		cmzyp(:,1,:) =  0.5 * (Escat_x (li:ui-1,uj,lk:uk-1) + Escat_x (li:ui-1,uj,lk+1:uk));
		cmxzp(:,:,1) =  0.5 * (Escat_y (li:ui-1,lj:uj-1,uk) + Escat_y (li+1:ui,lj:uj-1,uk));
		cmyzp(:,:,1) = -0.5 * (Escat_x (li:ui-1,lj:uj-1,uk) + Escat_x (li:ui-1,lj+1:uj,uk));
		cmyxm(1,:,:) = -0.5 * (Escat_z (li,lj:uj-1,lk:uk-1) + Escat_z (li,lj+1:uj,lk:uk-1));
		cmzxm(1,:,:) =  0.5 * (Escat_y (li,lj:uj-1,lk:uk-1) + Escat_y (li,lj:uj-1,lk+1:uk));
		cmxym(:,1,:) =  0.5 * (Escat_z (li:ui-1,lj,lk:uk-1) + Escat_z (li+1:ui,lj,lk:uk-1));
		cmzym(:,1,:) = -0.5 * (Escat_x (li:ui-1,lj,lk:uk-1) + Escat_x (li:ui-1,lj,lk+1:uk));
		cmxzm(:,:,1) = -0.5 * (Escat_y (li:ui-1,lj:uj-1,lk) + Escat_y (li+1:ui,lj:uj-1,lk));
		cmyzm(:,:,1) =  0.5 * (Escat_x (li:ui-1,lj:uj-1,lk) + Escat_x (li:ui-1,lj+1:uj,lk));

		cjyxp(1,:,:) = (Hscat_z (ui,lj:uj-1,lk:uk-1) + Hscat_z (ui,lj:uj-1,lk+1:uk));
		cjyxp(1,:,:) = -0.25 * (cjyxp(1,:,:) + Hscat_z (ui-1,lj:uj-1,lk:uk-1) + Hscat_z (ui-1,lj:uj-1,lk+1:uk));
		cjzxp(1,:,:) = (Hscat_y (ui,lj:uj-1,lk:uk-1) + Hscat_y (ui,lj+1:uj,lk:uk-1));
		cjzxp(1,:,:) =  0.25 * (cjzxp(1,:,:) + Hscat_y (ui-1,lj:uj-1,lk:uk-1) + Hscat_y (ui-1,lj+1:uj,lk:uk-1));
		cjzyp(:,1,:) = (Hscat_x (li:ui-1,uj,lk:uk-1) + Hscat_x (li+1:ui,uj,lk:uk-1));
		cjzyp(:,1,:) = -0.25 * (cjzyp(:,1,:) + Hscat_x (li:ui-1,uj-1,lk:uk-1) + Hscat_x (li+1:ui,uj-1,lk:uk-1));
		cjxyp(:,1,:) = (Hscat_z (li:ui-1,uj,lk:uk-1) + Hscat_z (li:ui-1,uj,lk+1:uk));
		cjxyp(:,1,:) = 0.25 * (cjxyp(:,1,:) + Hscat_z (li:ui-1,uj-1,lk:uk-1) + Hscat_z (li:ui-1,uj-1,lk+1:uk));
		cjyzp(:,:,1) = (Hscat_x (li:ui-1,lj:uj-1,uk) + Hscat_x (li+1:ui,lj:uj-1,uk));
		cjyzp(:,:,1) = 0.25 * (cjyzp(:,:,1) + Hscat_x (li:ui-1,lj:uj-1,uk-1) + Hscat_x (li+1:ui,lj:uj-1,uk-1));
		cjxzp(:,:,1) = (Hscat_y (li:ui-1,lj:uj-1,uk) + Hscat_y (li:ui-1,lj+1:uj,uk));
		cjxzp(:,:,1) = -0.25 * (cjxzp(:,:,1) + Hscat_y (li:ui-1,lj:uj-1,uk-1) + Hscat_y (li:ui-1,lj+1:uj,uk-1));
		cjyxm(1,:,:) = (Hscat_z (li,lj:uj-1,lk:uk-1) + Hscat_z (li,lj:uj-1,lk+1:uk));
		cjyxm(1,:,:) =  0.25 * (cjyxm(1,:,:) + Hscat_z (li-1,lj:uj-1,lk:uk-1) + Hscat_z (li-1,lj:uj-1,lk+1:uk));
		cjzxm(1,:,:) = (Hscat_y (li,lj:uj-1,lk:uk-1) + Hscat_y (li,lj+1:uj,lk:uk-1));
		cjzxm(1,:,:) =  -0.25 * (cjzxm(1,:,:) + Hscat_y (li-1,lj:uj-1,lk:uk-1) + Hscat_y (li-1,lj+1:uj,lk:uk-1));
		cjzym(:,1,:) = (Hscat_x (li:ui-1,lj,lk:uk-1) + Hscat_x (li+1:ui,lj,lk:uk-1));
		cjzym(:,1,:) =  0.25 * (cjzym(:,1,:) + Hscat_x (li:ui-1,lj-1,lk:uk-1) + Hscat_x (li+1:ui,lj-1,lk:uk-1));
		cjxym(:,1,:) = (Hscat_z (li:ui-1,lj,lk:uk-1) + Hscat_z (li:ui-1,lj,lk+1:uk));
		cjxym(:,1,:) = -0.25 * (cjxym(:,1,:) + Hscat_z (li:ui-1,lj-1,lk:uk-1) + Hscat_z (li:ui-1,lj-1,lk+1:uk));
		cjyzm(:,:,1) = (Hscat_x (li:ui-1,lj:uj-1,lk) + Hscat_x (li+1:ui,lj:uj-1,lk));
		cjyzm(:,:,1) = -0.25 * (cjyzm(:,:,1) + Hscat_x (li:ui-1,lj:uj-1,lk-1) + Hscat_x (li+1:ui,lj:uj-1,lk-1));
		cjxzm(:,:,1) = (Hscat_y (li:ui-1,lj:uj-1,lk) + Hscat_y (li:ui-1,lj+1:uj,lk));
		cjxzm(:,:,1) =  0.25 * (cjxzm(:,:,1) + Hscat_y (li:ui-1,lj:uj-1,lk-1) + Hscat_y (li:ui-1,lj+1:uj,lk-1));


end subroutine calculateJandM 
!==============================================================================
!****************************************************************************
! calculate far fields
subroutine calculatefarfields
	use global
	implicit none
	complex cff1,cff2;
	complex, dimension(361) :: Ntheta, Ltheta, Nphi, Lphi, rpr;
	integer ni,nj,nk,nf,i;
	real k,distance;

	j=(0.0,1.0);
	k = w*(mu_o*eps_o)**0.5;

	Ntheta = (0.0,0.0); Ltheta = (0.0,0.0); Nphi = (0.0,0.0); Lphi = (0.0,0.0); rpr = (0.0,0.0);
			! for +ax direction
			do nj =lj,uj-1,1
			do nk =lk,uk-1,1
			rpr = (ui - ci)*dx*sin(theta)*cos(phi) &
				  &     + (nj-cj+0.5)*dy*sin(theta)*sin(phi) &
				  &     + (nk-ck+0.5)*dz*cos(theta);
			Ntheta = Ntheta + dy*dz*(cjyxp(1,nj-lj+1,nk-lk+1)*cos(theta)*sin(phi) &
					 &  - cjzxp(1,nj-lj+1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Ltheta = Ltheta + dy*dz*(cmyxp(1,nj-lj+1,nk-lk+1)*cos(theta)*sin(phi) &
					 &  - cmzxp(1,nj-lj+1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Nphi = Nphi + dy*dz*(cjyxp(1,nj-lj+1,nk-lk+1)*cos(phi))*exp(j*k*rpr);
			Lphi = Lphi + dy*dz*(cmyxp(1,nj-lj+1,nk-lk+1)*cos(phi))*exp(j*k*rpr);
				end do
			end do
 

			! for -ax direction
			do nj =lj,uj-1,1
			do nk =lk,uk-1,1
			rpr = (li - ci)*dx*sin(theta)*cos(phi) &
				   &     + (nj-cj+0.5)*dy*sin(theta)*sin(phi) &
				   &     + (nk-ck+0.5)*dz*cos(theta);
			Ntheta = Ntheta + dy*dz*(cjyxm(1,nj-lj+1,nk-lk+1)*cos(theta)*sin(phi) &
					 &	 - cjzxm(1,nj-lj+1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Ltheta = Ltheta + dy*dz*(cmyxm(1,nj-lj+1,nk-lk+1)*cos(theta)*sin(phi) &
					 &   - cmzxm(1,nj-lj+1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Nphi = Nphi + dy*dz*(cjyxm(1,nj-lj+1,nk-lk+1)*cos(phi))*exp(j*k*rpr);
			Lphi = Lphi + dy*dz*(cmyxm(1,nj-lj+1,nk-lk+1)*cos(phi))*exp(j*k*rpr);
				end do
			end do

			! for +ay direction
			do ni =li,ui-1,1
			do nk =lk,uk-1,1
			rpr = (ni - ci + 0.5)*dx*sin(theta)*cos(phi) &
				  &      + (uj-cj)*dy*sin(theta)*sin(phi) &
				  &      + (nk-ck+0.5)*dz*cos(theta);
			Ntheta = Ntheta + dx*dz*(cjxyp(ni-li+1,1,nk-lk+1)*cos(theta)*cos(phi) &
					&    - cjzyp(ni-li+1,1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Ltheta = Ltheta + dx*dz*(cmxyp(ni-li+1,1,nk-lk+1)*cos(theta)*cos(phi) &
					&    - cmzyp(ni-li+1,1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Nphi = Nphi + dx*dz*(-cjxyp(ni-li+1,1,nk-lk+1)*sin(phi))*exp(j*k*rpr);
			Lphi = Lphi + dx*dz*(-cmxyp(ni-li+1,1,nk-lk+1)*sin(phi))*exp(j*k*rpr);
				end do
			end do

			! for -ay direction
			do ni =li,ui-1,1
			do nk =lk,uk-1,1
			rpr = (ni - ci + 0.5)*dx*sin(theta)*cos(phi) &
				  &      + (lj-cj)*dy*sin(theta)*sin(phi) &
				  &      + (nk-ck+0.5)*dz*cos(theta);
			Ntheta = Ntheta + dx*dz*(cjxym(ni-li+1,1,nk-lk+1)*cos(theta)*cos(phi) &
				  &	     - cjzym(ni-li+1,1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Ltheta = Ltheta + dx*dz*(cmxym(ni-li+1,1,nk-lk+1)*cos(theta)*cos(phi) &
				  &      - cmzym(ni-li+1,1,nk-lk+1)*sin(theta))*exp(j*k*rpr);
			Nphi = Nphi + dx*dz*(-cjxym(ni-li+1,1,nk-lk+1)*sin(phi))*exp(j*k*rpr);
			Lphi = Lphi + dx*dz*(-cmxym(ni-li+1,1,nk-lk+1)*sin(phi))*exp(j*k*rpr);
				end do
			end do


			! for +az direction
			do ni =li,ui-1,1
			do nj =lj,uj-1,1
			rpr = (ni-ci+0.5)*dx*sin(theta)*cos(phi) &
				  &      + (nj - cj + 0.5)*dy*sin(theta)*sin(phi) &
				  &      + (uk-ck)*dz*cos(theta);
			Ntheta = Ntheta + dx*dy*(cjxzp(ni-li+1,nj-lj+1,1)*cos(theta)*cos(phi) &
				  &      + cjyzp(ni-li+1,nj-lj+1,1)*cos(theta)*sin(phi))*exp(j*k*rpr);
			Ltheta = Ltheta + dx*dy*(cmxzp(ni-li+1,nj-lj+1,1)*cos(theta)*cos(phi) &
				  &      + cmyzp(ni-li+1,nj-lj+1,1)*cos(theta)*sin(phi))*exp(j*k*rpr);
			Nphi = Nphi + dx*dy*(-cjxzp(ni-li+1,nj-lj+1,1)*sin(phi)+cjyzp(ni-li+1,nj-lj+1,1)*cos(phi))*exp(j*k*rpr);
			Lphi = Lphi + dx*dy*(-cmxzp(ni-li+1,nj-lj+1,1)*sin(phi)+cmyzp(ni-li+1,nj-lj+1,1)*cos(phi))*exp(j*k*rpr);
				end do
			end do

			! for -az direction
			do ni =li,ui-1,1
			do nj =lj,uj-1,1
			rpr = (ni-ci+0.5)*dx*sin(theta)*cos(phi) &
				   &     + (nj - cj + 0.5)*dy*sin(theta)*sin(phi) &
				   &     + (lk-ck)*dz*cos(theta);
			Ntheta = Ntheta + dx*dy*(cjxzm(ni-li+1,nj-lj+1,1)*cos(theta)*cos(phi) &
				   &    + cjyzm(ni-li+1,nj-lj+1,1)*cos(theta)*sin(phi))*exp(j*k*rpr);
			Ltheta = Ltheta + dx*dy*(cmxzm(ni-li+1,nj-lj+1,1)*cos(theta)*cos(phi) &
				   &    + cmyzm(ni-li+1,nj-lj+1,1)*cos(theta)*sin(phi))*exp(j*k*rpr);
			Nphi = Nphi + dx*dy*(-cjxzm(ni-li+1,nj-lj+1,1)*sin(phi)+cjyzm(ni-li+1,nj-lj+1,1)*cos(phi))*exp(j*k*rpr);
			Lphi = Lphi + dx*dy*(-cmxzm(ni-li+1,nj-lj+1,1)*sin(phi)+cmyzm(ni-li+1,nj-lj+1,1)*cos(phi))*exp(j*k*rpr);
				end do
			end do

		distance = 1000*c/frequency;
!		farEtheta = abs((-j*k*exp(-j*k*distance))*(1/(4*pi*distance))*(Lphi+nu_o*Ntheta));
!		farEphi   = abs((j*k*exp(-j*k*distance))*(1/(4*pi*distance))*(Ltheta-nu_o*Nphi));
		sigmaThetaTheta = 4*pi*((distance)**2)*(abs((-j*k*exp(-j*k*distance))*(1/(4*pi*distance))*(Lphi+nu_o*Ntheta))**2);
		sigmaPhiTheta = 4*pi*((distance)**2)*(abs((j*k*exp(-j*k*distance))*(1/(4*pi*distance))*(Ltheta-nu_o*Nphi))**2);

end subroutine calculatefarfields
!****************************************************************************
!****************************************************************************
! write farfield results to a file
subroutine saveFarFieldsToFile
	use global
	implicit none
	integer i;

	result_filename = 'farfield_theta_' // trim(projectFileName(1:scan(projectFileName,'.')-1)) // '.m';

	open (1, FILE = result_filename);
	write(1,'("stt = [ ")') 
	do i=1,360,1 
		write(1,11) i-1 , sigmaThetaTheta(i)
	end do
	write(1,'(" ]; ")') 
	write(1,'("plot(stt(1:181,1),stt(1:181,2),''r--'',''linewidth'',1.5 ); ")') 
	write(1,'("grid on ")') 
	write(1,'("legend(''FDFD'');")') 
	write(1,'("xlabel(''\theta [degrees]'',''FontSize'',12);")') 
	write(1,'("ylabel(''\sigma_{\theta\theta}'',''FontSize'',14);")') 



	close(1);

	result_filename = 'farfield_phi_' // trim(projectFileName(1:scan(projectFileName,'.')-1)) // '.m';

	open (1, FILE = result_filename);
	write(1,'("spt = [ ")') 
	do i=1,360,1 
		write(1,11) i-1 , sigmaPhiTheta(i)
	end do
	write(1,'(" ]; ")') 
	write(1,'("plot(spt(1:181,1),spt(1:181,2),''r--'',''linewidth'',1.5 ); ")') 
	write(1,'("grid on ")') 
	write(1,'("legend(''FDFD'');")') 
	write(1,'("xlabel(''\theta [degrees]'',''FontSize'',12);")') 
	write(1,'("ylabel(''\sigma_{\phi\theta}'',''FontSize'',14);")') 

	close(1);
11  FORMAT(I5, ' ', E15.8)

end subroutine saveFarFieldsToFile

!=============================================================================
subroutine saveJandMtoFile
	use global
	implicit none
	integer mi, mj;

	open (1, FILE = "cjxzp.txt");
	write(1,*) ui-li;
	write(1,*) uj-lj;
	do mi=1,ui-li,1 
		do mj=1,uj-lj,1 
			write(1,*) abs(cjxzp(mi,mj,1))
		end do
	end do
	close(1);

	open (1, FILE = "cjyzp.txt");
	write(1,*) ui-li;
	write(1,*) uj-lj;
	do mi=1,ui-li,1 
		do mj=1,uj-lj,1 
			write(1,*) abs(cjyzp(mi,mj,1))
		end do
	end do
	close(1);

	open (1, FILE = "cmxzp.txt");
	write(1,*) ui-li;
	write(1,*) uj-lj;
	do mi=1,ui-li,1 
		do mj=1,uj-lj,1 
			write(1,*) abs(cmxzp(mi,mj,1))
		end do
	end do
	close(1);

	open (1, FILE = "cmyzp.txt");
	write(1,*) ui-li;
	write(1,*) uj-lj;
	do mi=1,ui-li,1 
		do mj=1,uj-lj,1 
			write(1,*) abs(cmyzp(mi,mj,1))
		end do
	end do
	close(1);


end subroutine saveJandMtoFile
!****************************************************************************
subroutine precond(n,rhs)
		use global
		implicit none
		integer n; complex*16 rhs(*);

!		call lusol(n, rhs, rhs, P_mat, JCOL, IROW, Pmat_diag_ptr)

		return;

end subroutine precond
!==============================================================================

program fdfd_3d
	use global
	implicit none

	call readProjectData
	call initialize

	call setIncidentField
	call createObjects
	call setPMLboundaries
	call setCoefficients
	call solveSparseMatrices
	call calculateJandM
	call calculatefarfields
	call createResultFile
	call saveFarFieldsToFile

end program fdfd_3d

