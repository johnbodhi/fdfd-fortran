subroutine saveCapturedFieldsToFile
	use global
	implicit none
	integer mi,mj,mk; 
    

	result_filename = 'captured_fields_' // trim(projectFileName(1:scan(projectFileName,'.')-1)) // '.m';

	open (1, FILE = result_filename);
    
    write(1,'("clc; clear all; close all;")') 
    
    write(1,*) 'Nx = ', Nx
	write(1,*) 'Ny = ', Ny
	write(1,*) 'Nz = ', Nz
	write(1,*) 'Nxy = ', Nxy
	write(1,*) 'ko = ', ko
	write(1,*) 'nu_o = ', nu_o
    write(1,*) 'dt = ', dt
    
	write(1,'("Einc_x = [ ")')   
    do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Einc_x(mi,mj,rp-1)), aimag(Einc_x(mi,mj,rp-1));
		end do
    end do 
    write(1,'(" ]; ")') 
    
    write(1,'("Einc_y = [ ")') 
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Einc_y(mi,mj,rp-1)), aimag(Einc_y(mi,mj,rp-1));
		end do
    end do    
    write(1,'(" ]; ")')
    
    write(1,'("Escat_x_Ref = [ ")') 
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Escat_x(mi,mj,rp)), aimag(Escat_x(mi,mj,rp));
		end do
    end do 
    write(1,'(" ]; ")')
       
	write(1,'("Escat_y_Ref = [ ")')
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Escat_y(mi,mj,rp)), aimag(Escat_y(mi,mj,rp));
		end do
    end do 	
    write(1,'(" ]; ")')
    
     write(1,'("Hscat_x_Ref = [ ")') 
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Hscat_x(mi,mj,rp)), aimag(Hscat_x(mi,mj,rp));
		end do
    end do 
    write(1,'(" ]; ")')
       
	write(1,'("Hscat_y_Ref = [ ")')
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Hscat_y(mi,mj,rp)), aimag(Hscat_y(mi,mj,rp));
		end do
    end do 	
    write(1,'(" ]; ")')
    
    write(1,'("Escat_x_Tra = [ ")')
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Escat_x(mi,mj,tp)), aimag(Escat_x(mi,mj,tp));
		end do
    end do 
    write(1,'(" ]; ")')
    
	write(1,'("Escat_y_Tra = [ ")')
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Escat_y(mi,mj,tp)), aimag(Escat_y(mi,mj,tp));
		end do
    end do   
    write(1,'(" ]; ")')
    
     write(1,'("Hscat_x_Tra = [ ")')
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Hscat_x(mi,mj,tp)), aimag(Hscat_x(mi,mj,tp));
		end do
    end do 
    write(1,'(" ]; ")')
    
	write(1,'("Hscat_y_Tra = [ ")')
	do mi = 1,Nx,1
		do mj = 1,Ny,1
            write(1,*) real(Hscat_y(mi,mj,tp)), aimag(Hscat_y(mi,mj,tp));
		end do
    end do   
    write(1,'(" ]; ")')    
    
    
    write(1,'("[ Gmco, Gmcr, Gpco, Gpcr, Tmco, Tmcr, Tpco, Tpcr ] = calculateGandT(  Einc_x, Einc_y, Escat_x_Ref, Escat_y_Ref, Hscat_x_Ref, Hscat_y_Ref, Escat_x_Tra, Escat_y_Tra, Hscat_x_Tra, Hscat_y_Tra, Nxy, ko, nu_o, dt );")') 

    close(1);
210 FORMAT(E15.8)   
end subroutine saveCapturedFieldsToFile