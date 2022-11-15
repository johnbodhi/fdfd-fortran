module global

	implicit none
	complex*16 j;
	character*500 projectFileName , initial_value_filename;
	character*500 result_filename, geometry_filename, time_filename;
	character*500 materials_filename;

	real Sx, Sy, Sz, dx, dy, dz, frequency, w, A, sp, rp, tp, ref, tra, dt;
    real, allocatable :: se(:,:,:), sm(:,:,:);
	complex*16,allocatable :: eps_r(:,:,:), mu_r(:,:,:);
	integer numberOfSpheres,numberOfMaterials,numberOfCubes,numberOfCylinders;

	real,dimension(1000) :: electric_conductivity, magnetic_conductivity;
    
	! sphere parameters
	real,dimension(1000) :: centerX, centerY, centerZ, radius
	real,dimension(1000) :: sphere_eps_r_real,sphere_mu_r_real, sphere_kappa_real;
	real,dimension(1000) :: sphere_eps_r_imag,sphere_mu_r_imag, sphere_kappa_imag;

	! cube parameters
	real,dimension(1000) :: cube_min_x, cube_min_y, cube_min_z
	real,dimension(1000) :: cube_max_x, cube_max_y, cube_max_z
	real,dimension(1000) :: cube_eps_r_real,cube_mu_r_real, cube_kappa_real;
	real,dimension(1000) :: cube_eps_r_imag,cube_mu_r_imag, cube_kappa_imag;

	! cylinder parameters
	real,dimension(1000) :: cylinder_center_x, cylinder_center_y, cylinder_min_z, cylinder_max_z,cylinder_radius
	real,dimension(1000) :: cylinder_eps_r_real,cylinder_mu_r_real, cylinder_kappa_real;
	real,dimension(1000) :: cylinder_eps_r_imag,cylinder_mu_r_imag, cylinder_kappa_imag;


	integer Nx, Ny, Nz, Nxyz, Nxy, Nu, Npml;
	integer Nxm1, Nym1, Nzm1
	integer Nxm2, Nym2, Nzm2
	real :: pi,c,eps_o,mu_o,nu_o, ko, phi;
	
	real,allocatable :: theta(:), farEtheta(:), farEphi(:),sigmaThetaTheta(:),sigmaPhiTheta(:); 
    
 
	integer, allocatable   :: Medium(:,:,:); 
	complex*16,allocatable :: Escat_x  (:,:,:), Escat_y  (:,:,:), Escat_z  (:,:,:);
	complex*16,allocatable :: Hscat_x  (:,:,:), Hscat_y  (:,:,:), Hscat_z  (:,:,:);
	complex*16,allocatable :: Einc_x  (:,:,:), Einc_y  (:,:,:), Einc_z  (:,:,:);
	complex*16,allocatable :: Hinc_x  (:,:,:), Hinc_y  (:,:,:), Hinc_z  (:,:,:);

	complex*16,allocatable :: Bex(:,:,:), Bey(:,:,:), Bez(:,:,:);
	complex*16,allocatable :: Bhx(:,:,:), Bhy(:,:,:), Bhz(:,:,:);
	complex*16,allocatable :: tmpx(:,:,:), tmpy(:,:,:), tmpz(:,:,:);

	complex*16,allocatable :: Cexhy  (:,:,:), Cexhz(:,:,:),Cexex(:,:,:);
	complex*16,allocatable :: Ceyhz  (:,:,:), Ceyhx(:,:,:),Ceyey(:,:,:);
	complex*16,allocatable :: Cezhx  (:,:,:), Cezhy(:,:,:),Cezez(:,:,:);

	complex*16,allocatable :: Chxey  (:,:,:), Chxez(:,:,:),Chxhx(:,:,:);
	complex*16,allocatable :: Chyez  (:,:,:), Chyex(:,:,:),Chyhy(:,:,:);
	complex*16,allocatable :: Chzex  (:,:,:), Chzey(:,:,:),Chzhz(:,:,:);
    
    complex*16,allocatable :: pex(:,:,:),pey(:,:,:),pez(:,:,:),phx(:,:,:),phy(:,:,:),phz(:,:,:);

	complex*16,allocatable :: Bvec(:), Xvec(:);

	complex*16,allocatable :: work(:,:);
	integer					  printstep, number_of_iterations;
	real*8                    tolerance;
	integer			    	  ell

	logical     print_residual
	logical     nonzero_x;
	logical     run_on_gpu;
	integer     info;
	

	integer li,lj,lk,ui,uj,uk,ci,cj,ck, JMdistanceFromPml;
    ! Arrays used in near to farfield transformation
	complex,allocatable :: cjxyp(:,:,:), cjxzp(:,:,:), cjyxp(:,:,:);
	complex,allocatable :: cjyzp(:,:,:), cjzxp(:,:,:), cjzyp(:,:,:);
	complex,allocatable :: cjxym(:,:,:), cjxzm(:,:,:), cjyxm(:,:,:);
	complex,allocatable :: cjyzm(:,:,:), cjzxm(:,:,:), cjzym(:,:,:);

	complex,allocatable :: cmxyp(:,:,:), cmxzp(:,:,:), cmyxp(:,:,:);
	complex,allocatable :: cmyzp(:,:,:), cmzxp(:,:,:), cmzyp(:,:,:);
	complex,allocatable :: cmxym(:,:,:), cmxzm(:,:,:), cmyxm(:,:,:);
	complex,allocatable :: cmyzm(:,:,:), cmzxm(:,:,:), cmzym(:,:,:);
	
	complex*16,allocatable :: eps_xy(:,:,:), eps_xz(:,:,:), eps_xi(:,:,:);
	complex*16,allocatable :: eps_yz(:,:,:), eps_yx(:,:,:), eps_yi(:,:,:);
	complex*16,allocatable :: eps_zx(:,:,:), eps_zy(:,:,:), eps_zi(:,:,:);

	complex*16,allocatable :: mu_xy(:,:,:), mu_xz(:,:,:), mu_xi(:,:,:);
	complex*16,allocatable :: mu_yz(:,:,:), mu_yx(:,:,:), mu_yi(:,:,:);
	complex*16,allocatable :: mu_zx(:,:,:), mu_zy(:,:,:), mu_zi(:,:,:);
    
    complex*16,allocatable :: se_xi(:,:,:),se_xy(:,:,:),se_xz(:,:,:);
    complex*16,allocatable :: se_yx(:,:,:),se_yi(:,:,:),se_yz(:,:,:);
    complex*16,allocatable :: se_zx(:,:,:),se_zy(:,:,:),se_zi(:,:,:);
    
    complex*16,allocatable :: sm_xi(:,:,:),sm_xy(:,:,:),sm_xz(:,:,:);
    complex*16,allocatable :: sm_yx(:,:,:),sm_yi(:,:,:),sm_yz(:,:,:);
    complex*16,allocatable :: sm_zx(:,:,:),sm_zy(:,:,:),sm_zi(:,:,:);   

end module global