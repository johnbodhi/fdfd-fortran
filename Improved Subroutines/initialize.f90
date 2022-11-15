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
    !dt = abs(cube_max_z(1) - cube_min_z(1));
    dt = dz*(tp-rp);
	Nx = floor(Sx/dx) + 1; Ny = floor(Sy/dy) + 1; Nz = floor(Sz/dz) + 1;
	Nx = (((Nx-1)/16)+1)*16;
	Ny = (((Ny-1)/16)+1)*16;
    !Nz = (((Nz-1)/16)+1)*16;
	
	Nxyz = Nx*Ny*Nz;
	Nu = 3*Nxyz;
	ell = 2;

    Nxm1=Nx-1; Nym1=Ny-1; Nzm1=Nz-1;
    Nxm2=Nx-2; Nym2=Ny-2; Nzm2=Nz-2;
    
	allocate(eps_r(Nxm1,Nym1,Nzm1),mu_r(Nxm1,Nym1,Nzm1),se(Nxm1,Nym1,Nzm1),sm(Nxm1,Nym1,Nzm1));
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
    
    allocate(se_xi(Nx,Ny,Nz),se_xy(Nx,Ny,Nz),se_xz(Nx,Ny,Nz));
    allocate(se_yx(Nx,Ny,Nz),se_yi(Nx,Ny,Nz),se_yz(Nx,Ny,Nz));
    allocate(se_zx(Nx,Ny,Nz),se_zy(Nx,Ny,Nz),se_zi(Nx,Ny,Nz));
    
    allocate(sm_xi(Nx,Ny,Nz),sm_xy(Nx,Ny,Nz),sm_xz(Nx,Ny,Nz));
    allocate(sm_yx(Nx,Ny,Nz),sm_yi(Nx,Ny,Nz),sm_yz(Nx,Ny,Nz));
    allocate(sm_zx(Nx,Ny,Nz),sm_zy(Nx,Ny,Nz),sm_zi(Nx,Ny,Nz));

	eps_xy = eps_o; eps_xz = eps_o; eps_xi = eps_o;
	eps_yz = eps_o; eps_yx = eps_o; eps_yi = eps_o;
	eps_zx = eps_o; eps_zy = eps_o; eps_zi = eps_o;

	mu_xy = mu_o; mu_xz = mu_o; mu_xi = mu_o;
	mu_yz = mu_o; mu_yx = mu_o; mu_yi = mu_o;
	mu_zx = mu_o; mu_zy = mu_o; mu_zi = mu_o;

end subroutine initialize