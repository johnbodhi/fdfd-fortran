subroutine setPMLboundaries
	use global
	implicit none
	integer mi,mj,mk,pml_type;
	real sigma,sigma_max,RO,dpml;
    
    real se_xi_pml,se_xy_pml,se_xz_pml,se_yx_pml,se_yi_pml,se_yz_pml,se_zx_pml,se_zy_pml,se_zi_pml;
    real sm_xi_pml,sm_xy_pml,sm_xz_pml,sm_yx_pml,sm_yi_pml,sm_yz_pml,sm_zx_pml,sm_zy_pml,sm_zi_pml;    
    
    real eps_xi_pml,eps_xy_pml,eps_xz_pml,eps_yx_pml,eps_yi_pml,eps_yz_pml,eps_zx_pml,eps_zy_pml,eps_zi_pml;
    real mu_xi_pml,mu_xy_pml,mu_xz_pml,mu_yx_pml,mu_yi_pml,mu_yz_pml,mu_zx_pml,mu_zy_pml,mu_zi_pml;    
    
    real k_xi,k_xy,k_xz,k_yx,k_yi,k_yz,k_zx,k_zy,k_zi;
    
    
! To set the boundaries PML
    pml_type = 2; ! pml_type = 0 constant cond., pml_type = 1 linear cond., pml_type = 2 parabolic cond.
    dpml = Npml*dx; RO = 1e-7; sigma_max = -eps_o*c*(pml_type+1)*log(RO)/(2*dpml);
    

! set the x=0 PML
do mi=0,Npml,1
			se_yx_pml = sigma_max * ((mi*dx)/dpml)**(pml_type+1);
			se_zx_pml = se_yx_pml;
			eps_yx (Npml-mi+1,:,:) = eps_o - j * se_yx_pml/ w; 
			eps_zx (Npml-mi+1,:,:) = eps_o - j * se_zx_pml/ w;     
			se_yx  (Npml-mi+1,:,:) = eps_o - j * se_yx_pml/ w;
			se_zx  (Npml-mi+1,:,:) = eps_o - j * se_zx_pml/ w;
end do
do mi=0,Npml-1,1
			sm_yx_pml = sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
			sm_zx_pml = sm_yx_pml;
			mu_yx  (Npml-mi,:,:) =  (eps_o - j * sm_yx_pml/ w) * mu_o / eps_o;
			mu_zx  (Npml-mi,:,:) =  (eps_o - j * sm_zx_pml/ w) * mu_o / eps_o;
			sm_yx  (Npml-mi,:,:) =  (eps_o - j * sm_yx_pml/ w) * mu_o / eps_o;
			sm_zx  (Npml-mi,:,:) =  (eps_o - j * sm_zx_pml/ w) * mu_o / eps_o;            
end do
do mi=0,Npml-1,1
			eps_yx_pml = eps_o - j * se_yx_pml/ w;
			eps_zx_pml = eps_o - j * se_zx_pml/ w;            
			mu_yx_pml  = eps_o - j * se_yx_pml/ w;
			mu_zx_pml  = eps_o - j * se_zx_pml/ w;            
			k_yx = real(w*sqrt(mu_yx_pml*eps_yx_pml/2)*sqrt(sqrt(1+(se_yx_pml/(eps_yx_pml*w))**2)-1));
			k_zx = real(w*sqrt(mu_zx_pml*eps_zx_pml/2)*sqrt(sqrt(1+(se_zx_pml/(eps_zx_pml*w))**2)-1));
			Einc_x (Npml-mi+1,:,:) =  exp(-j*k_yx*mi*dz); 
			Hinc_y (Npml-mi,:,:)   =  exp(-j*k_zx*(mi+0.5)*dz);   
end do

!  set the x=N boundary PML
do mi=0,Npml,1
			se_yx_pml = sigma_max * ((mi*dx)/dpml)**(pml_type+1);
			se_zx_pml = se_yx_pml;
			eps_yx (mi+Nx-Npml,:,:) = eps_o - j * se_yx_pml/ w;      
			eps_zx (mi+Nx-Npml,:,:) = eps_o - j * se_zx_pml/ w;  
			se_yx  (mi+Nx-Npml,:,:) = eps_o - j * se_yx_pml/ w;
			se_zx  (mi+Nx-Npml,:,:) = eps_o - j * se_zx_pml/ w;    
end do   
do mi=0,Npml-1,1
			sm_yx_pml = sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
			sm_zx_pml = sm_yx_pml;
			mu_yx  (mi+Nx-Npml,:,:) = (eps_o - j * sm_yx_pml/ w) * mu_o / eps_o;
			mu_zx  (mi+Nx-Npml,:,:) = (eps_o - j * sm_zx_pml/ w) * mu_o / eps_o;
			sm_yx  (mi+Nx-Npml,:,:) = (eps_o - j * sm_yx_pml/ w) * mu_o / eps_o;
			sm_zx  (mi+Nx-Npml,:,:) = (eps_o - j * sm_zx_pml/ w) * mu_o / eps_o;  
end do   
do mi=0,Npml-1,1
			eps_yx_pml = eps_o - j * se_yx_pml/ w;
			eps_zx_pml = eps_o - j * se_zx_pml/ w;            
			mu_yx_pml  = eps_o - j * se_yx_pml/ w;
			mu_zx_pml  = eps_o - j * se_zx_pml/ w;            
			k_yx = real(w*sqrt(mu_yx_pml*eps_yx_pml/2)*sqrt(sqrt(1+(se_yx_pml/(eps_yx_pml*w))**2)-1));
			k_zx = real(w*sqrt(mu_zx_pml*eps_zx_pml/2)*sqrt(sqrt(1+(se_zx_pml/(eps_zx_pml*w))**2)-1));
			Einc_x (mi+Nx-Npml,:,:) =  exp(-j*k_yx*mi*dz); 
			Hinc_y (mi+Nx-Npml,:,:) =  exp(-j*k_zx*(mi+0.5)*dz);   
end do

! set the y=0 PML
do mi=0,Npml,1
			se_xy_pml = sigma_max * ((mi*dx)/dpml)**(pml_type+1);
			se_zy_pml = se_xy_pml;
			eps_xy (:,Npml-mi+1,:) = eps_o - j * se_xy_pml/ w;     
			eps_zy (:,Npml-mi+1,:) = eps_o - j * se_zy_pml/ w;    
			se_xy  (:,Npml-mi+1,:) = eps_o - j * se_xy_pml/ w;
			se_zy  (:,Npml-mi+1,:) = eps_o - j * se_zy_pml/ w;   
end do
do mi=0,Npml-1,1
			sm_xy_pml = sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
			sm_zy_pml = sm_xy_pml;
			mu_xy  (:,Npml-mi,:) = (eps_o - j * sm_xy_pml/ w) * mu_o / eps_o;
			mu_zy  (:,Npml-mi,:) = (eps_o - j * sm_zy_pml/ w) * mu_o / eps_o;
			sm_xy  (:,Npml-mi,:) = eps_o - j * sm_xy_pml/ w;
			sm_zy  (:,Npml-mi,:) = eps_o - j * sm_zy_pml/ w;  
end do
do mi=0,Npml-1,1
			eps_xy_pml = eps_o - j * se_xy_pml/ w;
			eps_zy_pml = eps_o - j * se_zy_pml/ w;            
			mu_xy_pml  = eps_o - j * se_xy_pml/ w;
			mu_zy_pml  = eps_o - j * se_zy_pml/ w;            
			k_xy = real(w*sqrt(mu_xy_pml*eps_xy_pml/2)*sqrt(sqrt(1+(se_xy_pml/(eps_xy_pml*w))**2)-1));
			k_zy = real(w*sqrt(mu_zy_pml*eps_zy_pml/2)*sqrt(sqrt(1+(se_zy_pml/(eps_zy_pml*w))**2)-1));
			Einc_x (:,Npml-mi+1,:) =  exp(-j*k_xy*mi*dz); 
			Hinc_y (:,Npml-mi,:)   =  exp(-j*k_zy*(mi+0.5)*dz);   
end do

!  set the y=N boundary PML
do mi=0,Npml,1
			se_xy_pml = sigma_max * ((mi*dx)/dpml)**(pml_type+1);
			se_zy_pml = se_xy_pml;
			eps_xy (:,mi+Ny-Npml,:) = eps_o - j * se_xy_pml / w;     
			eps_zy (:,mi+Ny-Npml,:) = eps_o - j * se_zy_pml / w;   
			se_xy  (:,mi+Ny-Npml,:) = eps_o - j * se_xy_pml/ w;
			se_zy  (:,mi+Ny-Npml,:) = eps_o - j * se_zy_pml/ w;     
end do   
do mi=0,Npml-1,1
			sm_xy_pml = sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
			sm_xz_pml = sm_xy_pml;
			mu_xy  (:,mi+Ny-Npml,:) = (eps_o - j * sm_xy_pml/ w) * mu_o / eps_o;
			mu_zy  (:,mi+Ny-Npml,:) = (eps_o - j * sm_zy_pml/ w) * mu_o / eps_o;
			sm_xy  (:,mi+Ny-Npml,:) = eps_o - j * sm_xy_pml/ w;
			sm_zy  (:,mi+Ny-Npml,:) = eps_o - j * sm_zy_pml/ w;     
end do  
do mi=0,Npml-1,1
			eps_xy_pml = eps_o - j * se_xy_pml/ w;
			eps_zy_pml = eps_o - j * se_zy_pml/ w;            
			mu_xy_pml  = eps_o - j * se_xy_pml/ w;
			mu_zy_pml  = eps_o - j * se_zy_pml/ w;            
			k_xy = real(w*sqrt(mu_xy_pml*eps_xy_pml/2)*sqrt(sqrt(1+(se_xy_pml/(eps_xy_pml*w))**2)-1));
			k_zy = real(w*sqrt(mu_zy_pml*eps_zy_pml/2)*sqrt(sqrt(1+(se_zy_pml/(eps_zy_pml*w))**2)-1));
			Einc_x (:,mi+Ny-Npml,:) =  exp(-j*k_xy*mi*dz); 
			Hinc_y (:,mi+Ny-Npml,:) =  exp(-j*k_zy*(mi+0.5)*dz);   
end do

! set the z=0 PML
do mi=0,Npml,1
			se_xz_pml = sigma_max * ((mi*dx)/dpml)**(pml_type+1);
			se_yz_pml = se_xz_pml;
			eps_xz (:,:,Npml-mi+1) = eps_o - j * se_xz_pml/ w;     
			eps_yz (:,:,Npml-mi+1) = eps_o - j * se_yz_pml/ w;   
			se_xz  (:,:,Npml-mi+1) = eps_o - j * se_xz_pml/ w;
			se_yz  (:,:,Npml-mi+1) = eps_o - j * se_yz_pml/ w; 
end do
do mi=0,Npml-1,1
			sm_xz_pml = sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
			sm_yz_pml = sm_xz_pml;
			mu_xz  (:,:,Npml-mi) = (eps_o - j * sm_xz_pml/ w) * mu_o / eps_o;
			mu_yz  (:,:,Npml-mi) = (eps_o - j * sm_yz_pml/ w) * mu_o / eps_o;
			sm_xz  (:,:,Npml-mi) = eps_o - j * sm_xz_pml/ w;
			sm_yz  (:,:,Npml-mi) = eps_o - j * sm_yz_pml/ w;    
end do
do mi=0,Npml-1,1
			eps_xz_pml = eps_o - j * se_xz_pml/ w;
			eps_yz_pml = eps_o - j * se_yz_pml/ w;            
			mu_xz_pml  = eps_o - j * se_xz_pml/ w;
			mu_yz_pml  = eps_o - j * se_yz_pml/ w;            
			k_xz = real(w*sqrt(mu_xz_pml*eps_xz_pml/2)*sqrt(sqrt(1+(se_xz_pml/(eps_xz_pml*w))**2)-1));
			k_yz = real(w*sqrt(mu_yz_pml*eps_yz_pml/2)*sqrt(sqrt(1+(se_yz_pml/(eps_yz_pml*w))**2)-1));
			Einc_x (:,:,Npml-mi+1) =  exp(-j*k_xz*mi*dz); 
			Hinc_y (:,:,Npml-mi)   =  exp(-j*k_yz*(mi+0.5)*dz);   
end do

!  set the z=N boundary PML
do mi=0,Npml,1
			se_xz_pml = sigma_max * ((mi*dx)/dpml)**(pml_type+1);
			se_yz_pml = se_xz_pml;
			eps_xz (:,:,mi+Nz-Npml) = eps_o - j * se_xz_pml/ w;     
			eps_yz (:,:,mi+Nz-Npml) = eps_o - j * se_yz_pml/ w;   
			se_xz  (:,:,mi+Nz-Npml) = eps_o - j * se_xz_pml/ w;
			se_yz  (:,:,mi+Nz-Npml) = eps_o - j * se_yz_pml/ w;   
end do   
do mi=0,Npml-1,1
			sm_xz_pml = sigma_max * (((mi+0.5)*dx)/dpml)**(pml_type+1);
			sm_yz_pml = sm_xz_pml;
			mu_xz  (:,:,mi+Nz-Npml) = (eps_o - j * sm_xz_pml/ w) * mu_o / eps_o;
			mu_yz  (:,:,mi+Nz-Npml) = (eps_o - j * sm_yz_pml/ w) * mu_o / eps_o;
			sm_xz  (:,:,mi+Nz-Npml) = eps_o - j * sm_xz_pml/ w;
			sm_yz  (:,:,mi+Nz-Npml) = eps_o - j * sm_yz_pml/ w;      
end do 
do mi=0,Npml-1,1
			eps_xz_pml = eps_o - j * se_xz_pml/ w;
			eps_yz_pml = eps_o - j * se_yz_pml/ w;            
			mu_xz_pml  = eps_o - j * se_xz_pml/ w;
			mu_yz_pml  = eps_o - j * se_yz_pml/ w;            
			k_xz = real(w*sqrt(mu_xz_pml*eps_xz_pml/2)*sqrt(sqrt(1+(se_xz_pml/(eps_xz_pml*w))**2)-1));
			k_yz = real(w*sqrt(mu_yz_pml*eps_yz_pml/2)*sqrt(sqrt(1+(se_yz_pml/(eps_yz_pml*w))**2)-1));
			Einc_x (:,:,mi+Nz-Npml)   =  exp(-j*k_xz*mi*dz); 
			Hinc_y (:,:,mi+Nz-Npml)   =  exp(-j*k_yz*(mi+0.5)*dz);   
end do

end subroutine setPMLboundaries