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
    
    allocate(pex(Nx,Ny,Nz),pey(Nx,Ny,Nz),pez(Nx,Ny,Nz),phx(Nx,Ny,Nz),phy(Nx,Ny,Nz),phz(Nx,Ny,Nz));
	
	tmpx = (0.0,0.0); tmpy = (0.0,0.0); tmpz = (0.0,0.0);
	
    Cexhy = 1/(j*w*dz*eps_xz); Cexhz = 1/(j*w*dy*eps_xy); Cexex = -(eps_xi - eps_o) / eps_xi;  
    Ceyhz = 1/(j*w*dx*eps_yx); Ceyhx = 1/(j*w*dz*eps_yz); Ceyey = -(eps_yi - eps_o) / eps_yi;  
    Cezhx = 1/(j*w*dy*eps_zy); Cezhy = 1/(j*w*dx*eps_zx); Cezez = -(eps_zi - eps_o) / eps_zi;  
    
    Chxey = 1/(j*w*dz*mu_xz); Chxez = 1/(j*w*dy*mu_xy); Chxhx = (mu_xi - mu_o) / mu_xi;  
    Chyez = 1/(j*w*dx*mu_yx); Chyex = 1/(j*w*dz*mu_yz); Chyhy = (mu_yi - mu_o) / mu_yi;  
    Chzex = 1/(j*w*dy*mu_zy); Chzey = 1/(j*w*dx*mu_zx); Chzhz = (mu_zi - mu_o) / mu_zi;  

    ! divergence enforcing coefficients...        
    
    Cexhx = mu_xi/(j*w*dx*s*eps_xi); Cexhy = 1/(j*w*dz*s*eps_xz); Cexhz = 1/(j*w*dy*s*eps_xy); Cexex = -(eps_xi - eps_o) / eps_xi;  
    Ceyhy = mu_yi/(j*w*dy*s*eps_yi); Ceyhz = 1/(j*w*dx*s*eps_yx); Ceyhx = 1/(j*w*dz*s*eps_yz); Ceyey = -(eps_yi - eps_o) / eps_yi;  
    Cezhz = mu_zi/(j*w*dz*s*eps_zi); Cezhx = 1/(j*w*dy*s*eps_zy); Cezhy = 1/(j*w*dx*s*eps_zx); Cezez = -(eps_zi - eps_o) / eps_zi;  
    
    Chxex = eps_xi/(j*w*dx*s*mu_xi); Chxey = 1/(j*w*dz*s*mu_xz); Chxez = 1/(j*w*dy*s*mu_xy); Chxhx = (mu_xi - mu_o) / mu_xi;  
    Chyey = eps_yi/(j*w*dy*s*mu_yi); Chyez = 1/(j*w*dx*s*mu_yx); Chyex = 1/(j*w*dz*s*mu_yz); Chyhy = (mu_yi - mu_o) / mu_yi;  
    Chzez = eps_zi/(j*w*dz*s*mu_zi); Chzex = 1/(j*w*dy*s*mu_zy); Chzey = 1/(j*w*dx*s*mu_zx); Chzhz = (mu_zi - mu_o) / mu_zi;  
    
    
    ! totally broken symmetry divergence enforcing coefficients...
    
    Cexhy = 1/(j*w*dz*eps_xz); Cexhz = 1/(j*w*dy*eps_xy); 
    Ceyhz = 1/(j*w*dx*eps_yx); Ceyhx = 1/(j*w*dz*eps_yz); 
    Cezhx = 1/(j*w*dy*eps_zy); Cezhy = 1/(j*w*dx*eps_zx); 
    
    Cexex = ((eps_o - eps_xi) / eps_xi + &
        1 / (j*w*dx*eps_xi)*(eps_xi(2:nx,2:ny,2:nz)-eps_xi(1:nx-1,1:ny,1:nz))) / &
        (j*w*eps_xi)    
    Ceyey = ((eps_o - eps_yi) / eps_yi + &
        1 / (j*w*dy*eps_yi)*(eps_yi(2:nx,2:ny,2:nz)-eps_yi(1:nx,1:ny-1,1:nz))) / &
        (j*w*eps_yi)    
    Cezez = ((eps_o - eps_zi) / eps_zi + &
        1 / (j*w*dz*eps_zi)*(eps_zi(2:nx,2:ny,2:nz)-eps_zi(1:nx,1:ny,1:nz-1))) / &
        (j*w*eps_zi)
    
    Chxey = 1/(j*w*dz*mu_xz); Chxez = 1/(j*w*dy*mu_xy); 
    Chyez = 1/(j*w*dx*mu_yx); Chyex = 1/(j*w*dz*mu_yz); 
    Chzex = 1/(j*w*dy*mu_zy); Chzey = 1/(j*w*dx*mu_zx); 
    
    Chxhx = ((mu_xi - mu_o) / mu_xi + &
        1 / (j*w*dx*mu_xi)*(mu_xi(2:nx,2:ny,2:nz)-mu_xi(1:nx-1,1:ny,1:nz))) / &
        (j*w*mu_xi);
    Chyhy = ((mu_yi - mu_o) / mu_yi + &
        1 / (j*w*dy*mu_yi)*(mu_yi(2:nx,2:ny,2:nz)-mu_yi(1:nx,1:ny-1,1:nz))) / &
        (j*w*mu_yi);
    Chzhz = ((mu_zi - mu_o) / mu_zi + &
        1 / (j*w*dz*mu_zi)*(mu_zi(2:nx,2:ny,2:nz)-mu_zi(1:nx,1:ny,1:nz-1))) / &
        (j*w*mu_zi);
    
    ! partially broken symmetry divergence enforcing coefficients...
    
    Cexhy = 1/(j*w*dz*eps_xz); Cexhz = 1/(j*w*dy*eps_xy); 
    Ceyhz = 1/(j*w*dx*eps_yx); Ceyhx = 1/(j*w*dz*eps_yz); 
    Cezhx = 1/(j*w*dy*eps_zy); Cezhy = 1/(j*w*dx*eps_zx); 
    
    Cexex = ( (eps_o - eps_xi) / eps_xi + &
        1 / (j*w*dx*eps_xi)*( eps_xi(2:nx,2:ny,2:nz)-eps_xi(1:nx-1,1:ny,1:nz)) ) / &
        ( 1 + 1/(j*w*dx*eps_xi)*(eps_xi(2:nx,2:ny,2:nz)-eps_xi(1:nx-1,1:ny,1:nz)));
    Ceyey = ( (eps_o - eps_yi) / eps_yi + &
         1 / (j*w*dy*eps_yi)*(eps_yi(2:nx,2:ny,2:nz)-eps_yi(1:nx,1:ny-1,1:nz)) ) / &
         ( 1 + 1/(j*w*dy*eps_yi)*(eps_yi(2:nx,2:ny,2:nz)-eps_yi(1:nx,1:ny-1,1:nz)));
    Cezez = ( (eps_o - eps_zi) / eps_zi + &
         1 / (j*w*dz*eps_zi)*(eps_zi(2:nx,2:ny,2:nz)-eps_zi(1:nx,1:ny,1:nz-1)) ) / &
         ( 1 + 1/(j*w*dz*eps_zi)*(eps_zi(2:nx,2:ny,2:nz)-eps_zi(1:nx,1:ny,1:nz-1)));
    
    Chxey = 1/(j*w*dz*mu_xz); Chxez = 1/(j*w*dy*mu_xy); 
    Chyez = 1/(j*w*dx*mu_yx); Chyex = 1/(j*w*dz*mu_yz); 
    Chzex = 1/(j*w*dy*mu_zy); Chzey = 1/(j*w*dx*mu_zx); 
    
    Chxhx = ( (mu_xi - mu_o) / mu_xi + &
         1 / (j*w*dx*mu_xi)*(mu_xi(2:nx,2:ny,2:nz)-mu_xi(1:nx-1,1:ny,1:nz)) ) / &
         ( 1 - 1/(j*w*dx*mu_xi)*(mu_xi(2:nx,2:ny,2:nz)-mu_xi(1:nx-1,1:ny,1:nz)));
    Chyhy = ( (mu_yi - mu_o) / mu_yi + &
         1 / (j*w*dy*mu_yi)*(mu_yi(2:nx,2:ny,2:nz)-mu_yi(1:nx,1:ny-1,1:nz)) ) / &
         ( 1 - 1/(j*w*dy*mu_yi)*(mu_yi(2:nx,2:ny,2:nz)-mu_yi(1:nx,1:ny-1,1:nz)));
    Chzhz = ( (mu_zi - mu_o) / mu_zi + &
         1 / (j*w*dz*mu_zi)*(mu_zi(2:nx,2:ny,2:nz)-mu_zi(1:nx,1:ny,1:nz-1)) ) / &
         ( 1 - 1/(j*w*dz*mu_zi)*(mu_zi(2:nx,2:ny,2:nz)-mu_zi(1:nx,1:ny,1:nz-1)));
    
    
    ! totally symmetric divergence enforcing coefficients...
    
    Cexhy = 1/(( j * w * eps_xz + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx ) * dz );     
    Cexhz = 1/(( j * w * eps_xy + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx ) * dy );    
    Cexex = ( j * w * ( eps_o - eps_xi ) + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx ) / &
        & ( j * w * eps_xi + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx );          
    
    Ceyhz = 1/(( j * w * eps_yx + ( eps_yi(2:nx,2:ny,2:nz) - eps_yi(1:nx,1:ny-1,1:nz) ) / dy ) * dx );    
    Ceyhx = 1/(( j * w * eps_yz + ( eps_yi(2:nx,2:ny,2:nz) - eps_yi(1:nx,1:ny-1,1:nz) ) / dy ) * dz );   
    Ceyey = ( j * w * ( eps_o - eps_xi ) + ( eps_yi(2:nx,2:ny,2:nz) - eps_yi(1:nx,1:ny-1,1:nz) ) / dy ) / &
        & ( j * w * eps_yi + ( eps_yi(2:nx,2:ny,2:nz) - eps_xi(1:nx,1:ny-1,1:nz) ) / dy );  
    
    Cezhx = 1/(( j * w * eps_zy + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz ) * dy );     
    Cezhy = 1/(( j * w * eps_zx + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz ) * dx );     
    Cezez = ( j * w * ( eps_o - eps_zi ) + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz ) / &
        & ( j * w * eps_zi + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz );  
    
    Chxey = 1/(( j * w * mu_xz - ( mu_xi(2:nx,2:ny,2:nz) - mu_xi(1:nx-1,1:ny,1:nz) ) / dx ) * dz );     
    Chxez = 1/(( j * w * mu_xy - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx-1,1:ny,1:nz) ) / dx ) * dy );    
    Chxhx = ( j * w * ( mu_o - mu_xi ) + ( mu_xi(2:nx,2:ny,2:nz) - mu_xi(1:nx-1,1:ny,1:nz) ) / dx ) / &
        & ( j * w * mu_xi + ( mu_zi(2:nx,2:ny,2:nz) - mu_xi(1:nx-1,1:ny,1:nz) ) / dx );  
    
    Chyez = 1/(( j * w * mu_yx - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny-1,1:nz) ) / dy ) * dx );    
    Chyex = 1/(( j * w * mu_xz - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny-1,1:nz) ) / dy ) * dz );    
    Chyhy = ( j * w * ( mu_o - mu_yi ) + ( mu_yi(2:nx,2:ny,2:nz) - mu_yi(1:nx,1:ny-1,1:nz-1) ) / dy ) / &
        & ( j * w * mu_yi + ( mu_yi(2:nx,2:ny,2:nz) - mu_xi(1:nx,1:ny-1,1:nz) ) / dy );  
    
    Chzex = 1/(( j * w * mu_zy - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz ) * dy );     
    Chzey = 1/(( j * w * mu_zx - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz ) * dx );   
    Chzhz = ( j * w * ( mu_o - mu_zi ) + ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz ) / &
        & ( j * w * mu_zi + ( mu_zi(1:nx,1:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz ); 
    
    
    ! totally symmetric divergene enforcing coefficients with conductivity...    
    
    Cexhy = 1/(( j * w * eps_xz + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx + se_xz) * dz );     
    Cexhz = 1/(( j * w * eps_xy + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx + se_xy) * dy );     
    Cexex = ( j * w * ( eps_o - eps_xi ) + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz)) / dx + se_xi ) / &
        & ( j * w * eps_xi + ( eps_xi(2:nx,2:ny,2:nz) - eps_xi(1:nx-1,1:ny,1:nz) ) / dx + se_xi );   
    
    Ceyhz = 1/(( j * w * eps_yx + ( eps_yi(2:nx,2:ny,2:nz) - eps_yi(1:nx,1:ny-1,1:nz) ) / dy + se_yx ) * dx );    
    Ceyhx = 1/(( j * w * eps_yz + ( eps_yi(2:nx,2:ny,2:nz) - eps_yi(1:nx,1:ny-1,1:nz) ) / dy + se_yz ) * dz );     
    Ceyey = ( j * w * ( eps_o - eps_xi ) + ( eps_yi(2:nx,2:ny,2:nz) - eps_yi(1:nx,1:ny-1,1:nz) ) / dy + se_yi ) / &
        & ( j * w * eps_yi + ( eps_yi(2:nx,2:ny,2:nz) - eps_xi(1:nx,1:ny-1,1:nz) ) / dy + se_yi );  
    
    Cezhx = 1/(( j * w * eps_zy + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz + se_zy ) * dy );     
    Cezhy = 1/(( j * w * eps_zx + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz + se_zx ) * dx );     
    Cezez = ( j * w * ( eps_o - eps_zi ) + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny,1:nz-1) ) / dz + se_zi ) / &
        & ( j * w * eps_zi + ( eps_zi(2:nx,2:ny,2:nz) - eps_zi(1:nx,1:ny-1,1:nz) ) / dz + se_zi );  
    
    Chxey = 1/(( j * w * mu_xz - ( mu_xi(2:nx,2:ny,2:nz) - mu_xi(1:nx-1,1:ny,1:nz) ) / dx + sm_xz ) * dz );     
    Chxez = 1/(( j * w * mu_xy - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx-1,1:ny,1:nz) ) / dx + sm_xy ) * dy );      
    Chxhx = ( j * w * ( mu_xi - mu_o ) + ( mu_xi(2:nx,2:ny,2:nz) - mu_xi(1:nx-1,1:ny,1:nz-1) ) / dx - sm_xi ) / &
        & ( j * w * mu_xi - ( mu_zi(2:nx,2:ny,2:nz) - mu_xi(1:nx-1,1:ny,1:nz) ) / dx + sm_xi );  
    
    Chyez = 1/(( j * w * mu_yx - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny-1,1:nz) ) / dy + sm_yx ) * dx );    
    Chyex = 1/(( j * w * mu_xz - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny-1,1:nz) ) / dy + sm_yz ) * dz );       
    Chyhy = ( j * w * ( mu_yi - mu_o ) + ( mu_yi(2:nx,2:ny,2:nz) - mu_yi(1:nx,1:ny-1,1:nz-1) ) / dy - sm_yi ) / &
        & ( j * w * mu_yi - ( mu_yi(2:nx,2:ny,2:nz) - mu_xi(1:nx,1:ny-1,1:nz) ) / dy + sm_yi );  
    
    Chzex = 1/(( j * w * mu_zy - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz + sm_zy ) * dy );     
    Chzey = 1/(( j * w * mu_zx - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz + sm_zy ) * dx  );      
    Chzhz = ( j * w * ( mu_zi - mu_o ) + ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz - sm_zi ) / &
        & ( j * w * mu_zi - ( mu_zi(2:nx,2:ny,2:nz) - mu_zi(1:nx,1:ny,1:nz-1) ) / dz + sm_zi );     
    
    
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
    
      
    ! field divergence equations V1...
    
    Bex(2:nx,2:ny,2:nz) = Bex(2:nx,2:ny,2:nz) - &
                    (-Cexhi(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)+Cexhi(2:nx,2:ny,2:nz)*Bhx(1:nx-1,2:ny,2:nz) &
                     -Cexhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)+Cexhz(2:nx,2:ny,2:nz)*Bhz(2:nx,1:ny-1,2:nz) &
                     +Cexhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)-Cexhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,1:nz-1) &
                     -Cexhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)+Cexhx(2:nx,2:ny,2:nz)*Bhx(1:nx-1,2:ny,2:nz));  
    
    Bey(2:nx,2:ny,2:nz) = Bey(2:nx,2:ny,2:nz) - &
                    (-Ceyhi(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)+Ceyhi(2:nx,2:ny,2:nz)*Bhy(2:nx,1:ny-1,2:nz) &
                     -Ceyhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)+Ceyhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,1:nz-1) &
                     +Ceyhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)-Ceyhz(2:nx,2:ny,2:nz)*Bhz(1:nx-1,2:ny,2:nz) &
                     -Ceyhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)+Ceyhy(2:nx,2:ny,2:nz)*Bhy(2:nx,1:ny-1,2:nz));  
                     
    Bez(2:nx,2:ny,2:nz) = Bez(2:nx,2:ny,2:nz) - &
                    (-Cezhi(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)+Cezhi(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,1:nz-1) &
					 -Cezhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)+Cezhy(2:nx,2:ny,2:nz)*Bhy(1:nx-1,2:ny,2:nz) &
                     +Cezhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)-Cezhx(2:nx,2:ny,2:nz)*Bhx(2:nx,1:ny-1,2:nz) &
                     -Cezhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)+Cezhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,1:nz-1));  

    
     ! field divergence equations V2...
    
    Bex(2:nx,2:ny,2:nz) = Bex(2:nx,2:ny,2:nz) - &
                    (-Cexhi(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)-Cexhi(2:nx,2:ny,2:nz)*Bhx(1:nx-3,2:ny,2:nz) &
                     -Cexhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)+Cexhz(2:nx,2:ny,2:nz)*Bhz(2:nx,1:ny-1,2:nz) &
                     +Cexhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)-Cexhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,1:nz-1) &
                     -Cexhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)+Cexhx(2:nx,2:ny,2:nz)*Bhx(1:nx-3,2:ny,2:nz));  
    
    Bey(2:nx,2:ny,2:nz) = Bey(2:nx,2:ny,2:nz) - &
                    (-Ceyhi(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)-Ceyhi(2:nx,2:ny,2:nz)*Bhy(2:nx,1:ny-3,2:nz) &
                     -Ceyhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)+Ceyhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,1:nz-1) &
                     +Ceyhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)-Ceyhz(2:nx,2:ny,2:nz)*Bhz(1:nx-1,2:ny,2:nz) &
                     -Ceyhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)+Ceyhy(2:nx,2:ny,2:nz)*Bhy(2:nx,1:ny-3,2:nz));  
                     
    Bez(2:nx,2:ny,2:nz) = Bez(2:nx,2:ny,2:nz) - &
                    (-Cezhi(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)-Cezhi(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,1:nz-3) &
					 -Cezhy(2:nx,2:ny,2:nz)*Bhy(2:nx,2:ny,2:nz)+Cezhy(2:nx,2:ny,2:nz)*Bhy(1:nx-1,2:ny,2:nz) &
                     +Cezhx(2:nx,2:ny,2:nz)*Bhx(2:nx,2:ny,2:nz)-Cezhx(2:nx,2:ny,2:nz)*Bhx(2:nx,1:ny-1,2:nz) &
                     -Cezhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,2:nz)+Cezhz(2:nx,2:ny,2:nz)*Bhz(2:nx,2:ny,1:nz-3));  
    
    
    
  ! periodic boundary conditions in lieu of perfectly matched layers on transverse boundaries 
    
 !   Bex(1:nx,1,2:nz) = Bex(1:nx,1,2:nz) - &
 !                   (-Cexhz(1:nx,1,2:nz)*Bhz(1:nx,1,2:nz)+Cexhz(1:nx,1,2:nz)*Bhz(1:nx,nyp1,2:nz) &
 !                    +Cexhy(1:nx,1,2:nz)*Bhy(1:nx,1,2:nz)-Cexhy(1:nx,1,2:nz)*Bhy(1:nx,1,1:nz-1)*exp(j*ky*Py));  
 !   
 !   Bex(1:nx,nyp1,2:nz) = Bex(1:nx,1,2:nz)*exp(-j*ky*Py);
 !   
 !   Bey(1,1:ny,2:nz) = Bey(1,1:ny,2:nz) - &
 !                   (-Ceyhx(1,1:ny,2:nz)*Bhx(1,1:ny,2:nz)+Ceyhx(1,1:ny,2:nz)*Bhx(1,1:ny,1:nz-1) &
 !                    +Ceyhz(1,1:ny,2:nz)*Bhz(1,1:ny,2:nz)-Ceyhz(1,1:ny,2:nz)*Bhz(nxp1,1:ny,2:nz)*exp(j*kx*Px));  
 !   
 !   Bey(nxp1,1:ny,2:nz) = Bey(1,1:ny,2:nz)*exp(-j*kx*Px);
 !                   
 !   Bez(1,2:ny,1:nz) = Bez(1,2:ny,1:nz) - &
 !                   (-Cezhy(1,2:ny,1:nz)*Bhy(1,2:ny,1:nz)+Cezhy(1,2:ny,1:nz)*Bhy(nxp1,2:ny,1:nz)*exp(j*kx*Px) &
 !                    +Cezhx(1,2:ny,1:nz)*Bhx(1,2:ny,1:nz)-Cezhx(1,2:ny,1:nz)*Bhx(1,1:ny-1,1:nz));       
 !   
 !   Bez(nxp1,2:ny,1:nz) = Bez(1,2:ny,1:nz)*exp(-j*kx*Px);
 !   
 !   Bez(2:nx,1,1:nz) = Bez(2:nx,1,1:nz) - &
 !                   (-Cezhy(2:nx,1,1:nz)*Bhy(2:nx,1,1:nz)+Cezhy(2:nx,1,1:nz)*Bhy(1:nx-1,1,1:nz) &
 !                    +Cezhx(2:nx,1,1:nz)*Bhx(2:nx,1,1:nz)-Cezhx(2:nx,1,1:nz)*Bhx(2:nx,nyp1,1:nz)*exp(j*ky*Py));    
 !   
 !   Bez(1:nx,nyp1,1:nz) = Bez(2:nx,1,1:nz)*exp(-j*ky*Py);    
 !    
 !   Bez(1,1,1:nz) = Bez(1,1,1:nz) - &
 !                   (-Cezhy(1,1,1:nz)*Bhy(1,1,1:nz)+Cezhy(1,1,1:nz)*Bhy(nxp1,1,1:nz)*exp(j*kx*Px) &
 !                    +Cezhx(1,1,1:nz)*Bhx(1,1,1:nz)-Cezhx(1,1,1:nz)*Bhx(1,nyp1,1:nz)*exp(j*ky*Py));   
 !     
	!Bez(nxp1,1,1:nz) = Bez(1,1,1:nz)*exp(-j*kx*Px);
 !   
 !   Bez(1,nyp1,1:nz) = Bez(1,1,1:nz)*exp(-j*ky*Py);
 !   
 !   Bez(nxp1,nyp1,1:nz) = Bez(1,1,1:nz)*exp(-j*kx*Px)*exp(-j*ky*Py);
    
    
   ! floquet phase corrections for flat surfaces
    
   ! do mi = 1,nx,1
   !     do mj = 1,ny,1                    
			!pex(mi,mj,1) = exp( j * ko * ( mi + 0.5 ) * dx ) * exp( j * ko * ( mj - 1 ) * dy );
			!pey(mi,mj,1) = exp( j * ko * ( mi - 1 ) * dx ) * exp( j * ko * ( mj + 0.5 ) * dy );
			!pez(mi,mj,1) = exp( j * ko * ( mi - 1 ) * dx ) * exp( j * ko * ( mj - 1 ) * dy );
			!phz(mi,mj,1) = exp( j * ko * ( mi + 0.5 ) * dx ) * exp( j * ko * ( mj + 0.5 ) * dy );     
   !     enddo   
   ! enddo    
   ! 
   ! phx = pey;
   ! phy = pex;    
   
    ! reflection and transmission plane fields for tem waves...
    
    !do mi = 1,nx,1
    !    do mj = 1,ny,1             
    !        Bex(mi,mj,ref) = Bex(mi,mj,ref) * pex(mi,mj,1);
    !        Bey(mi,mj,ref) = Bey(mi,mj,ref) * pey(mi,mj,1);
    !        Bez(mi,mj,ref) = 0.5 * ( Bez(mi,mj,ref) + Bez(mi,mj,ref-1) ) * pez(mi,mj,1);
    !            
    !        Bhx(mi,mj,ref) = 0.5 * ( Bhx(mi,mj,ref) + Bhx(mi,mj,ref-1) ) * phx(mi,mj,1);
    !        Bhy(mi,mj,ref) = 0.5 * ( Bhy(mi,mj,ref) + Bhy(mi,mj,ref-1) ) * phy(mi,mj,1);
    !        Bhz(mi,mj,ref) = Bhz(mi,mj,ref) * phz(mi,mj,1);       
    !        
    !        
    !        Bex(mi,mj,tra) = Bex(mi,mj,tra) * pex(mi,mj,1);
    !        Bey(mi,mj,tra) = Bey(mi,mj,tra) * pey(mi,mj,1);
    !        Bez(mi,mj,tra) = 0.5 * ( Bez(mi,mj,tra) + Bez(mi,mj,tra-1) ) * pez(mi,mj,1);
    !            
    !        Bhx(mi,mj,tra) = 0.5 * ( Bhx(mi,mj,tra) + Bhx(mi,mj,tra-1) ) * phx(mi,mj,1);
    !        Bhy(mi,mj,tra) = 0.5 * ( Bhy(mi,mj,tra) + Bhy(mi,mj,tra-1) ) * phy(mi,mj,1);
    !        Bhz(mi,mj,tra) = Bhz(mi,mj,tra) * phz(mi,mj,1);    
    !        
    !        Bex(mi,mj,rp) = Bex(mi,mj,rp) * pex(mi,mj,1);
   !         Bey(mi,mj,rp) = Bey(mi,mj,rp) * pey(mi,mj,1);
   !         Bez(mi,mj,rp) = 0.5 * ( Bez(mi,mj,rp) + Bez(mi,mj,rp-1) ) * pez(mi,mj,1);
   !             
   !         Bhx(mi,mj,rp) = 0.5 * ( Bhx(mi,mj,rp) + Bhx(mi,mj,rp-1) ) * phx(mi,mj,1);
   !         Bhy(mi,mj,rp) = 0.5 * ( Bhy(mi,mj,rp) + Bhy(mi,mj,rp-1) ) * phy(mi,mj,1);
   !         Bhz(mi,mj,rp) = Bhz(mi,mj,rp) * phz(mi,mj,1);       
   !         
   !         
   !         Bex(mi,mj,tp) = Bex(mi,mj,tp) * pex(mi,mj,1);
   !         Bey(mi,mj,tp) = Bey(mi,mj,tp) * pey(mi,mj,1);
   !         Bez(mi,mj,tp) = 0.5 * ( Bez(mi,mj,tp) + Bez(mi,mj,tp-1) ) * pez(mi,mj,1);
   !             
   !         Bhx(mi,mj,tp) = 0.5 * ( Bhx(mi,mj,tp) + Bhx(mi,mj,tp-1) ) * phx(mi,mj,1);
   !         Bhy(mi,mj,tp) = 0.5 * ( Bhy(mi,mj,tp) + Bhy(mi,mj,tp-1) ) * phy(mi,mj,1);
   !         Bhz(mi,mj,tp) = Bhz(mi,mj,tp) * phz(mi,mj,1);    
   !         
   !         !Einc_x(mi,mj,sp) = Einc_x(mi,mj,sp) * pex(mi,mj,1);
   !         !Einc_y(mi,mj,sp) = Einc_y(mi,mj,sp) * pey(mi,mj,1);
   !         !Einc_z(mi,mj,sp) = 0.5 * ( Einc_z(mi,mj,sp) + Einc_z(mi,mj,sp-1) ) * pez(mi,mj,rp);
   !         
   !         !Hinc_x(mi,mj,sp) = Hinc_x(mi,mj,sp) * phx(mi,mj,1);
   !         !Hinc_y(mi,mj,sp) = Hinc_y(mi,mj,sp) * phy(mi,mj,1);
   !         !Hinc_z(mi,mj,sp) = 0.5 * ( Hinc_x(mi,mj,sp) + Hinc_z(mi,mj,sp-1) ) * phz(mi,mj,sp);
    !    enddo
    !enddo
    
    
    ! floquet phase corrections for curved surfaces
    !
    ! do mi = 1,nx,1
    !    do mj = 1,ny,1          
    !        do mk = 1,nz,1
				!pex(mi,mj,mk) = exp( j * ko * ( mi + 0.5 ) * dx ) * exp( j * ko * ( mj - 1 ) * dy );
				!pey(mi,mj,mk) = exp( j * ko * ( mi - 1 ) * dx ) * exp( j * ko * ( mj + 0.5 ) * dy );
				!pez(mi,mj,mk) = exp( j * ko * ( mi - 1 ) * dx ) * exp( j * ko * ( mj - 1 ) * dy );
				!phz(mi,mj,mk) = exp( j * ko * ( mi + 0.5 ) * dx ) * exp( j * ko * ( mj + 0.5 ) * dy );     
    !        end do
    !    enddo   
    !enddo    
    !
    !phx = pey;
    !phy = pex;    
    
    
    ! reflection and transmission plane fields for tem waves...
    
    !do mi = 1,nx,1
    !    do mj = 1,ny,1  
    !        do mk = 1,nz,1
				!Bex(mi,mj,mk) = Bex(mi,mj,mk) * pex(mi,mj,mk);
				!Bey(mi,mj,mk) = Bey(mi,mj,mk) * pey(mi,mj,mk);
				!Bez(mi,mj,mk) = 0.5 * ( Bez(mi,mj,mk) + Bez(mi,mj,mk) ) * pez(mi,mj,mk);
    !            
				!Bhx(mi,mj,mk) = 0.5 * ( Bhx(mi,mj,mk) + Bhx(mi,mj,mk) ) * phx(mi,mj,mk);
				!Bhy(mi,mj,mk) = 0.5 * ( Bhy(mi,mj,mk) + Bhy(mi,mj,mk) ) * phy(mi,mj,mk);
				!Bhz(mi,mj,mk) = Bhz(mi,mj,mk) * phz(mi,mj,mk);               
    !        
				!Bex(mi,mj,mk) = Bex(mi,mj,mk) * pex(mi,mj,mk);
				!Bey(mi,mj,mk) = Bey(mi,mj,mk) * pey(mi,mj,mk);
				!Bez(mi,mj,mk) = 0.5 * ( Bez(mi,mj,mk) + Bez(mi,mj,mk) ) * pez(mi,mj,mk);
    !            
				!Bhx(mi,mj,mk) = 0.5 * ( Bhx(mi,mj,mk) + Bhx(mi,mj,mk) ) * phx(mi,mj,mk);
				!Bhy(mi,mj,mk) = 0.5 * ( Bhy(mi,mj,mk) + Bhy(mi,mj,mk) ) * phy(mi,mj,mk);
				!Bhz(mi,mj,mk) = Bhz(mi,mj,mk) * phz(mi,mj,mk); 
    !            Bex(mi,mj,rp) = Bex(mi,mj,rp) * pex(mi,mj,1);
    !  
    !        end do
    !    enddo
    !enddo
    
    Bvec(0*Nxyz+1:1*Nxyz) = reshape(Bex, (/Nxyz/));
    Bvec(1*Nxyz+1:2*Nxyz) = reshape(Bey, (/Nxyz/));
    Bvec(2*Nxyz+1:3*Nxyz) = reshape(Bez, (/Nxyz/));

end subroutine setCoefficients