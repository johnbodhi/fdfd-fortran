subroutine matvec(n,x,y)
    use global
    implicit none
    integer n; complex*16 x(nx,ny,nz,3), y(nx,ny,nz,3);

	tmpx(5:Nx,1:Nym1,1:Nzm1) = Chxez(5:Nxm1,1:Nym1,1:Nzm1)*(x(5:Nxm1,2:Ny,1:Nzm1,3)-x(5:Nxm1,1:Nym1,1:Nzm1,3)) &
				+ Chxey(5:Nxm1,1:Nym1,1:Nzm1)*(-x(5:Nxm1,1:Nym1,2:Nz,2)+x(5:Nxm1,1:Nym1,1:Nzm1,2)) &
				+ Chxex(5:Nx,1:Nym1,1:Nzm1)*(x(5:Nx,1:Nym1,1:Nzm1,1)+x(1:Nxm1-3,1:Nym1,1:Nzm1,1)) &
				- Chxei(5:Nx,1:Nym1,1:Nzm1)*(-x(5:Nx,1:Nym1,1:Nzm1,1)-x(1:Nxm1-3,1:Nym1,1:Nzm1,1));
 
	tmpy(1:Nxm1,5:Ny,1:Nzm1) = Chyex(1:Nxm1,5:Ny,1:Nzm1)*(x(1:Nxm1,5:Ny,2:Nz,1)-x(1:Nxm1,5:Ny,1:Nzm1,1)) &
				+ Chyez(1:Nxm1,5:Ny,1:Nzm1)*(-x(2:Nx,5:Ny,1:Nzm1,3)+x(1:Nxm1,5:Ny,1:Nzm1,3)) & 
				+ Chyey(1:Nxm1,5:Ny,1:Nzm1)*(x(1:Nxm1,5:Ny,1:Nzm1,2)+x(1:Nxm1,1:Nym1-3,1:Nzm1,2)) &
				- Chyei(1:Nxm1,5:Ny,1:Nzm1)*(-x(1:Nxm1,5:Ny,1:Nzm1,2)-x(1:Nxm1,1:Nym1-3,1:Nzm1,2));
 
	tmpz(1:Nxm1,1:Nym1,5:Nz) = Chzey(1:Nxm1,1:Nym1,5:Nz)*(x(2:Nx,1:Nym1,5:Nz,2)-x(1:Nxm1,1:Nym1,5:Nz,2)) &
				+ Chzex(1:Nxm1,1:Nym1,5:Nz)*(-x(1:Nxm1,2:Ny,5:Nz,1)+x(1:Nxm1,1:Nym1,5:Nz,1)) &
				+ Chzez(1:Nxm1,1:Nym1,5:Nz)*(x(1:Nxm1,1:Nym1,5:Nz,3)+x(1:Nxm1,1:Nym1,1:Nzm1-3,3)) &
				- Chzei(1:Nxm1,1:Nym1,5:Nz)*(-x(1:Nxm1,1:Nym1,5:Nz,3)-x(1:Nxm1,1:Nym1,1:Nzm1-3,3));
 
	y(5:Nx,2:Ny,2:Nz,1) = x(5:Nx,2:Ny,2:Nz,1) &
				- (Cexhz(5:Nx,2:Ny,2:Nz)*(-tmpz(5:Nx,2:Ny,2:Nz)+tmpz(5:Nx,1:Nym1,2:Nz)) &
				+  Cexhy(5:Nx,2:Ny,2:Nz)*(tmpy(5:Nx,2:Ny,2:Nz)-tmpy(5:Nx,2:Ny,1:Nzm1)) &
				+  Cexhx(5:Nx,2:Ny,2:Nz)*(tmpy(5:Nx,2:Ny,2:Nz)+tmpy(1:Nxm1-3,2:Ny,2:Nz)) &
				-  Cexhi(5:Nx,2:Ny,2:Nz)*(-tmpy(5:Nx,2:Ny,2:Nz)-tmpy(1:Nxm1-3,2:Ny,2:Nz)));
      
	y(2:Nx,5:Ny,2:Nz,2) = x(2:Nx,5:Ny,2:Nz,2) &
				- (Ceyhx(2:Nx,5:Ny,2:Nz)*(-tmpx(2:Nx,5:Ny,2:Nz)+tmpx(2:Nx,5:Ny,1:Nzm1)) &
				+  Ceyhz(2:Nx,5:Ny,2:Nz)*(tmpz(2:Nx,5:Ny,2:Nz)-tmpz(1:Nxm1,5:Ny,2:Nz)) &
				+  Ceyhy(2:Nx,5:Ny,2:Nz)*(tmpz(2:Nx,5:Ny,2:Nz)+tmpz(2:Nx,1:Nym1-3,2:Nz)) &
				-  Ceyhi(2:Nx,5:Ny,2:Nz)*(-tmpz(2:Nx,5:Ny,2:Nz)-tmpz(2:Nx,1:Nym1-3,2:Nz)));
 
	y(2:Nx,2:Ny,5:Nz,3) = x(2:Nx,2:Ny,5:Nz,3) &
				- (Cezhy(2:Nx,2:Ny,5:Nz)*(-tmpy(2:Nx,2:Ny,5:Nz)+tmpy(1:Nxm1,2:Ny,5:Nz)) &
				+  Cezhx(2:Nx,2:Ny,5:Nz)*(tmpx(2:Nx,2:Ny,5:Nz)-tmpx(2:Nx,1:Nym1,5:Nz)) &
				+  Cezhz(2:Nx,2:Ny,5:Nz)*(tmpx(2:Nx,2:Ny,5:Nz)+tmpx(2:Nx,2:Ny,1:Nzm1-3)) &
				-  Cezhi(2:Nx,2:Ny,5:Nz)*(-tmpx(2:Nx,2:Ny,5:Nz)-tmpx(2:Nx,2:Ny,1:Nzm1-3)));

    
end subroutine matvec
