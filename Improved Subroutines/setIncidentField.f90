subroutine setIncidentField
	use global
	implicit none
	integer indx,indy,indz;

! traveling in the + z direction
	forall (indx=1:Nx,indy=1:Ny,indz=1:Nz) Einc_x(indx,indy,indz) = A*exp(-j*ko*indz*dz); 
    forall (indx=1:Nx,indy=1:Ny,indz=1:Nz) Hinc_y(indx,indy,indz) = A*exp(-j*ko*(indz+0.5)*dz); 
    !forall (indx=1:Nx,indy=1:Ny,indz=1:Nz) Einc_y(indx,indy,indz) = A*exp(-j*ko*indz*dz); 
	!forall (indx=1:Nx,indy=1:Ny,indz=1:Nz) Hinc_x(indx,indy,indz) = A*exp(-j*ko*(indz+0.5)*dz); 
    Hinc_y = Hinc_y / nu_o;
    !Hinc_x = Hinc_x / nu_o;

end subroutine setIncidentField