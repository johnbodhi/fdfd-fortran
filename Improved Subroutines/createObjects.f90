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
                        se(mi,mj,mk)    = electric_conductivity(r);
                        sm(mi,mj,mk)    = magnetic_conductivity(r);
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
                        se(mi,mj,mk)    = electric_conductivity(r);
                        sm(mi,mj,mk)    = magnetic_conductivity(r);
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
                        se(mi,mj,mk)    = electric_conductivity(r);
                        sm(mi,mj,mk)    = magnetic_conductivity(r);
					end if
				end do
			end do
		end do
    end do
    
	
    se_xi(2:Nxm1, 2:Nym1, 2:Nzm1) = (se(2:Nxm1, 2:Nym1, 2:Nzm1) + se(2:Nxm1, 1:Nym2, 2:Nzm1) &
								+ se(2:Nxm1, 1:Nym2, 1:Nzm2) + se(2:Nxm1, 2:Nym1, 1:Nzm2))/4;
	se_xy = se_xi;
	se_xz = se_xi;

	se_yi(2:Nxm1, 2:Nym1, 2:Nzm1) = (se(2:Nxm1, 2:Nym1, 2:Nzm1) + se(1:Nxm2, 2:Nym1, 2:Nzm1) &
								+ se(2:Nxm1, 2:Nym1, 1:Nzm2) + se(1:Nxm2, 2:Nym1, 1:Nzm2))/4;
	se_yz = se_yi;
	se_yx = se_yi;

	se_zi(2:Nxm1, 2:Nym1, 2:Nzm1) = (se(2:Nxm1, 2:Nym1, 2:Nzm1) + se(1:Nxm2, 2:Nym1, 2:Nzm1) &
								+ se(2:Nxm1, 1:Nym2, 2:Nzm1) + se(1:Nxm2, 1:Nym2, 2:Nzm1))/4;
	se_zx = se_zi;
	se_zy = se_zi;
    
    
    sm_xi(2:Nxm1, 2:Nym1, 2:Nzm1) = (sm(2:Nxm1, 2:Nym1, 2:Nzm1) + sm(1:Nxm2, 2:Nym1, 2:Nzm1))/2;
	sm_xy = sm_xi; 
	sm_xz = sm_xi; 

	sm_yi(2:Nxm1, 2:Nym1, 2:Nzm1) = (sm(2:Nxm1, 2:Nym1, 2:Nzm1) + sm(2:Nxm1, 1:Nym2, 2:Nzm1))/2;
	sm_yz = sm_yi; 
	sm_yx = sm_yi; 

	sm_zi(2:Nxm1, 2:Nym1, 2:Nzm1) = (sm(2:Nxm1, 2:Nym1, 2:Nzm1) + sm(2:Nxm1, 2:Nym1, 1:Nzm2))/2;
	sm_zx = sm_zi; 
	sm_zy = sm_zi;
    
    

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