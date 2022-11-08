! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none

      include "test_suite_extras_def.inc"      
      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"
      ! include "other_cgrav.inc"   
      
      !include 'standard_run_star_extras.inc'

      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         ! Custom routines
         s% other_kap_get => my_other_kap_get
         s% other_remove_surface => my_other_remove_surface

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr, k
         type (star_info), pointer :: s

         real(dp) :: X,Z,R,Mdot_edd,yd
         real(dp), dimension(:), allocatable :: Lrad,  ycol
         real(dp) :: LEdd, Lmax, L_over_LEdd, L_rel_err
         real(dp), parameter :: rel_err_tol = 1d-3
         logical :: file_exists
         integer :: percent_save
         real(dp) :: frac_save
         character (len=strlen) :: frac_save_str,save_model_filename
         integer :: k_sonic, k_phot

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         


         ! Various stopping conditions depending on which part of the simulation we are at
         
         if (s% x_logical_ctrl(1)) then
            ! check if eps_nuc larger in He layer than H layer, and if it's convective. This indicates "ignition"
            ! Boundary between layers is given by (Cumming and Bildsten 2000):
            X = s%accretion_species_xa(1)
            Z = s%accretion_species_xa(3)
            R = s%r(1)
            Mdot_edd = 2*mp*clight/(1+X)/R/sige * 4*pi*R**2 / Msun * secyer
            yd = 6.8d8 * s%mass_change/(0.1*Mdot_edd) * 0.01/Z * X/0.71

            allocate(ycol(s%nz))
            ycol = s%P/(10**s%log_surface_gravity) ! y=P/g hydrostatic equilibrium

            !if (ycol(s%max_eps_nuc_k) > yd) then
            if ( ANY( ycol(s%max_eps_nuc_k) .gt. yd .and. s% mixing_type == convective_mixing ) ) then
               extras_check_model = terminate
               s% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'Ignited!'
            end if

         else if (s% x_logical_ctrl(2)) then
            ! check if within 1e-3 relative tolerate of Lmax = x_ctrl(2)*LEdd
            ! at the surface. Retry with small timestep if overshoot

            ! radiative luminosity (see get_Lrad in star/private/star_utils.f90)
            allocate(Lrad(s%nz))
            Lrad = s%L - s%L_conv

            !LEdd = 4*pi*standard_cgrav*Msun*1.4*clight/s% opacity(1) 
            ! LEdd = pi4 * clight * s%cgrav(1) * s%m_grav(1) / s%opacity(1) ! ~same as line above, cgrav=G unless GR factors on, m_grav is mass_enclosed, might be slightly larger than 1.4
            LEdd = pi4 * clight * s%cgrav(1) * s%m_grav(1) / (0.2 * (1 + s%X(1)))  ! Ledd at infinity, so no T corrections
            L_over_LEdd = Lrad(1)/LEdd

            Lmax = s%x_ctrl(2)*LEdd
            L_rel_err = (Lrad(1) - Lmax)/Lmax

            if (abs(L_rel_err) .lt. 0.50) then
               write(*,*) "Lrad/LEdd=", L_over_LEdd
            end if
            
            !do percent_save = 50,90,10
            do percent_save = INT(s%x_ctrl(2)*100 - 50), 90, 10
               frac_save = real(percent_save)/100
               write(frac_save_str,'(f4.2)') frac_save
               save_model_filename = 'models/ns_env_'//trim(frac_save_str)//'Edd.mod'
               inquire(FILE=save_model_filename, EXIST=file_exists)
               if (L_over_LEdd .gt. frac_save .and. L_over_LEdd .lt. frac_save+0.05 .and. .not. file_exists) then
                  call star_write_model(s% id, save_model_filename, ierr)
                  write(*,*) 'saved to ', save_model_filename

                  ! force a profile save
                  s% need_to_save_profiles_now = .true.
                  s% save_profiles_model_priority = 10
               end if
            end do

            if (abs(L_rel_err) .lt. rel_err_tol) then
               extras_check_model = terminate
               s% termination_code = t_xtra2
               termination_code_str(t_xtra2) = 'Eddington!'
            else if (L_rel_err > 0) then
               extras_check_model = redo
               s%dt = 0.5d0 * s%dt 
               write(*,*) "Above Eddington tolerance, redoing with 1/2 timestep (dt=", s%dt, ")"  
            end if


         else if (s% x_logical_ctrl(3)) then
            ! check if timestep has dropped below x_ctrl(3) (years)
            ! but only if star_age > x_ctrl(4) (seconds)
            ! **  s%dt is in sec but s%star_age is in yrs!
            if (s%dt/secyer < s%x_ctrl(3) .and. s%star_age*secyer > s%x_ctrl(4)) then
               extras_check_model = terminate
               s% termination_code = t_xtra3
               termination_code_str(t_xtra3) = 'Critical timestep reached!'
            end if 


         ! else if (s% x_logical_ctrl(4)) then
         !    call smooth_outer_abundances(id)
         ! end if

         else if (s% x_logical_ctrl(4)) then
            ! Check if sonic point is past the photosphere
            ! but only if star_age > x_ctrl(4) (seconds)
            if (s%star_age*secyer > s%x_ctrl(4)) then 

               k_sonic = minloc(abs(s%v_div_csound - 1), 1)
               k_phot = minloc(abs( s%L/(4*pi*(s%r)**2) - boltz_sigma*(s%T)**4), 1)
<<<<<<< Updated upstream
=======

               write(*,*) "Sonic pt: k=",k_sonic," r(k)=", s%r(k_sonic)
               write(*,*) "Phot: k=",k_phot," r(k)=", s%r(k_phot)
>>>>>>> Stashed changes
               if (k_phot > k_sonic - 100) then
                  write(*,*) "ksonic ", k_sonic, "kphot ", k_phot
               end if
               
               if (k_phot > k_sonic) then
                  write(*,*) "Sonic point: k=", k_phot, ", r=", s%r(k_sonic)/1d5, " km, rho=", s%rho(k_sonic), " g/cm3"
                  write(*,*) "Photosphere: k=", k_phot, ", r=", s%r(k_phot)/1d5, " km, rho=", s%rho(k_phot), " g/cm3"
                  extras_check_model = terminate
                  s% termination_code = t_xtra4
                  termination_code_str(t_xtra4) = 'Sonic point past the photosphere'
               end if


            end if
         end if

      end function extras_check_model



      subroutine my_other_remove_surface(id, ierr, k)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer, intent(out) :: k

         type (star_info), pointer :: s
         real(dp) :: min_density
         integer :: j,j_remove

         ierr = 0
         include 'formats'
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_density: get_star_ptr ierr', ierr
            return
         end if

         min_density = s%x_ctrl(5)
         
         ! remove everything below any point that goes below minimum density
         if (minval(s%rho, mask=s%rho>0) > min_density) return
         do j=1, s%nz
            if (s%rho(j) < min_density) then
               j_remove = j
            end if
            if (s%rho(j) > s%x_ctrl(6)) then  ! stop searching at some density
               exit
            end if
         end do
         ! call do_remove_surface(id, k_remove, ierr) ! don't need to call, will be done automatically with output k

         ! write(*,*) 'do_remove_surface (custom) - model #',s%model_number,' - remove at index',j_remove,' - density (k-1)',s% rho(j_remove-1),' - cutoff ', min_density
         k = j_remove ! The cell to remove down to.
      end subroutine my_other_remove_surface


      ! Smoothing abundance curves
      ! This is copied code from star/private/adjust_mass.f90
      ! But there it's only called if the mass changes, i.e. during accretion
      ! Putting it here  so it can be called at will   
      subroutine smooth_outer_abundances(id)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr,j,k,m,l
         real(dp) :: partial_xa_mass, region_total_mass, frac
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% smooth_outer_xa_big > 0.0d0 .and. s% smooth_outer_xa_small > 0.0d0) then
            write(*,*) 'doing smooth_outer_abundance (run_star_extras)'
            write(*,*) 'big: ', s% smooth_outer_xa_big, 'small: ', s% smooth_outer_xa_small
            write(*,*) 'total mass: ', s% m(1)
            write(*,*) 'f_big * m = ', s% smooth_outer_xa_big *s% m(1)
            write(*,*) 'f_small * m = ', s% smooth_outer_xa_small *s% m(1)
            write(*,*) ''
            write(*,*) 'm and dm at indexes: 224, 238'
            write(*,*) s%m(224), s%dm(224)
            write(*,*) s%m(238), s%dm(238)
            write(*,*) ''
            m = 1
            do k = 1, s% nz
               if (s% m(1) - s% m(k) > s% smooth_outer_xa_big * s% m(1)) exit
               region_total_mass = 0
               m = k
               do l = k, s% nz
                  if (s% m(k) - s% m(l) > s% smooth_outer_xa_small * s% m(1) * &
                     (1-(s% m(1) - s% m(k))/(s% smooth_outer_xa_big * s% m(1)))) exit
                  m = l
                  region_total_mass = region_total_mass + s% dm(l)
               end do
               write(*,*) "region_total_mass: ", region_total_mass
               write(*,*) "k,m : ", k, m
               write(*,*) ""
               if (m == k) cycle
               do j=1,s% species
                  partial_xa_mass = 0
                  do l = k, m
                     partial_xa_mass = partial_xa_mass + s% dm(l) * s% xa(j,l)
                  end do
                  do l = k, m
                     s% xa(j,l) = partial_xa_mass / region_total_mass
                  end do
               end do
            end do
            ! fix small errors to ensure xa's add up to unity
            do k=1,m
               frac = 1d0/sum(s% xa(1:s% species,k))
               do j=1,s% species
                  s% xa(j,k) = s% xa(j,k)*frac
               end do
            end do
            ! call star_write_model(s% id, "models/foo.mod", ierr)
         end if
      end subroutine smooth_outer_abundances


      subroutine my_other_kap_get( &
         id, k, handle, zbar, X, Z, XC, XN, XO, XNe, &
         log10_rho, log10_T, species, chem_id, net_iso, xa, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

         use chem_def
         use kap_lib
      
         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ion charge
         real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition    
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons

         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
      
         integer i, Zi, Ai
         ! real(dp) :: alpha, xx, T, rho, Ye, YZ2, kap1, kap2, eps
         real(dp) :: alpha, xx, T, rho, Ye, kap1, kap2, eps
         real(dp) kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT
         real(dp) kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         
         ! OUTPUT
         real(dp), intent(out) :: frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
            
         type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)

               
            
         ! Calculate Ye and YZ2
         rho = 10.0**log10_rho
         T = 10.0**log10_T
         Ye = 0.0
         ! YZ2 = 0.0
         do i = 1,species
            Zi = chem_isos%Z(chem_id(i))
            Ai = chem_isos%N(chem_id(i)) + Zi
            Ye = Ye + xa(i) * Zi / Ai
            ! YZ2 = YZ2 + xa(i) * Zi * Zi / Ai         	
            !write (*,'(i5,i5,g12.5,i5,i5,i5,i5)') species,i,xa(i),chem_id(i),net_iso(chem_id(i)),chem_isos%Z(chem_id(i)),&
            !		chem_isos%N(chem_id(i))+chem_isos%Z(chem_id(i))
         end do

         frac_Type2 = 1.0; 
         
    
         ! Paczynski electron scattering temperature part only
         ! write(*,*) "using electron scattering only"
         alpha = 0.86
         xx = (T/4.5d8)**alpha       
         kap_rad = 0.4 * Ye / (1.0 + xx) 
         dlnkap_rad_dlnRho = 0; 
         dlnkap_rad_dlnT = -alpha * xx / (1+xx)
      

         ! use the MESA electron conduction
         call kap_get_elect_cond_opacity( &
               zbar, log10_rho, log10_T, &
               kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
            if (ierr /= 0) return
         !write(*,*) 'kap_rad', kap_rad, ' kap_ec', kap_ec

         ! combine radiative and conductive opacities
         kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)
         dln_kap_dlnRho =  (kap/kap_ec) * dlnkap_ec_dlnRho + (kap/kap_rad) * dlnkap_rad_dlnRho
         dln_kap_dlnT =  (kap/kap_ec) * dlnkap_ec_dlnT + (kap/kap_rad) * dlnkap_rad_dlnT
         
      end subroutine my_other_kap_get
      
      



      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 1
         ! how_many_extra_history_columns = 3
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         use num_lib, only: binary_search
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: check,frac,tau00,taup1,v00,vp1,bb
         integer :: cnt, cur_type
         real(dp), dimension(:), allocatable :: burn_type
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.

         !! Luminosity at tau=1
         names(1) = "luminosity_tau1"

         ! Following star/private/report.f90 check_tau()
         check = 1d0 ! optical depth to look for
         k = binary_search(s%nz, s%tau, 0, check)
         if (k < 1 .or. k > s%nz-2) then

            ! Return -1 if can't find tau=1
            vals(1) = -1

         else

            tau00 = s% tau(k)
            taup1 = s% tau(k+1)
            ! tau00 <= check < taup1
            frac = (check - tau00)/(taup1 - tau00)

            if (frac > 1 .or. frac < 0) then
               vals(1) = -1

            else
               v00 = (s% L(k+1) + s% L(k))/2
               vp1 = (s% L(k+2) + s% L(k+1))/2
               vals(1) = (v00 + (vp1 - v00)*frac)/Lsun
            end if

         end if


         !! Number of mixing regions
         ! names(2) = "N_mix_regions"

         ! ! Following star/private/hisory.f90 
         ! cnt = 1
         ! cur_type = s% mixing_type(s% nz)
         ! do k = s% nz -1 ,1,-1
         !    if (cur_type == s% mixing_type(k)) cycle ! cycle means skip the next lines in the loop and start again
         !    cur_type = s% mixing_type(k)
         !    cnt = cnt+1
         ! end do
         ! vals(2) = cnt


         ! !! Number of burning regions
         ! names(3) = "N_burn_regions"
         
         ! ! Same thing. Def of burn_type from set_burn_types() in history.f90
         ! allocate(burn_type(s% nz))
         ! do k=1,s% nz
         !    bb = s% eps_nuc(k) - s% non_nuc_neu(k)
         !    burn_type(k) = int(sign(1d0,bb)*max(0d0,1d0+safe_log10(abs(bb))))
         ! end do
         ! cnt = 1
         ! cur_type = burn_type(s% nz)
         ! do k = s% nz -1 ,1,-1
         !    if (cur_type == burn_type(k)) cycle
         !    cur_type = burn_type(k)
         !    cnt = cnt+1
         ! end do
         ! vals(3) = cnt

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: seconds_save_profile, f
         logical :: do_mixing_history
         integer :: mixing_history_interval
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
   
         seconds_save_profile = s% x_ctrl(7)
         s% xtra(1) = s% star_age*secyer

         f = 1/seconds_save_profile

         ! Evaluates to true if f times the age in seconds has crossed an integer
         ! For instance if x_ctrl(7) = 0.1, then f=10 and we will save at star_age = (0.1,0.2,0.3,...)
         if (floor(f * s% xtra(1)) - floor(f * s% xtra_old(1)) .ne. 0) then
            s% need_to_save_profiles_now = .true.
            s% save_profiles_model_priority = 10
         end if

         do_mixing_history = s% x_logical_ctrl(10)
         mixing_history_interval = s% x_integer_ctrl(10)
         if (do_mixing_history) then
            if (mod(s% model_number, mixing_history_interval) .eq. 0) then
               call write_mixing_zones_history_file(id)
            end if
         end if

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      
      ! The issue with the existing way to output a history of the locations of mixing regions is that you cannot know in advance the number of mixing regions
      ! Too small a number won't be reported either. In these simulations, convective zones can break up into many (hundreds!) of tiny convective zones, hence 
      ! the need for this function, which generates a separate history file which will output the location of every single mixing region when called (as well 
      ! as the model number and star age). The structure is the same as in the original history file, first column is mix_type, second is mix_qtop, then next region and so on..
      ! Also write column depth instead of q for the mass coordinate of the top of the region
      subroutine write_mixing_zones_history_file(id)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr,io,nz,k,cur_type
         real(dp) :: ytop
         character (len=strlen) :: fname, dbl_fmt, int_fmt, txt_fmt
         logical :: file_exists, write_header
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dbl_fmt = s% star_history_dbl_format
         int_fmt = s% star_history_int_format
         txt_fmt = s% star_history_txt_format

         fname = trim(s %log_directory) // '/mixing_history.data'
         inquire(FILE=fname, EXIST=file_exists)
         if (.not. file_exists) then
            open(newunit=io, file=trim(fname), action='write', iostat=ierr)
            write_header = .true.
         else
            open(newunit=io, file=trim(fname), action='write', position='append', iostat=ierr)
            write_header = .false.
         endif

         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            return
         end if

         if (write_header) then            
            write(io, fmt=txt_fmt, advance='no') 'model_number'
            write(io, fmt=txt_fmt, advance='no') 'star_age'
            write(io, fmt=txt_fmt, advance='no') 'lgybot'
            write(io, fmt=txt_fmt, advance='no') 'mix type 1'
            write(io, fmt=txt_fmt, advance='no') 'lgytop 1'
            write(io, fmt=txt_fmt, advance='no') 'mix type 2'
            write(io, fmt=txt_fmt, advance='no') 'lgytop 2'
            write(io, fmt=txt_fmt, advance='no') '...->'
            write(io,*)
         end if

         write(io, fmt=int_fmt, advance='no') s% model_number
         write(io, fmt=dbl_fmt, advance='no') s% star_age
         nz = s% nz

         ! Start with maximum column
         write(io, fmt=dbl_fmt, advance='no') safe_log10(s% xmstar*sum(s% dq(1:nz))/(4*pi*s% r(nz)*s% r(nz)))

         ! Loop
         cur_type = s% mixing_type(nz)
         do k=nz-1,1,-1
            if (cur_type == s% mixing_type(k)) cycle
            write(io, fmt=int_fmt, advance='no') cur_type
            write(io, fmt=dbl_fmt, advance='no') safe_log10(s% xmstar*sum(s% dq(1:k-1))/(4*pi*s% r(k)*s% r(k)))
            cur_type = s% mixing_type(k)
         end do
         write(io, fmt=int_fmt, advance='no') cur_type  ! The last entry is the convective type that goes to the top of the grid
         write(io,*)

         close(io)

      end subroutine write_mixing_zones_history_file


      end module run_star_extras
      
