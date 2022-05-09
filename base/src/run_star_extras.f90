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
         real(dp) :: LEdd, Lmax, L_rel_err
         real(dp), parameter :: rel_err_tol = 1d-3

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         
         ! print v_flag
         ! write(*,*) "v_flag:", s%v_flag, "u_flag:", s%u_flag


         ! if (s% x_logical_ctrl(1)) then
         !    ! check if convective above a certain column depth logy=x_ctrl(1)
            ! allocate(ycol(s%nz))
            ! ycol = s%P/(10**s%log_surface_gravity) ! y=P/g hydrostatic equilibrium
         !    if ( ANY( s% mixing_type == convective_mixing .and. ycol .gt. 10**s%x_ctrl(1) ) ) then
         !          extras_check_model = terminate 
         !          s% termination_code = t_xtra1
         !          termination_code_str(t_xtra1) = 'Convective!'
         !    end if      
         
         if (s% x_logical_ctrl(1)) then
            ! check if eps_nuc larger in He layer than H layer
            ! Boundary between layers is given by (Cumming and Bildsten 2000):
            X = s%accretion_species_xa(1)
            Z = s%accretion_species_xa(3)
            R = s%r(1)
            Mdot_edd = 2*mp*clight/(1+X)/R/sige * 4*pi*R**2 / Msun * secyer
            yd = 6.8d8 * s%mass_change/(0.1*Mdot_edd) * 0.01/Z * X/0.71

            allocate(ycol(s%nz))
            ycol = s%P/(10**s%log_surface_gravity) ! y=P/g hydrostatic equilibrium

            if (ycol(s%max_eps_nuc_k) > yd) then
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
            LEdd = pi4 * clight * s%cgrav(1) * s%m_grav(1) / (0.2 * (1 + s%X(1)))
            Lmax = s%x_ctrl(2)*LEdd

            L_rel_err = (Lrad(1) - Lmax)/Lmax

            if (abs(L_rel_err) .lt. 0.25) then
               write(*,*) "Lrad/LEdd=", Lrad(1)/LEdd
            end if

            if (abs(L_rel_err) .lt. rel_err_tol) then
               extras_check_model = terminate
               s% termination_code = t_xtra2
               termination_code_str(t_xtra2) = 'Eddington!'
            else if (L_rel_err > 0) then
               extras_check_model = redo
               s%dt = 0.5d0 * s%dt 
               write(*,*) "Above Eddington tolerance, redoing with 1/2 timestep (dt=", s%dt, ")"  
            endif


         else if (s% x_logical_ctrl(3)) then
            ! check if timestep has dropped below x_ctrl(3) (years)
            ! but only if star_age > x_ctrl(4) (seconds)
            ! write(*,*) s%dt, s%dt/secyer
            ! s%dt is in sec but s%star_age is in yrs!
            if (s%dt/secyer < s%x_ctrl(3) .and. s%star_age*secyer > s%x_ctrl(4)) then
               extras_check_model = terminate
               s% termination_code = t_xtra3
               termination_code_str(t_xtra3) = 'Critical timestep reached!'
            endif 


         end if


         ! indicate where MESA terminated
         ! if (extras_check_model == terminate) s% termination_code = t_extras_check_model
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

         ! This does not work
         ! do while(minval(s%rho) < density)
         !    ! write(*,*) 'min rho = ', minval(s%rho)
         !    do k=1, s%nz
         !       if (s%rho(k) < density) then
         !          write(*,2) 'do_remove_surface (v3)', k, s%rho(k), density
                  ! write(*,*) 'do_remove_surface (v2) - model #',s%model_number,' - remove at index',k,' - density ',s% rho(k),' - cutoff ',density
         !          call do_remove_surface(id, k+1, ierr)
         !       endif 
         !    end do
         ! end do
         ! return

         ! Make it so that do_remove_surface is only called once
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
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

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

      end module run_star_extras
      
