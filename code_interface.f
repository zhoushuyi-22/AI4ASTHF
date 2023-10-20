!-----------------------------------
      subroutine sfc_ocean                                              &
!...................................
!  ---  inputs:
     &     ( im, ps, u1, v1, t1, q1, tskin, cm, ch, rcl,                &
     &       prsl1, prslki, slimsk, ddvel, flag_iter,                   &
!  ---  outputs:
     &       qsurf, cmm, chh, evap, hflx                                &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_ocean                                                     !
!       inputs:                                                         !
!          ( im, ps, u1, v1, t1, q1, tskin, cm, ch, rcl,                !
!            prsl1, prslki, slimsk, ddvel, flag_iter,                   !
!       outputs:                                                        !
!            qsurf, cmm, chh, evap, hflx )                              !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!         xxxx  --             created                                  !
!    oct  2006  -- h. wei      modified (need description)              !
!    apr  2009  -- y.-t. hou   modified to match the modified gbphys_v.f!
!                     rmoved unused variable from argument list.        !
!                     reformatted code and added program documentation. !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     ps       - real, surface pressure                            im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     rcl      - real,                                             im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     slimsk   - real, sea/land/ice mask (=0/1/2)                  im   !
!     ddvel    - real,                                             im   !
!     flag_iter- logical,                                          im   !
!                                                                       !
!  outputs:                                                             !
!     qsurf    - real, specific humidity at sfc                    im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only : cp => con_cp, rd => con_rd, eps => con_eps,  &
     &                     epsm1 => con_epsm1, hvap => con_hvap,        &
     &                     rvrdm1 => con_fvirt
      use netcdf
!
      implicit none
!
!  ---  constant parameters:
      real (kind=kind_phys), parameter :: cpinv  = 1.0/cp
      real (kind=kind_phys), parameter :: hvapi  = 1.0/hvap
      real (kind=kind_phys), parameter :: elocp  = hvap/cp

!  ---  inputs:
      integer, intent(in) :: im

      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &      t1, q1, tskin, cm, ch, rcl, prsl1, prslki, slimsk, ddvel

      logical, intent(in) :: flag_iter(im)

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: qsurf,       &
     &       cmm, chh, evap, hflx

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: psurf, ps1, q0, qss,      &
     &       rch, rho, theta1, tv1, xrcl, wind

      real (kind=kind_phys) :: tem

      integer :: i

      logical :: flag(im)
!  --------------------------Shuyi Zhou-----------------------------
      INTEGER :: fidA, status, w1_hs_id, w2_hs_id, w3_hs_id, w4_hs_id
      INTEGER :: b1_hs_id, b2_hs_id, b3_hs_id, b4_hs_id, input_id

      INTEGER :: fidB, w1_hl_id, w2_hl_id, w3_hl_id, w4_hl_id
      INTEGER :: b1_hl_id, b2_hl_id, b3_hl_id, b4_hl_id

      CHARACTER(LEN=100) :: file_in_hs
      CHARACTER(LEN=100) :: file_in_hl

      REAL*4 :: b1_hs(1,32), b2_hs(1,16), b3_hs(1,8), b4_hs(1,1)
      REAL*4 :: input(1,5), w1_hs(5,32), w2_hs(32,16), w3_hs(16,8) 
      REAL*4 :: w4_hs(8,1), lh_ncar(1,1), sh_ncar(1,1)
      REAL*4 :: wb1_hs(1,32),wb2_hs(1,16) ,wb3_hs(1,8) ,wb4_hs(1,1) 

      REAL*4 :: b1_hl(1,32), b2_hl(1,16), b3_hl(1,8), b4_hl(1,1)
      REAL*4 :: w1_hl(5,32), w2_hl(32,16), w3_hl(16,8), w4_hl(8,1)
      REAL*4 :: wb1_hl(1,32),wb2_hl(1,16) ,wb3_hl(1,8) ,wb4_hl(1,1) 

      file_in_hs='/GFPS8p/cess2/zsy/fortran_test/new/LHF_pinn.nc'
      file_in_hl='/GFPS8p/cess2/zsy/fortran_test/new/LHF_pinn.nc' 
      status = NF90_OPEN(TRIM(file_in_hs),0,fidA)
      if (status /= nf90_noerr) print *, '1t',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"w1",w1_hs_id)
      if (status /= nf90_noerr) print *, '2t',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"b1",b1_hs_id)
      if (status /= nf90_noerr) print *, '3t',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"w2",w2_hs_id)
      status = NF90_INQ_VARID(fidA,"b2",b2_hs_id)
      status = NF90_INQ_VARID(fidA,"w3",w3_hs_id)
      status = NF90_INQ_VARID(fidA,"b3",b3_hs_id)
      status = NF90_INQ_VARID(fidA,"w4",w4_hs_id)
      status = NF90_INQ_VARID(fidA,"b4",b4_hs_id)

      ! Read hs variable values
      status = NF90_GET_VAR(fidA,w1_hs_id,w1_hs)
      status = NF90_GET_VAR(fidA,b1_hs_id,b1_hs)
      status = NF90_GET_VAR(fidA,w2_hs_id,w2_hs)
      status = NF90_GET_VAR(fidA,b2_hs_id,b2_hs)
      status = NF90_GET_VAR(fidA,w3_hs_id,w3_hs)
      status = NF90_GET_VAR(fidA,b3_hs_id,b3_hs)
      status = NF90_GET_VAR(fidA,w4_hs_id,w4_hs)
      status = NF90_GET_VAR(fidA,b4_hs_id,b4_hs)

      status = NF90_CLOSE(fidA)
  
      ! Read hl dimension IDs
      status = NF90_OPEN(TRIM(file_in_hl),0,fidB)

      ! Read hs dimension IDs
      status = NF90_INQ_VARID(fidB,"w1",w1_hl_id)
      status = NF90_INQ_VARID(fidB,"b1",b1_hl_id)
      status = NF90_INQ_VARID(fidB,"w2",w2_hl_id)
      status = NF90_INQ_VARID(fidB,"b2",b2_hl_id)
      status = NF90_INQ_VARID(fidB,"w3",w3_hl_id)
      status = NF90_INQ_VARID(fidB,"b3",b3_hl_id)
      status = NF90_INQ_VARID(fidB,"w4",w4_hl_id)
      status = NF90_INQ_VARID(fidB,"b4",b4_hl_id)

 
      ! Read hs variable values
      status = NF90_GET_VAR(fidB,w1_hl_id,w1_hl)
      status = NF90_GET_VAR(fidB,b1_hl_id,b1_hl)
      status = NF90_GET_VAR(fidB,w2_hl_id,w2_hl)
      status = NF90_GET_VAR(fidB,b2_hl_id,b2_hl)
      status = NF90_GET_VAR(fidB,w3_hl_id,w3_hl)
      status = NF90_GET_VAR(fidB,b3_hl_id,b3_hl)
      status = NF90_GET_VAR(fidB,w4_hl_id,w4_hl)
      status = NF90_GET_VAR(fidB,b4_hl_id,b4_hl)

      status = NF90_CLOSE(fidB)

      !print *, "wwww", w1_hs(1,1:3)
!
!===> ...  begin here
!
!  --- ...  flag for open water
      do i = 1, im
         flag(i) = ( slimsk(i)==0.0 .and. flag_iter(i) )
      enddo

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specified
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qss is sat. hum. at surface

      do i = 1, im
        if ( flag(i) ) then
          xrcl(i)  = sqrt( rcl(i) )
          psurf(i) = 1000.0 * ps(i)
          ps1(i)   = 1000.0 * prsl1(i)
          theta1(i) = t1(i) * prslki(i)

          wind(i)  = xrcl(i) * sqrt(u1(i)*u1(i) + v1(i)*v1(i))          &
     &             + max( 0.0, min( ddvel(i), 30.0 ) )
          wind(i)  = max( wind(i), 1.0 )

          q0(i) = max( q1(i), 1.0e-8 )
          tv1(i) = t1(i) * (1.0 + rvrdm1*q0(i))
          rho(i) = ps1(i) / (rd*tv1(i))

          qss(i) = fpvs( tskin(i) )
          qss(i) = eps*qss(i) / (psurf(i) + epsm1*qss(i))
        endif
      enddo

      do i = 1, im
        if ( flag(i) ) then
          evap(i) = 0.0
          hflx(i) = 0.0
        endif
      enddo

!  --- ...  rcp = rho cp ch v

      do i = 1, im
        if ( flag(i) ) then
          rch(i) = rho(i) * cp * ch(i) * wind(i)
          cmm(i) = cm(i) * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
        endif
      enddo

!  --- ...  sensible and latent heat flux over open water

      do i = 1, im
        if ( flag(i) ) then
          hflx(i) = rch(i) * (tskin(i) - theta1(i))
          evap(i) = elocp*rch(i) * (qss(i) - q0(i))
          qsurf(i) = qss(i)
        endif
      enddo

      do i = 1, im
        if ( flag(i) ) then
          tem     = 1.0 / rho(i)
!  -----------------Shuyi Zhou-----------------------------

          input(1,1) = wind(i)
          input(1,2) = tskin(i)-273.15
          input(1,3) = t1(i)-273.15
          input(1,4) = qss(i)*1000
          input(1,5) = q0(i)*1000
          
          wb1_hs = matmul(input, w1_hs)+b1_hs
          where (wb1_hs<0) wb1_hs = 0.3*wb1_hs
          wb2_hs = matmul(wb1_hs, w2_hs)+b2_hs
          where (wb2_hs<0) wb2_hs = 0.3*wb2_hs
          wb3_hs = matmul(wb2_hs, w3_hs)+b3_hs
          where (wb3_hs<0) wb3_hs = 0.3*wb3_hs
          wb4_hs = matmul(wb3_hs, w4_hs)+b4_hs

          wb1_hl = matmul(input, w1_hl)+b1_hl
          where (wb1_hl<0) wb1_hl = 0.3*wb1_hl
          wb2_hl = matmul(wb1_hl, w2_hl)+b2_hl
          where (wb2_hl<0) wb2_hl = 0.3*wb2_hl
          wb3_hl = matmul(wb2_hl, w3_hl)+b3_hl
          where (wb3_hl<0) wb3_hl = 0.3*wb3_hl
          wb4_hl = matmul(wb3_hl, w4_hl)+b4_hl

          sh_ncar(1,1) = hflx(i) 
          lh_ncar(1,1) = evap(i)
          
          hflx(i) = wb4_hs(1,1) * tem * cpinv
          evap(i) = wb4_hl(1,1) * tem * hvapi
!  ------------------------------------------------------------
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean
!-----------------------------------
