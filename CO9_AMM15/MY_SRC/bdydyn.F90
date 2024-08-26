MODULE bdydyn
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to velocities
   !!======================================================================
   !! History :  1.0  !  2005-02  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-07  (D. Storkey) Move Flather implementation to separate routine.
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.2  !  2008-04  (R. Benshila) consider velocity instead of transport 
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!   bdy_dyn        : split velocities into barotropic and baroclinic parts
   !!                    and call bdy_dyn2d and bdy_dyn3d to apply boundary
   !!                    conditions
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE bdydyn2d        ! open boundary conditions for barotropic solution
   USE bdydyn3d        ! open boundary conditions for baroclinic velocities
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  !
   USE domvvl          ! variable volume
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn    ! routine called in dyn_nxt

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn( kt, dyn3d_only )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn  ***
      !!
      !! ** Purpose : - Wrapper routine for bdy_dyn2d and bdy_dyn3d.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           ::   kt           ! Main time step counter
      LOGICAL, INTENT(in), OPTIONAL ::   dyn3d_only   ! T => only update baroclinic velocities
      !
      INTEGER ::   jk, ii, ij, ib_bdy, ib, igrd     ! Loop counter
      LOGICAL ::   ll_dyn2d, ll_dyn3d, ll_orlanski
      REAL(wp) ::  r1_2dt   !  reciprical time step
      REAL(wp), ALLOCATABLE,       DIMENSION(:,:,:) ::   ua_tmp, va_tmp     ! RDP: temporary save of 2d momentum
      REAL(wp), DIMENSION(jpi,jpj) :: pua2d, pva2d     ! after barotropic velocities
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zbdytrdu, zbdytrdv  ! bdy trends
      !!----------------------------------------------------------------------
      !
      ll_dyn2d = .true.
      ll_dyn3d = .true.
      !
      IF( PRESENT(dyn3d_only) ) THEN
         IF( dyn3d_only )   ll_dyn2d = .false.
      ENDIF
      !
      ! RDP
      IF ( l_trddyn ) THEN
         r1_2dt = 1._wp / (2. * rdt)        ! Euler or leap-frog time step 
         IF( neuler == 0 .AND. kt == nit000 )   r1_2dt = 1._wp / rdt
         IF ( ll_dyn3d ) THEN
             ALLOCATE( zbdytrdu(jpi,jpj,jpk), zbdytrdv(jpi,jpj,jpk) )
             zbdytrdu(:,:) = 0._wp
             zbdytrdv(:,:) = 0._wp
         ENDIF
      ENDIF
      ! RDP END

      ll_orlanski = .false.
      DO ib_bdy = 1, nb_bdy
         IF ( cn_dyn2d(ib_bdy) == 'orlanski' .OR. cn_dyn2d(ib_bdy) == 'orlanski_npo' &
     &   .OR. cn_dyn3d(ib_bdy) == 'orlanski' .OR. cn_dyn3d(ib_bdy) == 'orlanski_npo')   ll_orlanski = .true.
      END DO

      ! RDP
      IF ( l_trddyn ) THEN
         ALLOCATE( ua_tmp(jpi,jpj,jpk), va_tmp(jpi,jpj,jpk) )
         ua_tmp(:,:,:) = ua(:,:,:)
         va_tmp(:,:,:) = va(:,:,:)
      ENDIF
      !RDP END

      !-------------------------------------------------------
      ! Split velocities into barotropic and baroclinic parts
      !-------------------------------------------------------

      !                          ! "After" velocities: 
      pua2d(:,:) = 0._wp
      pva2d(:,:) = 0._wp     
      DO jk = 1, jpkm1
         pua2d(:,:) = pua2d(:,:) + e3u_a(:,:,jk) * ua(:,:,jk) * umask(:,:,jk)
         pva2d(:,:) = pva2d(:,:) + e3v_a(:,:,jk) * va(:,:,jk) * vmask(:,:,jk)
      END DO
      pua2d(:,:) = pua2d(:,:) * r1_hu_a(:,:)
      pva2d(:,:) = pva2d(:,:) * r1_hv_a(:,:)

      DO jk = 1 , jpkm1
         ua(:,:,jk) = ( ua(:,:,jk) - pua2d(:,:) ) * umask(:,:,jk)
         va(:,:,jk) = ( va(:,:,jk) - pva2d(:,:) ) * vmask(:,:,jk)
      END DO


      IF( ll_orlanski ) THEN     ! "Before" velocities (Orlanski condition only) 
         DO jk = 1 , jpkm1
            ub(:,:,jk) = ( ub(:,:,jk) - ub_b(:,:) ) * umask(:,:,jk)
            vb(:,:,jk) = ( vb(:,:,jk) - vb_b(:,:) ) * vmask(:,:,jk)
         END DO
      ENDIF

      !-------------------------------------------------------
      ! Apply boundary conditions to barotropic and baroclinic
      ! parts separately
      !-------------------------------------------------------

      IF( ll_dyn2d )   CALL bdy_dyn2d( kt, pua2d, pva2d, ub_b, vb_b, r1_hu_a(:,:), r1_hv_a(:,:), ssha )

      IF( ll_dyn3d )   CALL bdy_dyn3d( kt )

      !-------------------------------------------------------
      ! Recombine velocities
      !-------------------------------------------------------
      !
      DO jk = 1 , jpkm1
         ua(:,:,jk) = ( ua(:,:,jk) + pua2d(:,:) ) * umask(:,:,jk)
         va(:,:,jk) = ( va(:,:,jk) + pva2d(:,:) ) * vmask(:,:,jk)
      END DO
      !
      IF ( ll_orlanski ) THEN
         DO jk = 1 , jpkm1
            ub(:,:,jk) = ( ub(:,:,jk) + ub_b(:,:) ) * umask(:,:,jk)
            vb(:,:,jk) = ( vb(:,:,jk) + vb_b(:,:) ) * vmask(:,:,jk)
         END DO
      END IF
      !

      IF ( l_trddyn ) THEN
         zbdytrdu(:,:,:) = ( ua(:,:,:) - ua_tmp(:,:,:) ) * r1_2dt
         zbdytrdv(:,:,:) = ( va(:,:,:) - va_tmp(:,:,:) ) * r1_2dt
         CALL trd_dyn( zbdytrdu, zbdytrdv, jpdyn_bdy, kt )
         DEALLOCATE( ua_tmp, va_tmp, zbdytrdu, zbdytrdv )
      ENDIF

   END SUBROUTINE bdy_dyn

   !!======================================================================
END MODULE bdydyn
