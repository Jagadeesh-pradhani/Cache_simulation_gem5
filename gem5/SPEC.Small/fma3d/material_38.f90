!! MATERIAL(*)%PVAL(1:22) Usage
!!
!! 1:Density 2:Bulk_Ln 3:Bulk_Qd 4:HG_Viscosity 5:HG_Stiffness
!!
!! PVAL 1/2/3 6/7/8 10/20/30 11/21/31 22  25/35/45 32   33   36  17/27 38
!!                  40/50    41/51                                 37
!!  6:  K1    D1    E        E        E1Ln  E_aa  E          Ko        A0
!!  7:        D2    Nu       Nu       E1Qd  E_bb  Nu         Kinf      A1
!!  8:              Lambda   Lambda   E2Ln  E_cc  Lam        Kdec      A2
!!  9:              G        G        E2Qd  Nu_ba G     G    G0   Gftn G
!! 10:  Yield                Yield    Fric  Nu_ca Yield Bulk Ginf Bulk Bulk
!! 11:  Eplas                Eplas    Bulk  Nu_cb Ep_ki F_a  Gdec      Pftn
!! 12:  H                    H        Fcmp  G_ab  Ep_is F_b            2G
!! 13:  Beta                 Beta     Gmod  G_ac  pparm F-c            4G/3
!! 14:                                Gcof  G_bc  rparm C_r            Pfra
!! 15:                                                  C_p
!! 16:  p_exp                p_exp                      V_d
!! 17:  D_exp                D_exp                      V_w
!! 18:                                                  Ltyp
!! 19:                                                  Tmax
!!
      SUBROUTINE MATERIAL_38_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 26-AUG-1995 16:57:32.00
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE tabulated_function_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          STRESS(6),                                                     &
     &          TABLE_LOOK_UP           ! Interogates TABULATED_FUNCTION
      INTEGER                                                                  &
     &          Table_ID
      CHARACTER                                                                &
     &          SURFACE_SHAPE*12
      LOGICAL                                                                  &
     &          SIGN_POSITIVE,SIGN_CHANGE
!!
!! INITIALIZE MATERIAL CONSTANTS.
!!
      MATERIAL(MatID)%PVAL(12) = (2.0 * MATERIAL(MatID)%PVAL(9))
      MATERIAL(MatID)%PVAL(13) = (4.0 * MATERIAL(MatID)%PVAL(9)) / 3.0
!!
!! Check unloading bulk modulus and tabulated function for required properties.
!!
      BULK = MATERIAL(MatID)%PVAL(10)
      Table_ID = NINT (MATERIAL(MatID)%PVAL(11))
      WRITE (MSG1,'(I8)') MATERIAL(MatID)%MatID
      WRITE (MSG2,'(I8)') TABULATED_FUNCTION(Table_ID)%TFID
!!
!! Check for zero unloading bulk modulus BULK.
!!
      IF (BULK .LE. 0.0) THEN
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'MATERIAL_38_INIT.001.00'//                                    &
     &    MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                 &
     &    MSGL//'Has A Zero Or Negative Unloading Bulk Modulus (BULK).'        &
     &    )
      ENDIF
!!
!! Check to see that the first pair of entries in the the tabulated pressure
!! versus bulk strain table are (0.0,0.0).
!!
      IF (TABULATED_FUNCTION(Table_ID)%X(1) .NE. 0.0 .OR.                      &
     &      TABULATED_FUNCTION(Table_ID)%Y(1) .NE. 0.0) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MATERIAL_38_INIT.002.00'//                              &
     &          MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//           &
     &          MSGL//'References Tabulated Function ID:   '//MSG2//           &
     &          MSGL//'That Does Not Start With A (0.0,0.0)-Pair.'             &
     &          )
      ENDIF
!!
!! Check to see that all volumetric strains are either positive or negative.
!!
      SIGN_POSITIVE = TABULATED_FUNCTION(Table_ID)%X(2) .GT. 0.0
      DO i = 3,TABULATED_FUNCTION(Table_ID)%Number_Of_Pairs
        IF (SIGN_POSITIVE) THEN
          SIGN_CHANGE = TABULATED_FUNCTION(Table_ID)%X(i) .LT. 0.0
        ELSE
          SIGN_CHANGE = TABULATED_FUNCTION(Table_ID)%X(i) .GT. 0.0
        ENDIF
        IF (SIGN_CHANGE) THEN
          CALL USER_MESSAGE                                                    &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'MATERIAL_38_INIT.003.00'//                                    &
     &    MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                 &
     &    MSGL//'References Tabulated Function ID:   '//MSG2//                 &
     &    MSGL//'That Has Sign Changes In The Bulk Strain Entries.'            &
     &    )
        ENDIF
      ENDDO
!!
!! Force all pressures and volumetric strains to be positive.
!!
      DO i = 1,TABULATED_FUNCTION(Table_ID)%Number_Of_Pairs
        TABULATED_FUNCTION(Table_ID)%X(i) =                                    &
     &    ABS (TABULATED_FUNCTION(Table_ID)%X(i))
        TABULATED_FUNCTION(Table_ID)%Y(i) =                                    &
     &    ABS (TABULATED_FUNCTION(Table_ID)%Y(i))
      ENDDO
!!
!! Check the slope for eack segment of the pressure versus bulk strain
!! tabulated function to insure that the unloading bulk modulus BULK
!! has a larger slope. (The table has pressure positive (+) in compression
!! and bulk strain positive (+) in compression.) (The table lookup for
!! PRESSURE forces the calculation of .SLOPE(*) in case this is the first
!! reference to a tabulated function.)
!!
      PRESSURE = TABLE_LOOK_UP (Table_ID,0.0)
      DO i = 1,TABULATED_FUNCTION(Table_ID)%Number_Of_Pairs-1
        IF (BULK .LT. TABULATED_FUNCTION(Table_ID)%SLOPE(i)) THEN
          CALL USER_MESSAGE                                                    &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'MATERIAL_38_INIT.004.00'//                                    &
     &    MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                 &
     &    MSGL//'References Tabulated Function ID:   '//MSG2//                 &
     &    MSGL//'That Has A Slope Larger Than '                                &
     &        //'The Unloading Bulk Modulus.'                                  &
     &    )
        ENDIF
      ENDDO
!!
!! Evaluate yield surface profile with respect to the pressure axis.
!!
      A0 = MATERIAL(MatID)%PVAL(6)
      A1 = MATERIAL(MatID)%PVAL(7)
      A2 = MATERIAL(MatID)%PVAL(8)
      A3 = A1*A1-4.0*A0*A2

      IF (A1 .EQ. 0.0) GO TO 30
      IF (ABS(A3)/(A1*A1) .LT. 1.0D-14) GO TO 50

 30     CONTINUE
      IF (A3 .LT. 0.0) GO TO 40
      IF (A3 .EQ. 0.0) GO TO 50
      IF (A3 .GT. 0.0) GO TO 90

 40     CONTINUE
      CALL USER_MESSAGE                                                        &
     & (                                                                       &
     & MSGL//'FATAL'//                                                         &
     & MSGL//'MATERIAL_38_INIT.005.00'//                                       &
     & MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                    &
     & MSGL//'Has Yield Surface Constants A0,A1,A2 That'//                     &
     & MSGL//'Produce An Imaginary Yield Surface. The'//                       &
     & MSGL//'Quantity (A1*A1-4.0*A0*A2) Is Negative.'                         &
     & )

 50     CONTINUE
      IF (A2 .LT. 0.0) GO TO 40
      IF (A2 .EQ. 0.0) GO TO 60
      IF (A2 .GT. 0.0) GO TO 80

 60     CONTINUE
      IF (A0 .LT. 0.0) GO TO 40
      IF (A0 .GE. 0.0) GO TO 70

 70     SURFACE_SHAPE = ' CYLINDER'
      PFRAC = -100.0*BULK
      GO TO 130

 80     SURFACE_SHAPE = ' CONE'
      PFRAC = -SQRT(A0/A2)
      GO TO 130

 90     CONTINUE
      IF (A2 .LT. 0.0) GO TO 100
      IF (A2 .EQ. 0.0) GO TO 110
      IF (A2 .GT. 0.0) GO TO 120

 100    SURFACE_SHAPE = 'n ELLIPSOID'
      PFRAC = -0.5*(A1-SQRT(A3))/A2
      GO TO 130

 110    SURFACE_SHAPE = ' PARABOLOID'
      PFRAC = -A0/A1
      GO TO 130

 120    SURFACE_SHAPE = ' HYPERBOLOID'
      PFRAC = -0.5*(A1-SQRT(A3))/A2
      GO TO 130

 130    CONTINUE
      CALL USER_MESSAGE                                                        &
     & (                                                                       &
     & MSGL//'INFORM'//                                                        &
     & MSGL//'MATERIAL_38_INIT.006.00'//                                       &
     & MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                    &
     & MSGL//'Has Yield Surface Constants A0,A1,A2 That Produce'//             &
     & MSGL//'A Yield Surface That Is A'//SURFACE_SHAPE                        &
     & )

      IF (PFRAC .GT. 0.0) THEN
        CALL USER_MESSAGE                                                      &
     &   (                                                                     &
     &   MSGL//'FATAL'//                                                       &
     &   MSGL//'MATERIAL_38_INIT.007.00'//                                     &
     &   MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                  &
     &   MSGL//'Has Yield Surface Constants A0,A1,A2 That Produce'//           &
     &   MSGL//'A "Tensile Failure" In Compression.'                           &
     &   )
      ELSE
        MATERIAL(MatID)%PVAL(14) = PFRAC
        WRITE (MSGF,'(1PE12.4)') -PFRAC
        CALL USER_MESSAGE                                                      &
     &   (                                                                     &
     &   MSGL//'INFORM'//                                                      &
     &   MSGL//'MATERIAL_38_INIT.008.00'//                                     &
     &   MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                  &
     &   MSGL//'Has Yield Surface Constants A0,A1,A2 That Produce'//           &
     &   MSGL//'A Hydrostatic Tensile Failure At: '//MSGF                      &
     &   )
      ENDIF
!!
      RETURN
!!
      ENTRY MATERIAL_38_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!!
      Ioff = Isv - 1
      STATE_VARIABLES(Ioff+1) = 0.0   ! Initial bulk strain
      STATE_VARIABLES(Ioff+2) = 0.0   ! Initial minimum volume strain
      STATE_VARIABLES(Ioff+3) = 0.0   ! Initial maximum pressure
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_38                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 26-AUG-1995 16:57:45.00
!!
!! SOIL AND CRUSHABLE FOAM
!!
!! Module computes the current value of the stress STRESS(1:6), the
!! density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)

      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Generate the required rotation operator (and new stretching components
!! if Ipolard equals one or two).
!!
      Igenrot = 1
      Ipolard = CONTROL%POLARD
!!
      CALL SYMMETRIC_TENSOR_ROTATION                                           &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
!!
!! Rotate stress from the configuration at time n to the configuration
!! at time n+1 (Ipolard=0).
!!
      IF (Ipolard .EQ. 0) THEN
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
      ENDIF
!!
!! Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*DTnext) *                       &
     &          (                                                              &
     &           Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -        &
     &          (Dxx+Dyy+Dzz) * INTERNAL_ENERGY                                &
     &          )
!!
!! Update stress.
!!
      CALL MATERIAL_38_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,                                                        &
     &          STATE_VARIABLES,                                               &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
!! Internal energy increment from time n+1/2 to time n+1.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*DTnext) *                       &
     &          (                                                              &
     &           Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -        &
     &          (Dxx+Dyy+Dzz) * INTERNAL_ENERGY                                &
     &          )
!!
!! If the end-of-the-interval Rashid approximate polar decomposition is
!! being used (Ipolard=1,2), rotate stress from the configuration at
!! time n to the configuration at time n+1. Note: the rotation operator R
!! has been saved from the generation entry (Igenrot=1).
!!
      IF (Ipolard .NE. 0) THEN
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_38_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,                                                        &
     &          STATE_VARIABLES,                                               &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
!! Copyright (c) by KEY Associates; 26-AUG-1995 16:57:17.00
!!
!! SOIL AND CRUSHABLE FOAM MATERIAL MODEL:
!!
!!      Originally developed by Ray Krieg and Zelma Beisinger
!!             Sandia National Laboratories, 1973.
!!
!!          Ported and adapted to FMA-3D by Samuel Key
!!                 KEY Associates, August 1995.
!!
!! This module is an elastic-plastic isotropic constitutive relation for
!! the low pressure response of cellular concrete in particular. The
!! yield surface in principal stress space is a surface of revolution
!! centered about the hydrostat with the open end pointing into
!! compression. The open end is capped with a plane which is also normal
!! to the hydrostat. The deviatoric part is elastic-perfectly plastic so
!! the surface of revolution is stationary. The volumetric part has
!! variable strain hardening so the end plane moves outward during
!! volumetric yielding. Deviatoric strains produce no volume change. The
!! deviatoric yield function and associated flow rule are described as
!!
!!              PHI = J2 - (A0 + A1*P + A2*P*P)
!!
!! where J2 is the second invariant of the deviatoric stress, that is,
!!
!!              J2(S) = (S*S)/2 = (SijSij)/2
!!
!! MATERIAL PROPERTY INPUT:
!!
!!      SHEAR - Elastic shear loading modulus
!!      BULK  - Elastic bulk unloading modulus
!!      A0    - Constant used to describe yield function
!!      A1    - Constant used to describe yield function
!!      A2    - Constant used to describe yield function
!!      PVFTN - Pressure versus Volumetric strain FTN ID
!!
!! The tabulated function pairs are (volumetric strain, pressure),
!! where the volumetric strain (x-axis) is positive in compression
!! and the pressure (y-axis) is positive in compression. The first
!! pair must be the origin (0.0, 0.0). The slope of this curve is
!! the loading or compaction bulk modulus:
!!
!!      LOADING_BULK_MODULUS = dPRESSURE/dBULK_STRAIN
!!
!!      STATE_VARIABLES(1) = Current bulk strain (at time(n+1))
!!      STATE_VARIABLES(2) = Minimum volumetric strain seen
!!      STATE_VARIABLES(3) = Maximum pressure seen
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          STATE_VARIABLES(*),                                            &
     &          MINIMUM_STRAIN,MAXIMUM_PRESSURE
      INTEGER                                                                  &
     &          Table_ID
      LOGICAL                                                                  &
     &          TENSILE_FAILURE
      DATA                                                                     &
     &          THIRD /0.33333333333333/
!!
!! The vector storage sequence is STRESS(1:6) = (Sxx,Syy,Szz,Sxy,Sxz,Syz)
!! In this module, pressure is positive in compression. The tabulated function
!! has the compressive bulk strain as positive values (required by the table
!! look-up function; see MATERIAL_38_INIT for pre-processing).
!!
!! Retrieve material constants.
!!
      GG = MATERIAL(MatID)%PVAL(12)    ! 2 x shear modulus

      A0 = MATERIAL(MatID)%PVAL(6)     ! yield surface constant
      A1 = MATERIAL(MatID)%PVAL(7)     ! yield surface constant
      A2 = MATERIAL(MatID)%PVAL(8)     ! yield surface constant

      BULK  = MATERIAL(MatID)%PVAL(10) ! unloading bulk modulus
      PFRAC = MATERIAL(MatID)%PVAL(14) ! tensile failure pressure
!!
!! Decompose old stress into spherical and deviatoric parts.
!!
      POLD = -THIRD*(STRESS(1)+STRESS(2)+STRESS(3))
      TDO1 = STRESS(1) + POLD
      TDO2 = STRESS(2) + POLD
      TDO3 = STRESS(3) + POLD
      TDO4 = STRESS(4)
      TDO5 = STRESS(5)
      TDO6 = STRESS(6)
!!
!! PART I. VOLUMETRIC PROCESS.
!!
!! Compute current bulk strain.
!!
      BULK_STRAIN_RATE = (Dxx+Dyy+Dzz)
      BULK_STRAIN = STATE_VARIABLES(1) + DTnext * BULK_STRAIN_RATE
      STATE_VARIABLES(1) = BULK_STRAIN
!!
!! Retrieve minimum bulk strain and maximum pressure seen up to now.
!!
      MINIMUM_STRAIN   = STATE_VARIABLES(2)
      MAXIMUM_PRESSURE = STATE_VARIABLES(3)
!!
!! Check for volumetric unloading/loading.
!!
      IF (BULK_STRAIN .GE. MINIMUM_STRAIN) THEN
!!
!! Elastic volumetric unloading (plus tensile failure).
!!
        PRESSURE = MAXIMUM_PRESSURE +                                          &
     &    BULK * (MINIMUM_STRAIN - BULK_STRAIN)

        TENSILE_FAILURE = (PRESSURE .LT. PFRAC)
        IF (TENSILE_FAILURE) THEN
          STRESS(1) = -PFRAC
          STRESS(2) = -PFRAC
          STRESS(3) = -PFRAC
          STRESS(4) = 0.0
          STRESS(5) = 0.0
          STRESS(6) = 0.0
        ENDIF
!!
!! Plastic volumetric loading (further compaction).
!!
      ELSE
        STATE_VARIABLES(2) = BULK_STRAIN
        Table_ID = NINT (MATERIAL(MatID)%PVAL(11))
        PRESSURE = TABLE_LOOK_UP (Table_ID,-BULK_STRAIN)
        STATE_VARIABLES(3) = PRESSURE
        TENSILE_FAILURE = .FALSE.
      ENDIF
!!
!! PART II. DEVIATORIC PROCESS.
!!
!! Check for tensile failure (in which event no deviatoric processing
!! is required).
!!
      IF (.NOT.TENSILE_FAILURE) THEN
!!
!! Elastic deviatoric stress increment.
!!
        TDRE1 = (GG*DTnext)*(Dxx - (THIRD*BULK_STRAIN_RATE))
        TDRE2 = (GG*DTnext)*(Dyy - (THIRD*BULK_STRAIN_RATE))
        TDRE3 = (GG*DTnext)*(Dzz - (THIRD*BULK_STRAIN_RATE))
        TDRE4 = (GG*DTnext)*(Dxy)
        TDRE5 = (GG*DTnext)*(Dxz)
        TDRE6 = (GG*DTnext)*(Dyz)
!!
!! Compute deviatoric elastic stress components.
!!
        TDE1 = TDO1 + TDRE1
        TDE2 = TDO2 + TDRE2
        TDE3 = TDO3 + TDRE3
        TDE4 = TDO4 + TDRE4
        TDE5 = TDO5 + TDRE5
        TDE6 = TDO6 + TDRE6
!!
!! Evaluate square of deviatoric yield surface radius (using current pressure).
!!
        YIELD2 = A0 + PRESSURE*(A1 + A2*PRESSURE)
!!
!! Evaluate deviatoric elastic stress components for yielding.
!!
        PHIDE = 0.5*(TDE1*TDE1+TDE2*TDE2+TDE3*TDE3                             &
     &          + 2.0*(TDE4*TDE4+TDE5*TDE5+TDE6*TDE6))-YIELD2

        IF (PHIDE .LE. 0.0) THEN
!!
!! Elastic: define deviatoric stress components at the end of the step.
!!
          TD1 = TDE1
          TD2 = TDE2
          TD3 = TDE3
          TD4 = TDE4
          TD5 = TDE5
          TD6 = TDE6
        ELSE
!!
!! Plastic: but prior to doing a radial return calculation, check
!! to see if the stress state at the beginning of the step needs to
!! be returned to the yield surface due to change in the pressure.
!!
          PHIO = 0.5*(TDO1*TDO1+TDO2*TDO2+TDO3*TDO3                            &
     &           + 2.0*(TDO4*TDO4+TDO5*TDO5+TDO6*TDO6))-YIELD2
          IF (PHIO .GT. 0.0) THEN
            RATIO = SQRT (0.9999*YIELD2/(PHIO+YIELD2))
            TDO1 = RATIO*TDO1
            TDO2 = RATIO*TDO2
            TDO3 = RATIO*TDO3
            TDO4 = RATIO*TDO4
            TDO5 = RATIO*TDO5
            TDO6 = RATIO*TDO6
!!
!! and redefine deviatoric elastic stress components,
!!
            TDE1 = TDO1+TDRE1
            TDE2 = TDO2+TDRE2
            TDE3 = TDO3+TDRE3
            TDE4 = TDO4+TDRE4
            TDE5 = TDO5+TDRE5
            TDE6 = TDO6+TDRE6
!!
!! and check once again for yielding at the end of the step.
!!
            PHIDE = 0.5*(TDE1*TDE1+TDE2*TDE2+TDE3*TDE3                         &
     &              + 2.0*(TDE4*TDE4+TDE5*TDE5+TDE6*TDE6))-YIELD2
          ENDIF

          IF (PHIDE .LE. 0.0) THEN
!!
!! Elastic: define deviatoric stress components at the end of the step.
!!
            TD1 = TDE1
            TD2 = TDE2
            TD3 = TDE3
            TD4 = TDE4
            TD5 = TDE5
            TD6 = TDE6
!!
!! Location of the deviatoric trial stress wrt the deviatoric yield surface.
!! (from an earlier informative print; maybe again in the future.)
!!
!!            YIELD_RADII = -1.0
          ELSE
!!            YIELD_RADII = SQRT (1.0+PHIDE/YIELD2)
!!
!! Plastic: relocate points outside of yield surface back on yield surface.
!! (Now, we really do have a plastic flow to calculate!)
!!
            SIZE1 =  TDE1*TDRE1+TDE2*TDRE2+TDE3*TDRE3                          &
     &              + (TDE4*TDRE4+TDE5*TDRE5+TDE6*TDRE6)                       &
     &              + (TDE4*TDRE4+TDE5*TDRE5+TDE6*TDRE6)
            SIZE2 =  TDRE1*TDRE1+TDRE2*TDRE2+TDRE3*TDRE3                       &
     &              + (TDRE4*TDRE4+TDRE5*TDRE5+TDRE6*TDRE6)                    &
     &              + (TDRE4*TDRE4+TDRE5*TDRE5+TDRE6*TDRE6)
!!
!! Elastic fraction of step.
!!
            R = 1.0
            EFRAC = (SIZE1-SQRT(SIZE1*SIZE1-PHIDE*(SIZE2+SIZE2)))/SIZE2
            IF (EFRAC .EQ. 0.0) GO TO 210
            COSPHI = (SIZE1-EFRAC*SIZE2)/SQRT(SIZE2*(YIELD2+YIELD2))
            C1 = EFRAC*SQRT(SIZE2/(YIELD2+YIELD2))
            C2 = EXP(-C1)
            IF ((C1-7.0) .LE. 0.0) GO TO 200
            RMA = 0.0
            RMB = SQRT((YIELD2+YIELD2)/SIZE2)
            GO TO 230
!!
!! Reduction factor on step size.
!!
 200          R = (1.0-C2*C2+COSPHI*((1.0-C2)*(1.0-C2)))/                      &
     &            ((C1*C2)+(C1*C2))
 210          SIZE3 = (YIELD2+PHIDE)+(YIELD2+PHIDE)+                           &
     &            ((SIZE1+SIZE1)+(R-EFRAC)*SIZE2)*(R-EFRAC)
            IF (ABS(SIZE3) .LE. 1.0D-30*(YIELD2+YIELD2)) THEN
              WRITE (MSG1,'(I8)') MATERIAL(MatID)%MatID
              CALL USER_MESSAGE                                                &
     &  (                                                                      &
     &  MSGL//'FATAL'//                                                        &
     &  MSGL//'MATERIAL_38_INTEGRATION.001.00'//                               &
     &  MSGL//'Material (MATERIAL) Input Record ID:'//MSG1//                   &
     &  MSGL//'Radial Return Algorithm Is Trying To Divide By 0.0'             &
     &  )
            ENDIF
            RMA = SQRT ((YIELD2+YIELD2)/SIZE3)
            RMB = RMA*(R-EFRAC)
!!
!! Radial movement back to the yield surface.
!!
 230          CONTINUE
            TD1 = RMA*TDE1 + RMB*TDRE1
            TD2 = RMA*TDE2 + RMB*TDRE2
            TD3 = RMA*TDE3 + RMB*TDRE3
            TD4 = RMA*TDE4 + RMB*TDRE4
            TD5 = RMA*TDE5 + RMB*TDRE5
            TD6 = RMA*TDE6 + RMB*TDRE6
          ENDIF
        ENDIF

        STRESS(1) = TD1 - PRESSURE
        STRESS(2) = TD2 - PRESSURE
        STRESS(3) = TD3 - PRESSURE
        STRESS(4) = TD4
        STRESS(5) = TD5
        STRESS(6) = TD6
      ENDIF
!!
!! Sound speeds squared * RHO(t)
!!
      G43 = MATERIAL(MatID)%PVAL(13)
      SHEAR = MATERIAL(MatID)%PVAL(9)

      SOUND_SPEED%RCL2 = BULK + G43
      SOUND_SPEED%RCS2 = SHEAR
!!
      RETURN
      END
