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
      SUBROUTINE MATERIAL_D6_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 8-AUG-1993 14:41:06.36
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_D6_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D7_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 8-AUG-1993 14:44:24.11
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_D7_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D8_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 8-AUG-1993 14:44:29.18
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_D8_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D9_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 8-AUG-1993 14:44:29.18
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_D9_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D6                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:14:25
!!
!! LINEAR FORCE-VELOCITY VISCOUS BEHAVIOR
!!
!! This subroutine computes the current value of the force FORCE, and
!! the longitudinal stiffness as a function of the velocity Drr.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Save old force to compute incremental stiffness.
!!
      Fold = FORCE
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Update FORCE.
!!
      FORCE = MATERIAL(MatID)%PVAL(6) * Drr
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Current damper stiffness.
!!
      Finc = FORCE - Fold
      IF (ABS(DTnext*Drr) .GT. 1.0E-25) THEN
        SOUND_SPEED%RCL2 = ABS(Finc / (DTnext*Drr))
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = 1.0
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D7                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:14:25
!!
!! LINEAR BI_DIRECTIONAL FORCE-VELOCITY VISCOUS BEHAVIOR
!!
!! This subroutine computes the current value of the force FORCE, and
!! the longitudinal stiffness as a function of the velocity Drr.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)

!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Save old force to compute incremental stiffness.
!!
      Fold = FORCE
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Update FORCE%
!!
      IF (Drr .GE. 0.0) THEN
        FORCE = MATERIAL(MatID)%PVAL(6) * Drr
      ELSE
        FORCE = MATERIAL(MatID)%PVAL(7) * Drr
      ENDIF
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Current damper stiffness.
!!
      Finc = FORCE - Fold
      IF (ABS(DTnext*Drr) .GT. 1.0E-25) THEN
        SOUND_SPEED%RCL2 = ABS(Finc / (DTnext*Drr))
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = 1.0
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D8                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:14:25
!!
!! PIECEWISE LINEAR FORCE-VELOCITY VISCOUS BEHAVIOR
!!
!! This subroutine computes the current value of the force FORCE, and
!! the longitudinal stiffness as a function of the velocity Drr.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          FtnID           ! Used to convert R*4 to I*4
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          TABLE_LOOK_UP,                                                 &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Save old force to compute incremental stiffness.
!!
      Fold = FORCE
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Update FORCE by using a table look-up.
!!
      FtnID = NINT (MATERIAL(MatID)%PVAL(6))
      FORCE = TABLE_LOOK_UP (FtnID,Drr)
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Current spring stiffness
!!
      Finc = FORCE - Fold
      IF (ABS(DTnext*Drr) .GT. 1.0E-25) THEN
        SOUND_SPEED%RCL2 = ABS(Finc / (DTnext*Drr))
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = 1.0
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_D9                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 6-FEB-1992 20:04:28
!!
!! USER DEFINED DAMPER CHARACTERISTICS
!!
!! This subroutine is required to compute the current value of the damper
!! force FORCE, and the longitudinal stiffness as a function of the velocity
!! Drr. If Drr is positive, the damper is extending. If Drr is negative, the
!! damper is shortening.
!!
!! FORCE
!! Input:  Contains the damping force from the last time step.
!! Output: Contains the damping force for this time step.
!!
!! STATE_VARIABLES(1:15)
!! Input:  Contains the values of the damper's history variables, if any,
!!         from the last time step.
!! Output: Contains the values of the damper's history variables, if any,
!!         for this time step.
!!
!!         STATE_VARIABLES( 1) = user's chioce
!!         STATE_VARIABLES( 2) = user's chioce
!!         STATE_VARIABLES( 3) = user's chioce
!!         : : :
!!         STATE_VARIABLES(15) = user's chioce
!!
!! INTERNAL_ENERGY
!! Input:  Contains the internal energy disapated up to this time step.
!! Output: Contains the internal energy disapated after this time step.
!!
!! DTnext
!! Input:  The time step over which the damping is to be calculated.
!! Output: None. DO NOT MODIFY!
!!
!! MatID
!! Input:  Contains the user's damping constants. MATERIAL(MatID) is a
!!         "structure" The constants are accessed as MATERIAL(MatID)%PVAL(6)
!!         through MATERIAL(MatID)%PVAL(22) The sample implmentation shows
!!         how they are used in a Fortran statement.
!! Output: None.
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          FtnID
      REAL(KIND(0D0))                                                          &
     &          FORCE,                                                         &
     &          DTnext,                                                        &
     &          TABLE_LOOK_UP,                                                 &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(15)
      LOGICAL, SAVE:: FIRST
      DATA                                                                     &
     &          FIRST /.TRUE./
!!
!! This common block contains the unit vector R that defines the direction
!! in which the damper is acting, that is, the orientation of the damper
!! element. It also contains the relative velocity of the two nodal points
!! that define the element.
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! OPTIONAL: The following is an example of how to use the standard mess-
!! age module for placing messages in the printer output. Two other mess-
!! ages also exist: (1) fatal messages causing termination of the execution
!! and indicated by 'FATAL' and (2) warning messages NOT causing termin-
!! ation and indicated by 'WARN'. (The message is contructed as a single
!! character string. The module USER_MESSAGE breaks the string into lines
!! by looking for the message line character MSGL. The string may contain
!! any number of MSGL characters greater than 1.
!!
      IF (FIRST) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'MATERIAL_D9.000.00'//                                   &
     &          MSGL//'First entry into MATERIAL_D9.'                          &
     &          )
        FIRST = .FALSE.
      ENDIF
!!
!! REQUIRED: Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! OPTIONAL: Save old force to compute incremental stiffness.
!!
      Fold = FORCE
!!
!! EXAMPLE USE OF DAMPING CONSTANTS: Select appropriate tabulated function
!! for interogation.
!!
      IF (Drr .GE. 0.0) THEN
        FtnID = NINT (MATERIAL(MatID)%PVAL(6))
      ELSE
        FtnID = NINT (MATERIAL(MatID)%PVAL(7))
      ENDIF
!!
!! EXAMPLE USE OF TABLE LOOK-UP FUNCTION: Update FORCE by using a table
!! look-up.
!!
      FORCE = TABLE_LOOK_UP (FtnID,Drr)
!!
!! REQUIRED: Internal energy increment from time n+1/2 to time n.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! REQUIRED: Current effective stiffness for critical time step calculation.
!!
      Finc = FORCE - Fold
      IF (ABS(DTnext*Drr) .GT. 1.0E-25) THEN
        SOUND_SPEED%RCL2 = ABS(Finc / (DTnext*Drr))
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = 1.0
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
