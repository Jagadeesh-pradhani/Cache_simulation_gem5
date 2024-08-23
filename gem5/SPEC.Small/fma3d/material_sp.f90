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
      SUBROUTINE MATERIAL_S1_INIT (MatID)
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
      ENTRY MATERIAL_S1_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S2_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 8-AUG-1993 14:44:19.84
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
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H. (Note that MATERIAL%Ymod = MATERIAL%K1)
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_S2_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S3_INIT (MatID)
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
      ENTRY MATERIAL_S3_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S5_INIT (MatID)
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
      ENTRY MATERIAL_S5_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S1                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:14:25
!!
!! LINEAR FORCE-DISPLACEMENT BEHAVIOR
!!
!! This subroutine computes the current value of the force, and the
!! longitudinal stiffness as a function of the displacement rate Drr
!! and the old force value.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: FORCE
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Update FORCE.
!!
      FORCE = FORCE + MATERIAL(MatID)%PVAL(6) * (DTnext*Drr)
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! Current spring stiffness
!!
      SOUND_SPEED%RCL2 = MATERIAL(MatID)%PVAL(6)
      SOUND_SPEED%RCS2 = 0.0
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S2                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:28:00
!!
!! ELASTIC-PLASTIC FORCE-DISPLACEMENT BEHAVIOR.
!!
!! This subroutine computes the current value of the force FORCE, the
!! longitudinal stiffness as a function of the displacement rate Drr,
!! the old force value, and the past history data STATE_VARIABLES(1:2).
!!
!!      STATE_VARIABLES(1) = Yield surface center, Alpha.
!!      STATE_VARIABLES(2) = Effective plastic strain, Efps
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: FORCE
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*DTnext*Drr)*FORCE
!!
!! Update force.
!!
      CALL MATERIAL_11_INTEGRATION                                             &
     &          (                                                              &
     &          FORCE,                                                         &
     &          STATE_VARIABLES(1),                                            &
     &          STATE_VARIABLES(2),                                            &
     &          DTnext,                                                        &
     &          Drr,                                                           &
     &          MatID                                                          &
     &          )
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*DTnext*Drr)*FORCE
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S3                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:14:25
!!
!! PIECEWISE LINEAR FORCE-DISPLACEMENT ELASTIC BEHAVIOR
!!
!! This subroutine computes the current value of the force FORCE, the
!! longitudinal stiffness as a function of the displacement rate Drr
!! and the old force value.
!!
!!      STATE_VARIABLES(1) = Spring displacement
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: FORCE
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID
!!
!! Local variables.
      INTEGER         :: FtnID
      REAL(KIND(0D0)) :: TABLE_LOOK_UP
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Compute displacement increment.
!!
      Einc = DTnext * Drr
!!
!! Update spring displacement.
!!
      STATE_VARIABLES(1) = STATE_VARIABLES(1) + Einc
!!
!! Save old force to compute incremental stiffness.
!!
      Fold = FORCE
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*Einc)*FORCE
!!
!! Update FORCE by using a table look-up.
!!
      FtnID = NINT (MATERIAL(MatID)%PVAL(6))
      FORCE = TABLE_LOOK_UP (FtnID,STATE_VARIABLES(1))
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*Einc)*FORCE
!!
!! Current spring stiffness.
!!
      Finc = FORCE - Fold
      IF (ABS (Einc) .GT. 1.0D-25) THEN
        SOUND_SPEED%RCL2 = ABS (Finc / Einc)
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = ONE
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_S5                                                   &
     &          (                                                              &
     &          FORCE,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID             &
     &          )
!!
!! Copyright (c) by KEY Associates; 6-FEB-1992 21:28:39
!!
!! USER DEFINED SPRING CHARACTERISTICS
!!
!! This subroutine is required to compute the current value of the spring
!! force FORCE, and the longitudinal stiffness as a function of the axial
!! relative velocity Drr. If Drr is positive, the spring is extending. If
!! Drr is negative, the spring is shortening.
!!
!! FORCE
!! Input:  Contains the spring force from the last time step.
!! Output: Contains the spring force for this time step.
!!
!! STATE_VARIABLES(1:15)
!! Input:  Contains the values of the spring's history variables, if any,
!!         from the last time step.
!! Output: Contains the values of the spring's history variables, if any,
!!         for this time step.
!!
!!         (In this example user spring-stiffness module, the first state
!!         variable location is used to accumulate axial relative displace-
!!         ment increments thereby having a record of the springs total
!!         deflection.)
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
!! Input:  The time step for which the spring force is to be calculated.
!! Output: None. Do not modify.
!!
!! MatID
!! Input:  Contains the user's spring constants. MATERIAL(MatID) is a
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
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: FORCE
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(15)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID
!!
!! Local variables.
      INTEGER         :: FtnID
      REAL(KIND(0D0)) :: TABLE_LOOK_UP
!!
      LOGICAL :: FIRST = .TRUE.
!!
!! This common block contains the unit vector R that defines the direction
!! in which the spring is acting, that is, the orientation of the spring
!! element. It also contains the axial relative velocity Drr between the
!! two nodal points that define the element.
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
     &          MSGL//'MATERIAL_S5.000.00'//                                   &
     &          MSGL//'First entry into MATERIAL_S5.'                          &
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
!! EXAMPLE USE OF STATE VARIABLE ARRAY: Update the spring displacement by
!! integrating the axial relative velocity.
!!
      STATE_VARIABLES(1) = STATE_VARIABLES(1) + DTnext*Drr
      Axial_Displacement = STATE_VARIABLES(1)
!!
!! EXAMPLE USE OF SPRING CONSTANTS: Select appropriate tabulated function
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
      FORCE = TABLE_LOOK_UP (FtnID,Axial_Displacement)
!!
!! REQUIRED: Internal energy increment from time n+1/2 to time n.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*FORCE
!!
!! REQUIRED: Current spring stiffness for critical time step calculation.
!!
      Finc = FORCE - Fold
      IF (ABS(DTnext*Drr) .GT. 1.0D-25) THEN
        SOUND_SPEED%RCL2 = ABS(Finc / (DTnext*Drr))
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = ONE
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
