!! MATERIAL(*)%PVAL(1:22) Usage
!!
!! 1:Density 2:Bulk_Ln 3:Bulk_Qd 4:HG_Viscosity 5:HG_Stiffness
!!
!! PVAL 1/2/3 6/7/8 10/20/30 11/21/31 22  25/35/45 32   33   36  17/27 38
!!                  40/50    41/51                                 37
!!  6:  K1    D1    E        E        EALn  E_aa  E          Ko        A0
!!  7:        D2    Nu       Nu       EAQd  E_bb  Nu         Kinf      A1
!!  8:              Lambda   Lambda   EBLn  E_cc  Lam        Kdec      A2
!!  9:              G        G        EBQd  Nu_ba G     G    G0   Gftn G
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
      SUBROUTINE MATERIAL_22_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:56
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(3)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! If the bulk modulus has been supplied as zero, compute a bulk modulus
!! based on the average value of the two fiber linear stiffness coefficients.
!! Poisson's ratio is taken to be 0.3; thus, 3(1-2nu) = 1.2
!!
      Fbulk = MATERIAL(MatID)%PVAL(11)
      IF (Fbulk .LE. 0.0) THEN
        EALinear = MATERIAL(MatID)%PVAL(6)
        EBLinear = MATERIAL(MatID)%PVAL(8)
        Fbulk = 0.5 * (EALinear + EBLinear) / 1.2
        MATERIAL(MatID)%PVAL(11) = Fbulk
      ENDIF
!!
      RETURN
!!
      ENTRY MATERIAL_22_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! Initialize stress point w/ initial fiber orientation; global coordinates.
!!
      Ioff = Isv - 1
      STATE_VARIABLES(Ioff+1) = SECTION_2D(SecID)%Ax
      STATE_VARIABLES(Ioff+2) = SECTION_2D(SecID)%Ay
      STATE_VARIABLES(Ioff+3) = SECTION_2D(SecID)%Az
      STATE_VARIABLES(Ioff+4) = SECTION_2D(SecID)%Bx
      STATE_VARIABLES(Ioff+5) = SECTION_2D(SecID)%By
      STATE_VARIABLES(Ioff+6) = SECTION_2D(SecID)%Bz
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_22                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-SEP-1992 09:17:32.13
!!
!! TENSION-ONLY, PLANE STRESS, NONLINEAR ELASTIC FABRIC W/ FIBER FRICTION
!!
!! This subroutine computes the current value of the stress STRESS(1:3), the
!! longitudinal sound speed CL, and the shear sound speed CS in the current
!! state, as a function of the stretching Drr,Dss,Drs the spin Wrs, the past
!! history data STATE_VARIABLES(1:7), and the old stress state.
!!
!!      STATE_VARIABLES(1) = A fiber-vector local r-component
!!      STATE_VARIABLES(2) = A fiber-vector local s-component
!!      STATE_VARIABLES(3) = B fiber-vector local r-component
!!      STATE_VARIABLES(4) = B fiber-vector local s-component
!!      STATE_VARIABLES(5) = A fiber strain
!!      STATE_VARIABLES(6) = B fiber strain
!!      STATE_VARIABLES(7) = AB fiber shear strain
!!      STATE_VARIABLES(8) = Areal strain, Ln(Area(t)/Area(0))
!!
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!      *  NOTE: The plane-stress direction is along the t-axis.  *
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(3),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching, spin
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Gr,Gs,                  & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
!! Since this is a model which tracks fiber directions and uses total fiber
!! strain as a state variable, rotation of the stress state from time n to
!! time n+1 is not required. The fiber orientation must be updated instead.
!!
      IF (PARAMVALUE%FIBROT .EQ. 0.0) THEN
!!
!! This is a rotation of each fiber. The angle between the two fabric fiber
!! directions changes with time. The rotation increments R1 and R2 reflect the
!! rotation of the individual fibers with respect to the local element r,s,t-
!! coordinates. This rotation treatment tracks large rotations of one fiber
!! relative to the other fiber that can lead to the two fibers lining up with
!! each other and producing unwanted physical behavior, that is, an extreme
!! orthotropic element stiffness.
!!
        Vrs = (Drs + Wrs)
        Vsr = (Drs - Wrs)
!!
        Ar  = STATE_VARIABLES(1)
        As  = STATE_VARIABLES(2)
        Vnt = (Vsr*Ar + Dss*As)*Ar - (Drr*Ar + Vrs*As)*As
        R1  = -(DTnext * (Vnt-Vsr))
!!
        Ar  = STATE_VARIABLES(3)
        As  = STATE_VARIABLES(4)
        Vnt = (Vsr*Ar + Dss*As)*Ar - (Drr*Ar + Vrs*As)*As
        R2  = -(DTnext * (Vnt-Vsr))
!!
      ELSE
!!
!! This is a "rigid body rotation." The angle between the two fabric fiber
!! directions remains fixed. The rotation increments R1 and R2 reflect the
!! average rotation of the fabric with respect to the local element r,s,t-
!! coordinates. This rotation treatment ignores large rotations of one fiber
!! relative to the other fiber that can lead to the two fibers lining up with
!! each other and producing unwanted physical behavior, that is, an extreme
!! orthotropic element stiffness.
!!
        Vsr = Drs-Wrs
!!
        R1 = (DTnext * (Wrs+Vsr))
        R2 = (DTnext * (Wrs+Vsr))
!!
      ENDIF
!!
!! Based on the incremental rotation values R1 and R2, rotate each fiber.
!! Note: the unit vectors defining the fiber directions are re-normalized.
!!
      Q1 = STATE_VARIABLES(1) + R1*STATE_VARIABLES(2)
      Q2 = STATE_VARIABLES(2) - R1*STATE_VARIABLES(1)
      Qim = 1.0 / SQRT (Q1*Q1 + Q2*Q2)
      STATE_VARIABLES(1) = Q1 * Qim
      STATE_VARIABLES(2) = Q2 * Qim
!!
      Q1 = STATE_VARIABLES(3) + R2*STATE_VARIABLES(4)
      Q2 = STATE_VARIABLES(4) - R2*STATE_VARIABLES(3)
      Qim = 1.0 / SQRT (Q1*Q1 + Q2*Q2)
      STATE_VARIABLES(3) = Q1 * Qim
      STATE_VARIABLES(4) = Q2 * Qim
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5 * DTnext) *                     &
     &          (Drr*STRESS(1) + Dss*STRESS(2) + (Drs+Drs)*STRESS(3))
!!
!! Update stress.
!!
      CALL MATERIAL_22_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,                                                        &
     &          STATE_VARIABLES(1),                                            &
     &          STATE_VARIABLES(3),                                            &
     &          STATE_VARIABLES(5),                                            &
     &          STATE_VARIABLES(6),                                            &
     &          STATE_VARIABLES(7),                                            &
     &          STATE_VARIABLES(8),                                            &
     &          DTnext,                                                        &
     &          Drr,Dss,Drs,                                                   &
     &          MatID                                                          &
     &          )
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5 * DTnext) *                     &
     &          (Drr*STRESS(1) + Dss*STRESS(2) + (Drs+Drs)*STRESS(3))
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_22_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,A,B,Strain_A,Strain_B,Strain_C,Akk,                     &
     &          DTnext,Drr,Dss,Drs,MatID                                       &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-SEP-1992 09:17:49.90
!!
!! TENSION-ONLY, PLANE STRESS, NONLINEAR ELASTIC FABRIC W/ FIBER FRICTION
!!
!! This is a model for bi-directional, woven fabric. Each fiber direction is
!! tracked with a unit vector parallel to the fiber; the strain is a state
!! variable for the fiber to simplify implementation of the tension-only
!! limitation. The shear carrying capacity is either based on the friction
!! between the woven fibers as they shear, or is based on the out-of-plane
!! curvature changes as the fibers distort to accomodate a new geometry
!! during shear.
!!
!! If the friction model for shear is used, maximum shear stress is given
!! by a friction coefficient times the mean value of the direct stresses:
!!  Shear = COF * SQRT (MAX (Sigma_A, 0.0) * MAX (Sigma_B, 0.0)).
!!
!! If the elastic model for shear is used, shear stress is given by a shear
!! modulus G(P) linearly dependent on the mean value of the direct stresses:
!!  G(P) = dG/dP x SQRT (MAX (Sigma_A, 0.0) * MAX (Sigma_B, 0.0)),
!! and
!!  Shear = 2G(P) x Shear_C.
!!
!!   MATERIAL%K1a   = 1st fiber Young's modulus, linear coefficient, E1A
!!   MATERIAL%K2a   = 1st fiber Young's modulus, quadratic coefficient, E2A
!!   MATERIAL%K1b   = 2nd fiber Young's modulus, linear coefficient, E1B
!!   MATERIAL%K2b   = 2nd fiber Young's modulus, quadratic coefficient, E2B
!!   MATERIAL%CoF   = Shear coefficient of friction, COF
!!   MATERIAL%Fblk  = Bulk modulus in compression, Blk
!!   MATERIAL%Fcomp = Fully compacted areal strain, Afc
!!   MATERIAL%Smodel= Shear model, (0,1=friction,elastic)
!!   MATERIAL%G1c   = Shear modulus coefficient, dG/dP, P=SQRT(Sigma_A*Sigma_B)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          STRESS(3),A(2),B(2)
!!
!! Retrieve linear and quadratic fiber coefficients, fiber weave cofficient of
!! friction, compressive bulk modulus, and the fully compacted areal strain.
!!
      EA1 = MATERIAL(MatID)%PVAL(6)
      EA2 = MATERIAL(MatID)%PVAL(7)
      EB1 = MATERIAL(MatID)%PVAL(8)
      EB2 = MATERIAL(MatID)%PVAL(9)
      COF = MATERIAL(MatID)%PVAL(10)
      Blk = MATERIAL(MatID)%PVAL(11)
      Afc = MATERIAL(MatID)%PVAL(12)
      SMD = MATERIAL(MatID)%PVAL(13)
      G1C = MATERIAL(MatID)%PVAL(14)
!!
      A11 = A(1)*A(1)
      A22 = A(2)*A(2)
      A12 = A(1)*A(2)
      B11 = B(1)*B(1)
      B22 = B(2)*B(2)
      B12 = B(1)*B(2)
!!
!! Integrate strain rate aligned with fibers. The resulting strain is
!! logarithmic, Strain = LOGe(L(t)/L(0)).
!!
      Strain_A = Strain_A + DTnext * (Drr*A11+Dss*A22+(Drs+Drs)*A12)
      Strain_B = Strain_B + DTnext * (Drr*B11+Dss*B22+(Drs+Drs)*B12)
!!
!! Compute fiber non-linear elastic stresses.
!!
      IF (Strain_A .GT. 0.0) THEN
        Sigma_A = (EA1 + EA2 * Strain_A) * Strain_A
      ELSE
        Sigma_A = 0.0
      ENDIF

      IF (Strain_B .GT. 0.0) THEN
        Sigma_B = (EB1 + EB2 * Strain_B) * Strain_B
      ELSE
        Sigma_B = 0.0
      ENDIF
!!
!! Construct stress due to the direct (axial) stress in each fiber.
!!
      STRESS(1) = Sigma_A*A11 + Sigma_B*B11
      STRESS(2) = Sigma_A*A22 + Sigma_B*B22
      STRESS(3) = Sigma_A*A12 + Sigma_B*B12
!!
!! Compute shear stress due to friction/elasticity in the weave.
!!
      PAB = Sigma_A * Sigma_B
      IF (PAB .GT. 0.0) THEN
!!
!! Compute the unit vector C bisecting the vectors A and B. Construct the
!! unit vector D orthogonal to C.
!!
        C1 = A(1) + B(1)
        C2 = A(2) + B(2)
        Cim = 1.0 / (C1*C1 + C2*C2)
        C11 = C1*C1 * Cim
        C22 = C2*C2 * Cim
        C12 = C1*C2 * Cim
        D11 = +C22
        D22 = +C11
        D12 = -C12
!!
!! Compute stretchings aligned with the vectors C and D and average them to get
!! the shearing Dcc. (A positive value for Dcc indicates the directions A and B
!! are coming together; a negative value for Dcc indicates the directions A and
!! B are spreading apart.)
!!
        Dcc = 0.5 * (Drr*(C11-D11)+Dss*(C22-D22)+(Drs+Drs)*(C12-D12))
!!
!! Select shear model. The shear models assume a simple "orthogonal" weave,
!! the fibers alternately pass over and then under each other.
!!
        IF (SMD .EQ. 0.0) THEN
!!
!! Friction Model For Shear (SMD equal to zero).
!! The shear stress arises from the fibers rubbing together in the weave
!! as the fiber directions A and B change relative to one another. The
!! greater the tension in both fiber directions, the greater is the
!! pressure that the fibers exert on each other as they cross in the weave
!! and the greater is the shear stress needed to cause shearing in the weave.
!! If either fiber direction is in compression, only minimal or zero shear
!! stress due to friction will occur; the model here gives zero. This model
!! is, in effect, a "rigid-plastic" formulation for shear.
!!
          Shear = SIGN (COF * SQRT (PAB), Dcc)
!!
        ELSE
!!
!! Stiffness Model For Shear (SMD not equal to zero).
!! The shear stress arises from the increase in bearing area between the
!! fiber contacts in the weave and increased out-of-plane flexing of the
!! fibers as the fiber directions A and B change relative to one another.
!! The greater the tension in both fiber directions, the the greater is
!! the pressure that the fibers exert on each other as they cross in the
!! weave and the greater is the shear stress needed to cause shearing in
!! the weave. If either fiber direction is in compression, only minimal or
!! zero shear stress due to fiber distortion will occur; the model here
!! gives zero.
!!
          Strain_C = Strain_C + DTnext * Dcc
          Gmod = G1C * SQRT (PAB)
          Shear = (Gmod + Gmod) * Strain_C
!!
        ENDIF
!!
!! Acummulate stress due to shearing in the weave.
!!
        STRESS(1) = STRESS(1) + (Shear * (C11 - D11))
        STRESS(2) = STRESS(2) + (Shear * (C22 - D22))
        STRESS(3) = STRESS(3) + (Shear * (C12 - D12))
      ENDIF
!!
!! Update bulk/volume strain. For a membrane element this is the areal strain
!! given by LOGe(area(t)/area(0)). For example, if the current area of an
!! element is one half of its original area, Akk = LOGe(0.5) = -0.693
!!
      Akk = Akk + DTnext * (Drr + Dss)
!!
!! In the event the areal strain of the element drops below the "fully
!! compacted" value Afc, compute a pressure (positive in tension).
!! Note: this is as much a numerical procedure to prevent the "collapse"
!! of tension-only fabric elements as it is a representation of the physical
!! compaction of fabric materials.
!!
      IF (Akk .LE. Afc) THEN
        Pressure = Blk * (Akk - Afc)
!!
!! Acummulate pressure due to compaction of the frabric.
!!
        STRESS(1) = STRESS(1) + Pressure
        STRESS(2) = STRESS(2) + Pressure
      ENDIF
!!
!! Sound speeds squared * rho(t)
!!
      Q1 = MAX (EA1, (EA1 + 2.0*EA2*Strain_A))
      Q2 = MAX (EB1, (EB1 + 2.0*EB2*Strain_B))
      SOUND_SPEED%RCL2 = MAX (Q1, Q2)
      SOUND_SPEED%RCS2 = SOUND_SPEED%RCL2 * 0.385
!!
      RETURN
      END
