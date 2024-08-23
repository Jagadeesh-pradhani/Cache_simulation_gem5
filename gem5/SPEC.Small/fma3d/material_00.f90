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
      SUBROUTINE INITIALIZE_MATERIALS
!!
!! Copyright (c) by KEY Associates; 19-FEB-1991 22:17:07
!!
!! Purpose: Initialize material models and the state variable storage. In
!! particular, initial stresses are modified if required to assure they
!! are admissable stress states for the constitutive model specified. For
!! example, an initial stress state which exceeds the yield condition in
!! an elastic-plastic material model will be modified accordingly.
!!
      USE shared_common_data
      USE material_
      USE hexah_
      USE penta_
      USE tetra_
      USE lsold_
      USE membt_
      USE membq_
      USE truss_
      USE platt_
      USE platq_
      USE beam_
      USE spring_
      USE damper_
      USE spring_bc_
      USE damper_bc_
      USE section_2d_
      USE section_1d_
      USE layering_
      USE tabulated_function_
      USE stress_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SecID
!!
!! Initialize material constants.
!!
      DO M = 1,NUMMT
      MatID = M
      MType = MATERIAL(M)%Type
        IF (MType .EQ.  1) CALL MATERIAL_S1_INIT(MatID)
        IF (MType .EQ.  2) CALL MATERIAL_S2_INIT(MatID)
        IF (MType .EQ.  3) CALL MATERIAL_S3_INIT(MatID)
        IF (MType .EQ.  5) CALL MATERIAL_S5_INIT(MatID)
        IF (MType .EQ.  6) CALL MATERIAL_D6_INIT(MatID)
        IF (MType .EQ.  7) CALL MATERIAL_D7_INIT(MatID)
        IF (MType .EQ.  8) CALL MATERIAL_D8_INIT(MatID)
        IF (MType .EQ.  9) CALL MATERIAL_D9_INIT(MatID)
        IF (MType .EQ. 10) CALL MATERIAL_10_INIT(MatID)
        IF (MType .EQ. 11) CALL MATERIAL_11_INIT(MatID)
        IF (MType .EQ. 12) CALL MATERIAL_12_INIT(MatID)
        IF (MType .EQ. 17) CALL MATERIAL_17_INIT(MatID)
        IF (MType .EQ. 20) CALL MATERIAL_20_INIT(MatID)
        IF (MType .EQ. 21) CALL MATERIAL_21_INIT(MatID)
        IF (MType .EQ. 22) CALL MATERIAL_22_INIT(MatID)
        IF (MType .EQ. 27) CALL MATERIAL_27_INIT(MatID)
        IF (MType .EQ. 25) CALL MATERIAL_25_INIT(MatID)
        IF (MType .EQ. 30) CALL MATERIAL_30_INIT(MatID)
        IF (MType .EQ. 31) CALL MATERIAL_31_INIT(MatID)
        IF (MType .EQ. 32) CALL MATERIAL_32_INIT(MatID)
        IF (MType .EQ. 33) CALL MATERIAL_33_INIT(MatID)
        IF (MType .EQ. 34) CALL MATERIAL_34_INIT(MatID)
        IF (MType .EQ. 35) CALL MATERIAL_35_INIT(MatID)
        IF (MType .EQ. 36) CALL MATERIAL_36_INIT(MatID)
        IF (MType .EQ. 37) CALL MATERIAL_37_INIT(MatID)
        IF (MType .EQ. 38) CALL MATERIAL_38_INIT(MatID)
        IF (MType .EQ. 39) CALL MATERIAL_39_INIT(MatID)
        IF (MType .EQ. 40) CALL MATERIAL_40_INIT(MatID)
        IF (MType .EQ. 41) CALL MATERIAL_41_INIT(MatID)
        IF (MType .EQ. 42) CALL MATERIAL_42_INIT(MatID)
        IF (MType .EQ. 47) CALL MATERIAL_47_INIT(MatID)
        IF (MType .EQ. 45) CALL MATERIAL_45_INIT(MatID)
        IF (MType .EQ. 50) CALL MATERIAL_50_INIT(MatID)
        IF (MType .EQ. 51) CALL MATERIAL_51_INIT(MatID)
        IF (MType .EQ. 52) CALL MATERIAL_52_INIT(MatID)
        IF (MType .EQ. 55) CALL MATERIAL_55_INIT(MatID)
        IF (MType .EQ. 57) CALL MATERIAL_57_INIT(MatID)
      ENDDO
!!
!! Initialize auxillary storage.
!!
!! Hexahedron, pentahedron and tetrahedron elements.
!!
      IF (NUMHX .GT. 0) THEN
        DO N = 1,NUMHX
          MatID = HEXAH(N)%PAR%MatID
          LupID = HEXAH(N)%PAR%LupID
          Isv = HEXAH(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 30) THEN
            CALL MATERIAL_30_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 31) THEN
            CALL MATERIAL_31_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 32) THEN
            CALL MATERIAL_32_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 33) THEN
            CALL MATERIAL_33_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 34) THEN
            CALL MATERIAL_34_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 35) THEN
            CALL MATERIAL_35_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 36) THEN
            CALL MATERIAL_36_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 37) THEN
            CALL MATERIAL_37_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 38) THEN
            CALL MATERIAL_38_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 39) THEN
            CALL MATERIAL_39_INI2                                              &
     &          (HEXAH(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      IF (NUMPX .GT. 0) THEN
        DO N = 1,NUMPX
          MatID = PENTA(N)%PAR%MatID
          LupID = PENTA(N)%PAR%LupID
          Isv = PENTA(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 30) THEN
            CALL MATERIAL_30_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 31) THEN
            CALL MATERIAL_31_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 32) THEN
            CALL MATERIAL_32_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 33) THEN
            CALL MATERIAL_33_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 34) THEN
            CALL MATERIAL_34_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 35) THEN
            CALL MATERIAL_35_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 36) THEN
            CALL MATERIAL_36_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 37) THEN
            CALL MATERIAL_37_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 38) THEN
            CALL MATERIAL_38_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 39) THEN
            CALL MATERIAL_39_INI2                                              &
     &          (PENTA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      IF (NUMTX .GT. 0) THEN
        DO N = 1,NUMTX
          MatID = TETRA(N)%PAR%MatID
          LupID = TETRA(N)%PAR%LupID
          Nsv = MATERIAL(MatID)%Nsv
          Isv = TETRA(N)%PAR%Isv
          IF (MATERIAL(MatID)%Type .EQ. 30) THEN
            CALL MATERIAL_30_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 31) THEN
            CALL MATERIAL_31_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 32) THEN
            CALL MATERIAL_32_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 33) THEN
            CALL MATERIAL_33_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 34) THEN
            CALL MATERIAL_34_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 35) THEN
            CALL MATERIAL_35_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 36) THEN
            CALL MATERIAL_36_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 37) THEN
            CALL MATERIAL_37_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 38) THEN
            CALL MATERIAL_38_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 39) THEN
            CALL MATERIAL_39_INI2                                              &
     &          (TETRA(N)%RES%Stress,MatID,LupID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
!! Layered solid 8-node hexahedron and 4-node membrane sub-elements.
!!
!! Of the 40-series material models used in the LSOLD layered solid
!! LSHEX hexahedral sub-elements, (40, 41, 45, 47), only material
!! model 45 requires the state variables to be initialized. Rather
!! than write a customized "_45_INI2" module, the initialization
!! is coded here explicitly.
!!
!! Of the 20-series material models used in the LSOLD layered solid
!! LSMBQ hexahedral sub-elements, (20, 21, 22, 25, 27), only material
!! models 22 and 25 require the state variables to be initialized. Rather
!! than write a customized "_45_INI2" module, the initialization
!! is coded here explicitly.
!!
      IF (NUMLS .GT. 0) THEN
        DO N = 1,NUMLS
          LupID = LSOLD(N)%PAR%LupID
          DO i = 1,LAYERING(LupID)%Number_of_Layers
            IF (LAYERING(LupID)%Ltype(i) .EQ. 0) THEN
              Isv   = LSHEX(LSOLD(N)%PAR%ID(i))%PAR%Isv
              MatID = LSHEX(LSOLD(N)%PAR%ID(i))%PAR%MatID
              IF (MATERIAL(MatID)%Type .EQ. 45) THEN
                STATE_VARIABLES(Isv  ) = LAYERING(LupID)%Ax(i)
                STATE_VARIABLES(Isv+1) = LAYERING(LupID)%Ay(i)
                STATE_VARIABLES(Isv+2) = LAYERING(LupID)%Az(i)
                STATE_VARIABLES(Isv+3) = LAYERING(LupID)%Bx(i)
                STATE_VARIABLES(Isv+4) = LAYERING(LupID)%By(i)
                STATE_VARIABLES(Isv+5) = LAYERING(LupID)%Bz(i)
              ENDIF
            ELSE
              Isv   = LSMBQ(LSOLD(N)%PAR%ID(i))%PAR%Isv
              MatID = LSMBQ(LSOLD(N)%PAR%ID(i))%PAR%MatID
              IF (MATERIAL(MatID)%Type .EQ. 22) THEN
                STATE_VARIABLES(Isv  ) = LAYERING(LupID)%Ax(i)
                STATE_VARIABLES(Isv+1) = LAYERING(LupID)%Ay(i)
                STATE_VARIABLES(Isv+2) = LAYERING(LupID)%Az(i)
                STATE_VARIABLES(Isv+3) = LAYERING(LupID)%Bx(i)
                STATE_VARIABLES(Isv+4) = LAYERING(LupID)%By(i)
                STATE_VARIABLES(Isv+5) = LAYERING(LupID)%Bz(i)
              ELSE IF (MATERIAL(MatID)%Type .EQ. 25) THEN
                STATE_VARIABLES(Isv  ) = LAYERING(LupID)%Ax(i)
                STATE_VARIABLES(Isv+1) = LAYERING(LupID)%Ay(i)
                STATE_VARIABLES(Isv+2) = LAYERING(LupID)%Az(i)
                STATE_VARIABLES(Isv+3) = LAYERING(LupID)%Bx(i)
                STATE_VARIABLES(Isv+4) = LAYERING(LupID)%By(i)
                STATE_VARIABLES(Isv+5) = LAYERING(LupID)%Bz(i)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
!! Membrane elements, 3-node and 4-node.
!!
      IF (NUMM3 .GT. 0) THEN
        DO N = 1,NUMM3
          MatID = MEMBT(N)%PAR%MatID
          SecID = MEMBT(N)%PAR%SecID
          Isv = MEMBT(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 20) THEN
            CALL MATERIAL_20_INI2                                              &
     &          (MEMBT(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 21) THEN
            CALL MATERIAL_21_INI2                                              &
     &          (MEMBT(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 22) THEN
            CALL MATERIAL_22_INI2                                              &
     &          (MEMBT(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 27) THEN
            CALL MATERIAL_27_INI2                                              &
     &          (MEMBT(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 25) THEN
            CALL MATERIAL_25_INI2                                              &
     &          (MEMBT(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      IF (NUMM4 .GT. 0) THEN
        DO N = 1,NUMM4
          MatID = MEMBQ(N)%PAR%MatID
          SecID = MEMBQ(N)%PAR%SecID
          Isv = MEMBQ(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 20) THEN
            CALL MATERIAL_20_INI2                                              &
     &          (MEMBQ(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 21) THEN
            CALL MATERIAL_21_INI2                                              &
     &          (MEMBQ(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 22) THEN
            CALL MATERIAL_22_INI2                                              &
     &          (MEMBQ(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 27) THEN
            CALL MATERIAL_27_INI2                                              &
     &          (MEMBQ(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 25) THEN
            CALL MATERIAL_25_INI2                                              &
     &          (MEMBQ(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
!! Truss elements, 2-node.
!!
      IF (NUMTR .GT. 0) THEN
        DO N = 1,NUMTR
          MatID = TRUSS(N)%PAR%MatID
          SecID = TRUSS(N)%PAR%SecID
          Isv = TRUSS(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 10) THEN
            CALL MATERIAL_10_INI2                                              &
     &          (TRUSS(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 11) THEN
            CALL MATERIAL_11_INI2                                              &
     &          (TRUSS(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 12) THEN
            CALL MATERIAL_12_INI2                                              &
     &          (TRUSS(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 17) THEN
            CALL MATERIAL_17_INI2                                              &
     &          (TRUSS(N)%RES%Stress,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
!! Plate elements, 3-node and 4-node.
!!
      IF (NUMP3 .GT. 0) THEN
        DO N = 1,NUMP3
          NPL = N
          Ist = PLATT(N)%PAR%Ist
          Ipts = Ipts_PLATT(NPL)
          MatID = PLATT(N)%PAR%MatID
          SecID = PLATT(N)%PAR%SecID
          Isv = PLATT(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 40) THEN
            CALL MATERIAL_40_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 41) THEN
            CALL MATERIAL_41_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 42) THEN
            CALL MATERIAL_42_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 45) THEN
            CALL MATERIAL_45_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 47) THEN
            CALL MATERIAL_47_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      IF (NUMP4 .GT. 0) THEN
        DO N = 1,NUMP4
          NPL = N
          Isv = PLATQ(N)%PAR%Isv
          Ipts = Ipts_PLATQ(NPL)
          MatID = PLATQ(N)%PAR%MatID
          SecID = PLATQ(N)%PAR%SecID
          Ist = PLATQ(N)%PAR%Ist
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 40) THEN
            CALL MATERIAL_40_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 41) THEN
            CALL MATERIAL_41_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 42) THEN
            CALL MATERIAL_42_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 45) THEN
            CALL MATERIAL_45_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 47) THEN
            CALL MATERIAL_47_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
!! Beam elements.
!!
      IF (NUMBM .GT. 0) THEN
        DO N = 1,NUMBM
          Ipts = 16
          MatID = BEAM(N)%PAR%MatID
          SecID = BEAM(N)%PAR%SecID
          Isv = BEAM(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 50) THEN
            CALL MATERIAL_50_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 51) THEN
            CALL MATERIAL_51_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 52) THEN
            CALL MATERIAL_52_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 55) THEN
            CALL MATERIAL_55_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 57) THEN
            CALL MATERIAL_57_INI2                                              &
     &          (Ipts,STRESS(1,Ist),MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
!! Springs and Dampers
!!
      IF (NUMSP .GT. 0) THEN
        DO N = 1,NUMSP
          MatID = SPRING(N)%PAR%MatID
          SecID = 1
          Isv = SPRING(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 1) THEN
            CALL MATERIAL_S1_INI2                                              &
     &          (SPRING(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 2) THEN
            CALL MATERIAL_S2_INI2                                              &
     &          (SPRING(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 3) THEN
            CALL MATERIAL_S3_INI2                                              &
     &          (SPRING(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 5) THEN
            CALL MATERIAL_S5_INI2                                              &
     &          (SPRING(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      IF (NUMDM .GT. 0) THEN
        DO N = 1,NUMDM
          MatID = DAMPER(N)%PAR%MatID
          SecID = 1
          Isv = DAMPER(N)%PAR%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 6) THEN
            CALL MATERIAL_D6_INI2                                              &
     &          (DAMPER(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 7) THEN
            CALL MATERIAL_D7_INI2                                              &
     &          (DAMPER(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 8) THEN
            CALL MATERIAL_D8_INI2                                              &
     &          (DAMPER(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 9) THEN
            CALL MATERIAL_D9_INI2                                              &
     &          (DAMPER(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
!! Spring BC and Damper BC
!!
      IF (NUMSC .GT. 0) THEN
        DO N = 1,NUMSC
          MatID = SPRING_BC(N)%MatID
          SecID = 1
          Isv = SPRING_BC(N)%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 1) THEN
            CALL MATERIAL_S1_INI2                                              &
     &          (SPRING_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 2) THEN
            CALL MATERIAL_S2_INI2                                              &
     &          (SPRING_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 3) THEN
            CALL MATERIAL_S3_INI2                                              &
     &          (SPRING_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 5) THEN
            CALL MATERIAL_S5_INI2                                              &
     &          (SPRING_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      IF (NUMVC .GT. 0) THEN
        DO N = 1,NUMVC
          MatID = DAMPER_BC(N)%MatID
          SecID = 1
          Isv = DAMPER_BC(N)%Isv
          Nsv = MATERIAL(MatID)%Nsv
          IF (MATERIAL(MatID)%Type .EQ. 6) THEN
            CALL MATERIAL_D6_INI2                                              &
     &          (DAMPER_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 7) THEN
            CALL MATERIAL_D7_INI2                                              &
     &          (DAMPER_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 8) THEN
            CALL MATERIAL_D8_INI2                                              &
     &          (DAMPER_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ELSE IF (MATERIAL(MatID)%Type .EQ. 9) THEN
            CALL MATERIAL_D9_INI2                                              &
     &          (DAMPER_BC(N)%RES%Force,MatID,SecID,Isv,Nsv)
          ENDIF
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_12_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:54
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
      ENTRY MATERIAL_12_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_39_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:05:00
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          LupID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_39_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_42_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:05:01
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Ipts,                                                          &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6,Ipts)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_42_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_52_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:05:03
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Ipts,                                                          &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          Lmod,                                                          &
     &          STRESS(6,Ipts)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lame parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)
 
      Lmod = Prat * Ymod / ((ONE + Prat) * (ONE - 2.0D0*Prat))
      Gmod = Ymod / (2.0D0 + 2.0D0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H.
!!
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_52_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_39                                                   &
     &          (STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID)
!!
!! Copyright (c) by KEY Associates; 3-APR-1992 09:02:02
!!
!! Stubb MATERIAL MODEL.
!!
!! Module computes the current value of the stress STRESS(1:6), the
!! density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
!! Note: This material model does not require any state variables. The
!! array STATE_VARIABLES is present to maintain the standard material
!! interface format.
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
!! Note: Until such time that a computation of the polar decomposition
!! of the deformation gradient, F = VR = RU is implemented to obtain
!! the orthogonal rotation R, the Jaumann stress flux based on the spin
!! W will be used. The Jaumann rate is good for shear strains up to 40%.
!!
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)
!!
!! Rotate stress from the configuration at time n to the configuration at
!! time n+1/2.
!!
!!      HDT  = 0.5D0 * DTnext
!!      dWxy = HDT * Wxy
!!      dWxz = HDT * Wxz
!!      dWyz = HDT * Wyz
!!      SXX =  ((dWxy+dWxy)*STRESS(4)) + ((dWxz+dWxz)*STRESS(5))
!!      SYY = -((dWxy+dWxy)*STRESS(4)) + ((dWyz+dWyz)*STRESS(6))
!!      SZZ = -((dWxz+dWxz)*STRESS(5)) - ((dWyz+dWyz)*STRESS(6))
!!      SXY =  dWxy*(STRESS(2)-STRESS(1)) + dWxz*STRESS(6) + dWyz*STRESS(5)
!!      SXZ =  dWxz*(STRESS(3)-STRESS(1)) + dWxy*STRESS(6) - dWyz*STRESS(4)
!!      SYZ =  dWyz*(STRESS(3)-STRESS(2)) - dWxy*STRESS(5) - dWxz*STRESS(4)
!!
!!      STRESS(1) = STRESS(1) + SXX
!!      STRESS(2) = STRESS(2) + SYY
!!      STRESS(3) = STRESS(3) + SZZ
!!      STRESS(4) = STRESS(4) + SXY
!!      STRESS(5) = STRESS(5) + SXZ
!!      STRESS(6) = STRESS(6) + SYZ
!!
!! Internal energy increment from time n to time n+1/2.
!!
!!      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D0*DTnext) *
!!      2       (
!!      2        Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +
!!      2       (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +
!!      2       (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -
!!      2       (Dxx+Dyy+Dzz) * INTERNAL_ENERGY
!!      2       )
!!
!! Update stress.
!!
!!      CALL MATERIAL_30_INTEGRATION
!!      2       (
!!      2       STRESS,MatID,
!!      2       DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz
!!      2       )
!!
!! Internal energy increment from time n+1/2 to time n+1.
!!
!!      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D0*DTnext) *
!!      2       (
!!      2        Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +
!!      2       (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +
!!      2       (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -
!!      2       (Dxx+Dyy+Dzz) * INTERNAL_ENERGY
!!      2       )
!!
!! Rotate stress from the configuration at time n+1/2 to the configuration at
!! time n+1.
!!
!!      SXX =  ((dWxy+dWxy)*STRESS(4)) + ((dWxz+dWxz)*STRESS(5))
!!      SYY = -((dWxy+dWxy)*STRESS(4)) + ((dWyz+dWyz)*STRESS(6))
!!      SZZ = -((dWxz+dWxz)*STRESS(5)) - ((dWyz+dWyz)*STRESS(6))
!!      SXY =  dWxy*(STRESS(2)-STRESS(1)) + dWxz*STRESS(6) + dWyz*STRESS(5)
!!      SXZ =  dWxz*(STRESS(3)-STRESS(1)) + dWxy*STRESS(6) - dWyz*STRESS(4)
!!      SYZ =  dWyz*(STRESS(3)-STRESS(2)) - dWxy*STRESS(5) - dWxz*STRESS(4)
!!
!!      STRESS(1) = STRESS(1) + SXX
!!      STRESS(2) = STRESS(2) + SYY
!!      STRESS(3) = STRESS(3) + SZZ
!!      STRESS(4) = STRESS(4) + SXY
!!      STRESS(5) = STRESS(5) + SXZ
!!      STRESS(6) = STRESS(6) + SYZ
!!
      RETURN
      END
!!_
      SUBROUTINE LOADING_MODULI (RCL2,RCS2,P1,P2,Sinc,Einc)
!!
!! Copyright (c) by KEY Associates; 29-NOV-1991 13:11:51
!!
!! Purpose: Compute loading bulk and shear moduli from stress and strain
!! increments. If the strain increment is zero, the elastic moduli P1 and
!! P2 are used.
!!
!!        Input:  P1 = Elastic Lame parameter lambda
!!                P2 = Elastic Lame parameter mu, shear modulus
!!                Sinc = Stress increment
!!                Einc = Strain increment
!!
!!        Output: RCL2 = Longitudinal sound speed squared x density
!!                     = Bulk + 4/3 Shear
!!                RCS2 = Shear sound speed squared x density
!!                     = Shear
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(OUT) :: RCL2
      REAL(KIND(0D0)), INTENT(OUT) :: RCS2
      REAL(KIND(0D0)), INTENT(IN)  :: P1
      REAL(KIND(0D0)), INTENT(IN)  :: P2
      REAL(KIND(0D0)), INTENT(IN)  :: Sinc(6)
      REAL(KIND(0D0)), INTENT(IN)  :: Einc(6)
!!
!! Local variables.
      REAL(KIND(0D0)), PARAMETER :: OneUpon3  = (1.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: TwoUpon3  = (2.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: FourUpon3 = (4.0D+0 / 3.0D+0)
!!
!! Compute bulk strain increment.
!!
      Bulk_Strain_Inc = Einc(1) + Einc(2) + Einc(3)
!!
!! Compute mean stress increment.
!!
      Pressure_Inc = OneUpon3 * (Sinc(1) + Sinc(2) + Sinc(3))
!!
!! Compute bulk modulus.
!!
      IF (Bulk_Strain_Inc .NE. 0.0) THEN
        Bulk_Modulus = Pressure_Inc/Bulk_Strain_Inc
      ELSE
        Bulk_Modulus = P1 + TwoUpon3 * P2
      ENDIF
!!
!! Compute magnitude of deviatoric strain increment.
!!
      DEmag2 = -Bulk_Strain_Inc*Bulk_Strain_Inc
      DO i = 1,6
        DEmag2 = DEmag2 + Einc(i)*Einc(i)
      ENDDO
      DEmag2 = DEmag2 + Einc(4)*Einc(4)+Einc(5)*Einc(5)+Einc(6)*Einc(6)
!!
!! Compute shear modulus.
!!
      IF (DEmag2 .GT. 0.0) THEN
!!
!! Compute magnitude of deviatoric stress increment.
!!
        DSmag2 = -9.0D0 * Pressure_Inc * Pressure_Inc
        DO i = 1,6
          DSmag2 = DSmag2 + Sinc(i)*Sinc(i)
        ENDDO
        DSmag2 = DSmag2+Sinc(4)*Sinc(4)+Sinc(5)*Sinc(5)+Sinc(6)*Sinc(6)
!!
        IF (DSmag2 .GT. 0.0 ) THEN
          Shear_Modulus = 0.5D0 * SQRT (DSmag2/DEmag2)
        ELSE
          Shear_Modulus = P2
        ENDIF
      ELSE
        Shear_Modulus = P2
      ENDIF

      RCL2 = Bulk_Modulus + FourUpon3 * Shear_Modulus
      RCS2 = Shear_Modulus
!!
      RETURN
      END
!!_
      SUBROUTINE SYMMETRIC_TENSOR_ROTATION                                     &
     &  (TENSOR,Igenrot,Ipolard,DTnext)
!!
!! Copyright (c) by KEY Associates; 21-JUN-1992 09:45:54.59
!!
!! Purpose: Either generate rotation operator (Igenrot=1) or perform
!! rotation of symmetric, second-order tensor (Igenrot=0). Based on the
!! values of Ipolard, one of the following will occur:
!!
!! (1) Rotation used in the Jaumann stress flux based on the spin
!!     Wij = (1/2)*(Vi,j - Vj,i), (Ipolard=0).
!!
!! (2) Rotation used in the Green-McInnis stress flux based on Rashid's
!!     approximate Polar Decomposition of the deformation gradient, F = RU,
!!     and a stretching D approximating Ln(dU)/dt, (Ipolard=1).
!!
!! (3) Rotation used in the Green-McInnis stress flux based on a very
!!     nearly exact Polar Decomposition of the deformation gradient, F = RU,
!!     and a stretching D approximating Ln(dU)/dt, (Ipolard=2).
!!
!! The components of the symmetric, second-order tensor are stored as follows:
!! TENSOR(1:6) = (Txx,Tyy,Tzz,Txy,Txz,Tyz)
!!
!!
      USE shared_common_data

      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      INTEGER                                                                  &
     &        Igenrot         ! I/O Flag to indicate if this is a first
                              !     entry to generate the rotation oper-
                              !     ator (Igenrot=1) or if this is a
                              !     subsequent entry to perform the
                              !     indicated rotation (Igenrot=0).
      INTEGER                                                                  &
     &        Ipolard         ! I/- Flag to indicate if this is a Jaumann
                              !     or spin based rotation (Ipolard=0), or
                              !     or if this is an approximate D and R
                              !     based on Rashid's work (Ipolard=1), or
                              !     if this is a Polar Decomposition based
                              !     on a more exact rotation (Ipolard=2),
      REAL(KIND(0D0))                                                          &
     &        TENSOR(6)       ! I/O Tensor to be rotated.

      REAL(KIND(0D0))         &
     &        T(6),           & ! -/- Local scratch storage for TENSOR
     &        DTinc,          & ! -/- Local scratch storage for DTnext
     &        F(3,3),         & ! -/- Deformation gradient
     &        R(3,3),         & ! -/- Rotation from F = RU
     &        dWxy,dWxz,dWyz, & ! -/- Jaumann rotation quantities
     &        Txx,Tyy,Tzz,Txy,Txz,Tyz

      COMMON /SOLID/                                                           &
     &        Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),           &
     &        Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,           &
     &        Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)

      SAVE                                                                     &
     &        DTinc,R,dWxy,dWxz,dWyz
!!
!! Select rotation operator (Jaumann spin W, Polar Decomposition R, or
!! Rashid's "strongly" objective Polar Decomposition).
!!
      IF (Ipolard .EQ. 0) THEN
!!
!! JAUMANN ROTATION.
!!
        IF (Igenrot .EQ. 1) THEN
!!
          DTinc = DTnext
!!
!! Generate rotation factors.
!!
          dWxy = DTinc * Wxy
          dWxz = DTinc * Wxz
          dWyz = DTinc * Wyz
!!
!! Shut-off rotation operator generation in-case routine is called again to
!! perform additional transformations.
!!
          Igenrot = 0
        ELSE
!!
!! Put tensor in local scratch storage.
!!
          DO i = 1,6
            T(i) = TENSOR(i)
          ENDDO
!!
          Txx =  ((dWxy+dWxy)*T(4)) + ((dWxz+dWxz)*T(5))
          Tyy = -((dWxy+dWxy)*T(4)) + ((dWyz+dWyz)*T(6))
          Tzz = -((dWxz+dWxz)*T(5)) - ((dWyz+dWyz)*T(6))
          Txy =  dWxy*(T(2)-T(1)) + dWxz*T(6) + dWyz*T(5)
          Txz =  dWxz*(T(3)-T(1)) + dWxy*T(6) - dWyz*T(4)
          Tyz =  dWyz*(T(3)-T(2)) - dWxy*T(5) - dWxz*T(4)
!!
          TENSOR(1) = T(1) + Txx
          TENSOR(2) = T(2) + Tyy
          TENSOR(3) = T(3) + Tzz
          TENSOR(4) = T(4) + Txy
          TENSOR(5) = T(5) + Txz
          TENSOR(6) = T(6) + Tyz
        ENDIF
!!
      ELSE IF (Ipolard .EQ. 1) THEN
!!
!! MARK RASHID'S APPROXIMATE POLAR DECOMPOSITION ROTATION AND STRONGLY
!! OBJECTIVE STRETCHING constructed at the END of the interval. Note: The
!! rotation is to be applied at the end of the interval AFTER constitutive
!! integration.
!!
        IF (Igenrot .EQ. 1) THEN
!!
          DTinc = DTnext
!!
!! Construct A = Cinv - 1. (The purpose is to construct a "strongly"
!! objective stretching D.)
!!
          Uxx = DTinc * Vxx
          Uxy = DTinc * Vxy
          Uxz = DTinc * Vxz
          Uyx = Dtinc * Vyx
          Uyy = DTinc * Vyy
          Uyz = DTinc * Vyz
          Uzx = DTinc * Vzx
          Uzy = DTinc * Vzy
          Uzz = DTinc * Vzz

          Axx = Uxx*Uxx + Uxy*Uxy + Uxz*Uxz - Uxx - Uxx
          Axy = Uxx*Uyx + Uxy*Uyy + Uxz*Uyz - Uxy - Uyx
          Axz = Uxx*Uzx + Uxy*Uzy + Uxz*Uzz - Uxz - Uzx
          Ayy = Uyx*Uyx + Uyy*Uyy + Uyz*Uyz - Uyy - Uyy
          Ayz = Uyx*Uzx + Uyy*Uzy + Uyz*Uzz - Uyz - Uzy
          Azz = Uzx*Uzx + Uzy*Uzy + Uzz*Uzz - Uzz - Uzz
          Bxx = 0.25D0 * Axx - 0.5D0
          Bxy = 0.25D0 * Axy
          Bxz = 0.25D0 * Axz
          Byy = 0.25D0 * Ayy - 0.5D0
          Byz = 0.25D0 * Ayz
          Bzz = 0.25D0 * Azz - 0.5D0
!!
!! Re-define stretching components.
!!
          IF (DTinc .EQ. 0.0) THEN
            DTinv = ONE
          ELSE
            DTinv = ONE / DTinc
          ENDIF
          Dxx = (Bxx*Axx + Bxy*Axy + Bxz*Axz) * DTinv
          Dxy = (Bxx*Axy + Bxy*Ayy + Bxz*Ayz) * DTinv
          Dxz = (Bxx*Axz + Bxy*Ayz + Bxz*Azz) * DTinv
          Dyy = (Bxy*Axy + Byy*Ayy + Byz*Ayz) * DTinv
          Dyz = (Bxy*Axz + Byy*Ayz + Byz*Azz) * DTinv
          Dzz = (Bxz*Axz + Byz*Ayz + Bzz*Azz) * DTinv
!!
!! Construct approximate but a strictly proper orthogonal rotation.
!!
          Ax = Uyz - Uzy
          Ay = Uzx - Uxz
          Az = Uxy - Uyx
          Q = 0.25D0 * (Ax*Ax + Ay*Ay + Az*Az)
          S = 2.0D0 - (Uxx + Uyy + Uzz)
          P = 0.25D0 * S * S
          Y = ONE / ((Q+P)*(Q+P)*(Q+P))

          C1 = SQRT (P * (ONE + (P*(Q+Q+(Q+P))) * (ONE-(Q+P)) * Y))
          C2 = 0.125D0 + 0.03125D0 * Q
!!        C2 = (ONE-C1)/(Ax*Ax + Ay*Ay + Az*Az)
          C3 = 0.5D0 *                                                           &
     &      SQRT (MAX(1.0D-30,((Q+P)*(Q-P)*(ONE-P)+(P*(Q+Q+(Q+P))))*Y))

          R(1,1) = C1 + (C2*Ax)*Ax
          R(1,2) =      (C2*Ax)*Ay - (C3*Az)
          R(1,3) =      (C2*Ax)*Az + (C3*Ay)
          R(2,1) =      (C2*Ay)*Ax + (C3*Az)
          R(2,2) = C1 + (C2*Ay)*Ay
          R(2,3) =      (C2*Ay)*Az - (C3*Ax)
          R(3,1) =      (C2*Az)*Ax - (C3*Ay)
          R(3,2) =      (C2*Az)*Ay + (C3*Ax)
          R(3,3) = C1 + (C2*Az)*Az
!!
!! Save rotation operator for hourglass force rotation.
!!
          DO i = 1,9
            Rotation(i,1) = R(i,1)
          ENDDO
!!
!! Shut-off rotation operator generation. ALL subsequent calls will perform
!! transformations.
!!
          Igenrot = 0
        ELSE
!!
!! Put tensor in local scratch storage.
!!
          DO i = 1,6
            T(i) = TENSOR(i)
          ENDDO
!!
!! Rotate tensor. This is a rotation forward from time n to the end of the
!! interval and thus, uses the inverse of the rotation R just obtained. That
!! is, Rt * T * R is used in place of R * T * Rt, the transformation we would
!! have used if we had been able to construct the "forward" deformation grad-
!! ient above.
!!
!! Rotate tensor.
!!                  Rt        T        R
!!               11 21 31   1 4 5   11 12 13
!!               12 22 32 * 4 2 6 * 21 22 23
!!               13 23 33   5 6 3   31 32 33
!!
          T11 = R(1,1)*T(1) + R(2,1)*T(4) + R(3,1)*T(5)
          T12 = R(1,1)*T(4) + R(2,1)*T(2) + R(3,1)*T(6)
          T13 = R(1,1)*T(5) + R(2,1)*T(6) + R(3,1)*T(3)
          T21 = R(1,2)*T(1) + R(2,2)*T(4) + R(3,2)*T(5)
          T22 = R(1,2)*T(4) + R(2,2)*T(2) + R(3,2)*T(6)
          T23 = R(1,2)*T(5) + R(2,2)*T(6) + R(3,2)*T(3)
          T31 = R(1,3)*T(1) + R(2,3)*T(4) + R(3,3)*T(5)
          T32 = R(1,3)*T(4) + R(2,3)*T(2) + R(3,3)*T(6)
          T33 = R(1,3)*T(5) + R(2,3)*T(6) + R(3,3)*T(3)
!!
          TENSOR(1) = T11*R(1,1) + T12*R(2,1) + T13*R(3,1)
          TENSOR(4) = T11*R(1,2) + T12*R(2,2) + T13*R(3,2)
          TENSOR(5) = T11*R(1,3) + T12*R(2,3) + T13*R(3,3)
          TENSOR(2) = T21*R(1,2) + T22*R(2,2) + T23*R(3,2)
          TENSOR(6) = T21*R(1,3) + T22*R(2,3) + T23*R(3,3)
          TENSOR(3) = T31*R(1,3) + T32*R(2,3) + T33*R(3,3)
        ENDIF
!!
      ELSE IF (Ipolard .EQ. 2) THEN
!!
!! VERY NEARLY EXACT POLAR DECOMPOSITION ROTATION AND STRONGLY OBJECTIVE
!! STRETCHING constructed at the END of the interval. Note: The rotation
!! is to be applied at the end of the interval AFTER constitutive integration.
!!
        IF (Igenrot .EQ. 1) THEN
!!
          DTinc = DTnext
!!
!! Construct A = Cinv - 1. (The purpose is to construct a "strongly"
!! objective stretching D.)
!!
          Uxx = DTinc * Vxx
          Uxy = DTinc * Vxy
          Uxz = DTinc * Vxz
          Uyx = Dtinc * Vyx
          Uyy = DTinc * Vyy
          Uyz = DTinc * Vyz
          Uzx = DTinc * Vzx
          Uzy = DTinc * Vzy
          Uzz = DTinc * Vzz

          Axx = Uxx*Uxx + Uxy*Uxy + Uxz*Uxz - Uxx - Uxx
          Axy = Uxx*Uyx + Uxy*Uyy + Uxz*Uyz - Uxy - Uyx
          Axz = Uxx*Uzx + Uxy*Uzy + Uxz*Uzz - Uxz - Uzx
          Ayy = Uyx*Uyx + Uyy*Uyy + Uyz*Uyz - Uyy - Uyy
          Ayz = Uyx*Uzx + Uyy*Uzy + Uyz*Uzz - Uyz - Uzy
          Azz = Uzx*Uzx + Uzy*Uzy + Uzz*Uzz - Uzz - Uzz
          Bxx = 0.25D0 * Axx - 0.5D0
          Bxy = 0.25D0 * Axy
          Bxz = 0.25D0 * Axz
          Byy = 0.25D0 * Ayy - 0.5D0
          Byz = 0.25D0 * Ayz
          Bzz = 0.25D0 * Azz - 0.5D0
!!
!! Re-define stretching components.
!!
          IF (DTinc .EQ. 0.0) THEN
            DTinv = ONE
          ELSE
            DTinv = ONE / DTinc
          ENDIF
          Dxx = (Bxx*Axx + Bxy*Axy + Bxz*Axz) * DTinv
          Dxy = (Bxx*Axy + Bxy*Ayy + Bxz*Ayz) * DTinv
          Dxz = (Bxx*Axz + Bxy*Ayz + Bxz*Azz) * DTinv
          Dyy = (Bxy*Axy + Byy*Ayy + Byz*Ayz) * DTinv
          Dyz = (Bxy*Axz + Byy*Ayz + Byz*Azz) * DTinv
          Dzz = (Bxz*Axz + Byz*Ayz + Bzz*Azz) * DTinv
!!
!! Construct the deformation gradient F from the end of the interval (time n+1
!! if Imidint equals 0 or time n+1/2 if Imidint equals 1) to the start of the
!! interval at time n. This is, in fact, the inverse of the deformation grad-
!! ient one would naturally expect to construct. But, because we have the
!! velociy gradient based on the geometry at the "end of the interval," we
!! have to construct the inverse deformation gradient.
!!
          F(1,1) = ONE - Uxx
          F(2,2) = ONE - Uyy
          F(3,3) = ONE - Uzz
          F(1,2) = -Uxy
          F(1,3) = -Uxz
          F(2,3) = -Uyx
          F(2,1) = -Uyz
          F(3,1) = -Uzx
          F(3,2) = -Uzy
!!
!! Obtain rotation R from Polar Decomposition, F = RU.
!!
          CALL POLAR_DECOMPOSITION (F,R)
!!
!! Save rotation operator for hourglass force rotation.
!!
          DO i = 1,9
            Rotation(i,1) = R(i,1)
          ENDDO
!!
!! Shut-off rotation operator generation in-case routine is called again to
!! perform additional transformations.
!!
          Igenrot = 0
        ELSE
!!
!! Put tensor in local scratch storage.
!!
          DO i = 1,6
            T(i) = TENSOR(i)
          ENDDO
!!
!! Rotate tensor. This is a rotation forward from time n to the end of the
!! interval and thus, uses the inverse of the rotation R just obtained. That
!! is, Rt * T * R is used in place of R * T * Rt, the transformation we would
!! have used if we had been able to construct the "forward" deformation grad-
!! ient above.
!!
!! Rotate tensor.
!!                  Rt        T        R
!!               11 21 31   1 4 5   11 12 13
!!               12 22 32 * 4 2 6 * 21 22 23
!!               13 23 33   5 6 3   31 32 33
!!
          T11 = R(1,1)*T(1) + R(2,1)*T(4) + R(3,1)*T(5)
          T12 = R(1,1)*T(4) + R(2,1)*T(2) + R(3,1)*T(6)
          T13 = R(1,1)*T(5) + R(2,1)*T(6) + R(3,1)*T(3)
          T21 = R(1,2)*T(1) + R(2,2)*T(4) + R(3,2)*T(5)
          T22 = R(1,2)*T(4) + R(2,2)*T(2) + R(3,2)*T(6)
          T23 = R(1,2)*T(5) + R(2,2)*T(6) + R(3,2)*T(3)
          T31 = R(1,3)*T(1) + R(2,3)*T(4) + R(3,3)*T(5)
          T32 = R(1,3)*T(4) + R(2,3)*T(2) + R(3,3)*T(6)
          T33 = R(1,3)*T(5) + R(2,3)*T(6) + R(3,3)*T(3)
!!
          TENSOR(1) = T11*R(1,1) + T12*R(2,1) + T13*R(3,1)
          TENSOR(4) = T11*R(1,2) + T12*R(2,2) + T13*R(3,2)
          TENSOR(5) = T11*R(1,3) + T12*R(2,3) + T13*R(3,3)
          TENSOR(2) = T21*R(1,2) + T22*R(2,2) + T23*R(3,2)
          TENSOR(6) = T21*R(1,3) + T22*R(2,3) + T23*R(3,3)
          TENSOR(3) = T31*R(1,3) + T32*R(2,3) + T33*R(3,3)
        ENDIF
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE POLAR_DECOMPOSITION (F,R)  !  ,U,Uinv,DetF)
!!
!! Copyright (c) by KEY Associates; 20-JUN-1992 14:56:49.70
!!
!! Purpose: Compute the polar decomposition of the deformation gradient.
!! F = RU where R is an orthogonal rotation and U is a pure stretch.
!! The "packed" component ordering in C(1:6) is (Cxx,Cyy,Czz,Cxy,Cxz,Cyz).
!!
      USE shared_common_data

      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      REAL(KIND(0D0))       &
     &        F(3,3),       & ! I/- Deformation gradient
     &        R(3,3),       & ! -/O Rotation obtained f/ polar decomposition
     &        U(6),         & ! -/O Right stretch tensor
     &        Uinv(6),      & ! -/O Inverse of the right stretch tensor
     &        DetF,         & ! -/O Determinant of the deformation gradient
     &        C(6),         & ! -/- Right Cauchy-Green tensor, C = Ft * F
     &        CC(6),        & ! -/- Square of C, CC = C * C
     &        EV(3),        & ! -/- Eigenvalues of A
     &        I1U,I2U,I3U,  & ! -/- Invariants of U
     &        a1,b1,c1,     & ! -/- Scratch constants
     &        a2,b2,c2,d2     ! -/- Scratch constants
!!
!! Compute the right Cauchy-Green tensor, C = transpose(F) * F
!!
!!                   1 4 5   11 21 31   11 12 13
!!      C = Ft * F,  4 2 6 = 12 22 32 * 21 22 23
!!                   5 6 3   13 23 33   31 32 33
!!
      C(1) = F(1,1)*F(1,1) + F(2,1)*F(2,1) + F(3,1)*F(3,1)  !  Cxx
      C(4) = F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,1)*F(3,2)  !  Cxy
      C(5) = F(1,1)*F(1,3) + F(2,1)*F(2,3) + F(3,1)*F(3,3)  !  Cxz
      C(2) = F(1,2)*F(1,2) + F(2,2)*F(2,2) + F(3,2)*F(3,2)  !  Cyy
      C(6) = F(1,2)*F(1,3) + F(2,2)*F(2,3) + F(3,2)*F(3,3)  !  Cyz
      C(3) = F(1,3)*F(1,3) + F(2,3)*F(2,3) + F(3,3)*F(3,3)  !  Czz
!!
!! Compute the square CC of the right Cauchy-Green tensor C.
!!
!!                   1 4 5   1 4 5   1 4 5
!!      CC = C * C,  4 2 6 = 4 2 6 * 4 2 6
!!                   5 6 3   5 6 3   5 6 3
!!
      CC(1) = C(1)*C(1) + C(4)*C(4) + C(5)*C(5)  !  CCxx
      CC(4) = C(1)*C(4) + C(4)*C(2) + C(5)*C(6)  !  CCxy
      CC(5) = C(1)*C(5) + C(4)*C(6) + C(5)*C(3)  !  CCxz
      CC(2) = C(4)*C(4) + C(2)*C(2) + C(6)*C(6)  !  CCyy
      CC(6) = C(4)*C(5) + C(2)*C(6) + C(6)*C(3)  !  CCyz
      CC(3) = C(5)*C(5) + C(6)*C(6) + C(3)*C(3)  !  CCzz
!!
!! Compute the eigenvalues of C. Note: C is positive definite and should not
!! have any negative eigenvalues.
!!
      CALL EIGENVALUES (C,EV)
!!
      EV(1) = SQRT (MAX(1.0D-15,EV(1)))
      EV(2) = SQRT (MAX(1.0D-15,EV(2)))
      EV(3) = SQRT (MAX(1.0D-15,EV(3)))
!!
!! Compute the invariants of U and the determinant of F, DetF.
!!
      I1U = EV(1) + (EV(2) + EV(3))
      I2U = EV(1) * (EV(2) + EV(3)) + (EV(2) * EV(3))
      I3U = EV(1) * (EV(2) * EV(3))
!!
!!      DetF = I3U
!!
!! Compute the right stretch tensor U.
!!
!!      a1 = ONE / ((I1U*I2U) - I3U)
!!      b1 = I1U*I3U
!!      c1 = (I1U*I1U) - I2U
!!
!!      U(1) = a1 * ( b1 + c1*C(1) - CC(1) )
!!      U(2) = a1 * ( b1 + c1*C(2) - CC(2) )
!!      U(3) = a1 * ( b1 + c1*C(3) - CC(3) )
!!      U(4) = a1 * (      c1*C(4) - CC(4) )
!!      U(5) = a1 * (      c1*C(5) - CC(5) )
!!      U(6) = a1 * (      c1*C(6) - CC(6) )
!!
!! Compute the inverse of the right stretch tensor Uinv.
!!
      a2 =  ONE / (I3U*((I1U*I2U) - I3U))
      b2 =  (I1U*I2U)*I2U - I3U*((I1U*I1U) + I2U)
      c2 = -I3U - I1U*((I1U*I1U) - I2U - I2U)
      d2 =  I1U
!!
      Uinv(1) = a2 * ( b2 + c2*C(1) + d2*CC(1) )
      Uinv(2) = a2 * ( b2 + c2*C(2) + d2*CC(2) )
      Uinv(3) = a2 * ( b2 + c2*C(3) + d2*CC(3) )
      Uinv(4) = a2 * (      c2*C(4) + d2*CC(4) )
      Uinv(5) = a2 * (      c2*C(5) + d2*CC(5) )
      Uinv(6) = a2 * (      c2*C(6) + d2*CC(6) )
!!
!! Compute rotation matrix R = F * Uinv.
!!
      R(1,1) = F(1,1)*Uinv(1) + F(1,2)*Uinv(4) + F(1,3)*Uinv(5)
      R(1,2) = F(1,1)*Uinv(4) + F(1,2)*Uinv(2) + F(1,3)*Uinv(6)
      R(1,3) = F(1,1)*Uinv(5) + F(1,2)*Uinv(6) + F(1,3)*Uinv(3)
      R(2,1) = F(2,1)*Uinv(1) + F(2,2)*Uinv(4) + F(2,3)*Uinv(5)
      R(2,2) = F(2,1)*Uinv(4) + F(2,2)*Uinv(2) + F(2,3)*Uinv(6)
      R(2,3) = F(2,1)*Uinv(5) + F(2,2)*Uinv(6) + F(2,3)*Uinv(3)
      R(3,1) = F(3,1)*Uinv(1) + F(3,2)*Uinv(4) + F(3,3)*Uinv(5)
      R(3,2) = F(3,1)*Uinv(4) + F(3,2)*Uinv(2) + F(3,3)*Uinv(6)
      R(3,3) = F(3,1)*Uinv(5) + F(3,2)*Uinv(6) + F(3,3)*Uinv(3)
!!
      RETURN
      END
!!_
      SUBROUTINE EIGENVALUES (C,EV)
!!
!! Copyright (c) by KEY Associates; 20-JUN-1992 14:56:49.70
!!
!! Purpose: Compute eigenvalues of 3x3 symmetric matrix stored in packed format.
!! The "packed" component ordering in C(1:6) is (Cxx,Cyy,Czz,Cxy,Cxz,Cyz).
!!
      USE shared_common_data

      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      REAL(KIND(0D0))       &
     &        C(6),         & ! I/- Symmetric matrix in packed format
     &        EV(3),        & ! -/O Eigenvalues of C
     &        A(6),         & ! -/- Local scratch array for C
     &        B(3),         & ! -/- Local scratch array for eigenvalues
     &        JACtol,SQtol
      LOGICAL                                                                  &
     &        Converged
      DATA                                                                     &
     &        Max_Sweep /15/                                                   &
     &        JACtol    /1.0D-06/                                              &
     &        SQtol     /1.0D-12/

!!
!! Initialize eigenvalues EV and sweep parameters.
!!
      DO i = 1,6
        A(i) = C(i)
      ENDDO
      Isweep = 0
      EV(1) = A(1)
      EV(2) = A(2)
      EV(3) = A(3)
      Threshold = ONE
!!
!! Scale A to avoid problems with exponent overflow/underflow.
!!
      Amin = MIN (A(1),A(2),A(3))
      Amax = MAX (A(1),A(2),A(3))
      Iexp = -INT (0.25D0 * (LOG10(Amin) + LOG10(Amax)))
      Scale = (10.0D0)**Iexp
      B(1) = Scale
      B(2) = Scale
      B(3) = Scale
      A(1) = A(1) * Scale
      A(2) = A(2) * Scale
      A(3) = A(3) * Scale
      A(4) = A(4) * Scale
      A(5) = A(5) * Scale
      A(6) = A(6) * Scale
!!
!! Top of iteration loop.
!!
 100    Isweep = Isweep + 1
!!
      Threshold = MAX (SQtol, (1.0D-4)*Threshold)
!!
!! Work on lower triangle only (i>j). Rows are done from top to bottom.
!! Columns are done from left to right. Skip when already within tolerance.
!!
!! Row 2 and column 1.
!!
      Aratio = (A(4)*A(4))/(A(2)*A(1))
      IF (Aratio .GT. Threshold) THEN
        Abari = -B(2)*A(4)
        Abarj = -B(1)*A(4)
        Abar  =  A(2)*B(1) - A(1)*B(2)
        Arad  =  (0.5D0*Abar)*(0.5D0*Abar) + Abari*Abarj
        X = (0.5D0*Abar) + SIGN (SQRT (Arad), Abar)
        IF ( ABS(X) .LT. JACtol * MAX(ABS(Abari),ABS(Abarj)) ) THEN
          Alpha =  0.0
          Gamma = -A(4)/A(2)
        ELSE
          Alpha =  Abarj * (ONE / X)
          Gamma = -Abari * (ONE / X)
        ENDIF
        Ai = A(6)
        Aj = A(5)
        A(6) = Ai + Gamma * Aj
        A(5) = Aj + Alpha * Ai
        Aj = A(1)
        Bj = B(1)
        Ai = A(2)
        Bi = B(2)
        A(1) = Aj + (Alpha * Alpha) * Ai + (Alpha + Alpha) * A(4)
        B(1) = Bj + (Alpha * Alpha) * Bi
        A(2) = Ai + (Gamma * Gamma) * Aj + (Gamma + Gamma) * A(4)
        B(2) = Bi + (Gamma * Gamma) * Bj
        A(4) = 0.0
      ENDIF
!!
!! Row 3 and column 1.
!!
      Aratio = (A(5)*A(5))/(A(3)*A(1))
      IF (Aratio .GT. Threshold) THEN
        Abari = -B(3)*A(5)
        Abarj = -B(1)*A(5)
        Abar  =  A(3)*B(1) - A(1)*B(3)
        Arad  =  (0.5D0*Abar)*(0.5D0*Abar) + Abari*Abarj
        X = (0.5D0*Abar) + SIGN (SQRT (Arad), Abar)
        IF ( ABS(X) .LT. JACtol * MAX(ABS(Abari),ABS(Abarj)) ) THEN
          Alpha =  0.0
          Gamma = -A(5)/A(3)
        ELSE
          Alpha =  Abarj * (ONE / X)
          Gamma = -Abari * (ONE / X)
        ENDIF
        Ai = A(6)
        Aj = A(4)
        A(6) = Ai + Gamma * Aj
        A(4) = Aj + Alpha * Ai
        Aj = A(1)
        Bj = B(1)
        Ai = A(3)
        Bi = B(3)
        A(1) = Aj + (Alpha * Alpha) * Ai + (Alpha + Alpha) * A(5)
        B(1) = Bj + (Alpha * Alpha) * Bi
        A(3) = Ai + (Gamma * Gamma) * Aj + (Gamma + Gamma) * A(5)
        B(3) = Bi + (Gamma * Gamma) * Bj
        A(5) = 0.0
      ENDIF
!!
!! Row 3 and column 2.
!!
      Aratio = (A(6)*A(6))/(A(3)*A(2))
      IF (Aratio .GT. Threshold) THEN
        Abari = -B(3)*A(6)
        Abarj = -B(2)*A(6)
        Abar  =  A(3)*B(2) - A(2)*B(3)
        Arad  =  (0.5D0*Abar)*(0.5D0*Abar) + Abari*Abarj
        X = (0.5D0*Abar) + SIGN (SQRT (Arad), Abar)
        IF ( ABS(X) .LT. JACtol * MAX(ABS(Abari),ABS(Abarj)) ) THEN
          Alpha =  0.0
          Gamma = -A(6)/A(3)
        ELSE
          Alpha =  Abarj * (ONE / X)
          Gamma = -Abari * (ONE / X)
        ENDIF
        Ai = A(5)
        Aj = A(4)
        A(5) = Ai + Gamma * Aj
        A(4) = Aj + Alpha * Ai
        Aj = A(2)
        Bj = B(2)
        Ai = A(3)
        Bi = B(3)
        A(2) = Aj + (Alpha * Alpha) * Ai + (Alpha + Alpha) * A(6)
        B(2) = Bj + (Alpha * Alpha) * Bi
        A(3) = Ai + (Gamma * Gamma) * Aj + (Gamma + Gamma) * A(6)
        B(3) = Bi + (Gamma * Gamma) * Bj
        A(6) = 0.0
      ENDIF
!!
!! End of iteration loop. Update eigenvalues.
!!
      EV(1) = A(1) / B(1)
      EV(2) = A(2) / B(2)
      EV(3) = A(3) / B(3)
!!
!! Check off-diagonal entries for convergence.
!!
      Converged = (A(4)*A(4))/(A(2)*A(1)) .LE. SQtol                           &
     &      .AND. (A(5)*A(5))/(A(3)*A(1)) .LE. SQtol                           &
     &      .AND. (A(6)*A(6))/(A(3)*A(2)) .LE. SQtol

      IF (.NOT.Converged .AND. Isweep .LT. Max_Sweep) THEN
        GO TO 100
      ENDIF
!!
!! Reorder eigenvalues and return.
!!
!!      EVmin = MIN (EV(1),EV(2),EV(3))
!!      EVmax = MAX (EV(1),EV(2),EV(3))
!!      EV(2) = EV(1)+EV(2)+EV(3)-EVmin-EVmax
!!      EV(1) = EVmin
!!      EV(3) = EVmax
!!
      RETURN
      END
