!!
!! 54. Mean_Stresses (= Stress_Resultants / Area)
!!
      MODULE mean_stress_                                               
                                                                        
      TYPE :: mean_stress_type                                          
        REAL(KIND(0D0))  Nrr     ! Mean axial stress                               
        REAL(KIND(0D0))  Nrs     ! Mean shear stress                               
        REAL(KIND(0D0))  Nrt     ! Mean shear stress                               
        REAL(KIND(0D0))  Mrr     ! Mean torsional moment                           
        REAL(KIND(0D0))  Mrs     ! Mean bending moment                             
        REAL(KIND(0D0))  Mrt     ! Mean bending moment                             
      END TYPE                                                          
                                                                        
      TYPE (mean_stress_type) :: MEAN                                   
                                                                        
      END MODULE mean_stress_
