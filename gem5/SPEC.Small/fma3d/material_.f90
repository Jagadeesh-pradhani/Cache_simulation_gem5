!!
!! 10. MATERIAL PROPERTY DATA: Element material property structure.
!!
      MODULE material_                                                  
                                                                        
      INTEGER, PARAMETER :: NPNMP = 22  ! Number of entries in array representation
                                                                        
      TYPE :: material_type                                             
       INTEGER           MatID   ! Material ID (User defined value retained)
       INTEGER           Type    ! Material type (Constitutive model number)
       CHARACTER(32)     Label   ! Material description                    
       INTEGER           Nsv     ! Number of state variables required by Type
       INTEGER           Set     ! Property name set ID                    
       REAL(KIND(0D0))   Mass    ! Total mass of this material in body     
       REAL(KIND(0D0))   Xcm     ! X-coordinate of center of mass          
       REAL(KIND(0D0))   Ycm     ! Y-coordinate of center of mass          
       REAL(KIND(0D0))   Zcm     ! Z-coordinate of center of mass          
       REAL(KIND(0D0))   Bxx     ! Component of inertia tensor             
       REAL(KIND(0D0))   Byy     ! Component of inertia tensor             
       REAL(KIND(0D0))   Bzz     ! Component of inertia tensor             
       REAL(KIND(0D0))   Bxy     ! Component of inertia tensor             
       REAL(KIND(0D0))   Bxz     ! Component of inertia tensor             
       REAL(KIND(0D0))   Byz     ! Component of inertia tensor             
       REAL(KIND(0D0))PVAL(NPNMP)! Array representation for input.         
      END TYPE                                                          
                                                                        
      TYPE (material_type), DIMENSION(:), ALLOCATABLE :: MATERIAL       
                                                                        
      END MODULE material_
