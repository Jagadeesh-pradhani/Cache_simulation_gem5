c CLASS = S
c  
c  
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
c  

c full problem size
        integer isiz1, isiz2, isiz3
        parameter (isiz1=220, isiz2=220, isiz3=220)

c number of iterations and how often to print the norm
        integer itmax_default, inorm_default
        parameter (itmax_default=350, inorm_default=50)
        double precision dt_default
        parameter (dt_default = 1.5d-03)
        logical  convertdouble
        parameter (convertdouble = .false.)
        character*3 npbversion
        parameter (npbversion='2.3')
        character*9 cs1
        parameter (cs1='')
        character*9 cs2
        parameter (cs2='')
        character*6 cs3
        parameter (cs3='(none)')
        character*6 cs4
        parameter (cs4='(none)')
        character*4 cs5
        parameter (cs5='')
        character*23 cs6
        parameter (cs6='')
        character*6 cs7
        parameter (cs7='randi8')
