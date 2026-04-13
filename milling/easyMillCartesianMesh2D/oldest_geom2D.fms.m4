// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
include(./system/geometry_constants)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
11
(
rotor 
wall

box1_wall
wall

outlet
outletId

box2_wall 
wall 

box3_wall 
wall 

stator1_wall
wall   

stator2_wall
wall  

stator3_wall
wall  

nozzle_wall
wall  

inlet 
inletId  

axis 
wall  
)


22
(
(0   mygap 0)
(2   mygap 0)
(2.5 mygap 0)
(2.5 .5 0)
(2   .5 0)
(2 0 0)
(1 0 0)
(xstator ystator 0)
(.25 1 0)
(.25 2 0)
(0   2 0)

(0   mygap 0.1)
(2   mygap 0.1)
(2.5 mygap 0.1)
(2.5 .5 0.1)
(2 .5 0.1)
(2 0 0.1)
(1 0 0.1)
(xstator ystator 0.1)
(.25 1 0.1)
(.25 2 0.1)
(0   2 0.1)
)


22
(
((0 1 12) 0)
((0 12 11) 0)

((1 2 13) 1)
((1 13 12) 1)

((2 3 14) 2)
((2 14 13) 2)

((3 4 15) 3)
((3 15 14) 3)

((4 5 16) 4)
((4 16 15) 4)

((5 6 17) 5)
((5 17 16) 5)

((6 7 18) 6)
((6 18 17) 6)

((7 8 19) 7)
((7 19 18) 7)

((8 9 20) 8)
((8 20 19) 8)

((9 10 21) 9)
((9 21 20) 9)

((10 0 11) 10)
((10 11 21) 10)
)

0()

0
(
)


0
(
)


0
(
)
