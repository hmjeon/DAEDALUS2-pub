!
! =============================================================================
!
! Module - Exam_DAEDALUS2
! Last Updated : 01/28/2019, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of DAEDALUS2, which allows scientists to build and solve
! the sequence design of complex DNA nanostructures.
! Copyright 2019 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! DAEDALUS2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! DAEDALUS2 is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
module Exam_DAEDALUS2

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_DAEDALUS2_Tetrahedron           ! 1. V=4,  E=6,  F=4
    public Exam_DAEDALUS2_Cube                  ! 2. V=8,  E=12, F=6
    public Exam_DAEDALUS2_Octahedron            ! 3. V=6,  E=12, F=8
    public Exam_DAEDALUS2_Dodecahedron          ! 4. V=20, E=30, F=12
    public Exam_DAEDALUS2_Icosahedron           ! 5. V=12, E=30, F=20
    public Exam_DAEDALUS2_Penta_Bipyramid       ! V=7,  E=15, F=10

    ! Irregular; Arbitrary vertex angle but equal edge length
    public Exam_DAEDALUS2_Penta_Pyramid
    public Exam_DAEDALUS2_Penta_Cupola
    public Exam_DAEDALUS2_Prism_Hexa
    public Exam_DAEDALUS2_Anti_Prism_Penta
    public Exam_DAEDALUS2_Chiral_Object
    public Exam_DAEDALUS2_Cubeocta

    ! Irregular; Arbitrary vertex angle and edge length
    public Exam_DAEDALUS2_Asym_Tetra
    public Exam_DAEDALUS2_Twisted_Cube
    public Exam_DAEDALUS2_Twisted_Tri_Prism
    public Exam_DAEDALUS2_Bi_Propello_Tetra
    public Exam_DAEDALUS2_Bi_Propello_Cube

    ! End-open highly irregular geometries
    public Exam_DAEDALUS2_End_Cylinder_Tri
    public Exam_DAEDALUS2_Hemisphere_Quad
    public Exam_DAEDALUS2_Hemisphere_Tri
    public Exam_DAEDALUS2_Hourglass_Quad
    public Exam_DAEDALUS2_Hourglass_Tri

    ! Etc
    public Exam_DAEDALUS2_Anti_Prism_Toroid_Penta_Trapezo
    public Exam_DAEDALUS2_Isohedral_Toroid
    public Exam_DAEDALUS2_Bi_Propello_Octa
    public Exam_DAEDALUS2_Bi_Snub_Cube
    public Exam_DAEDALUS2_Bi_Penta_Icositetra
    public Exam_DAEDALUS2_Tetra_63_63_63_73_52_42
    public Exam_DAEDALUS2_Tetra_105_63_84_84_63_42
    public Exam_DAEDALUS2_Tri_Bipyramid_42_63_84_105

    ! Prism
    public Exam_DAEDALUS2_Prism_Triangle            ! V= 6, E= 9, F= 5
    public Exam_DAEDALUS2_Prism_Pentagon            ! V=10, E=15, F= 7
    public Exam_DAEDALUS2_Prism_Heptagon            ! V=14, E=21, F= 9
    public Exam_DAEDALUS2_Prism_Octagon             ! V=16, E=24, F=10
    public Exam_DAEDALUS2_Prism_Enneagon            ! V=18, E=27, F=11
    public Exam_DAEDALUS2_Prism_Decagon             ! V=20, E=30, F=12

    ! Antiprism
    public Exam_DAEDALUS2_Anti_Prism_Square         ! V= 8, E=16, F=10
    public Exam_DAEDALUS2_Anti_Prism_Hexagon        ! V=12, E=24, F=14
    public Exam_DAEDALUS2_Anti_Prism_Heptagon       ! V=14, E=28, F=16
    public Exam_DAEDALUS2_Anti_Prism_Octagon        ! V=16, E=32, F=18
    public Exam_DAEDALUS2_Anti_Prism_Enneagon       ! V=18, E=36, F=20
    public Exam_DAEDALUS2_Anti_Prism_Decagon        ! V=20, E=40, F=22

    ! Quad mesh
    public Exam_DAEDALUS2_End_Tri_Prism_Quad        ! Open end triangular prism with quad mesh
    public Exam_DAEDALUS2_End_Cube_Quad             ! Open end cube with quad mesh
    public Exam_DAEDALUS2_End_Penta_Prism_Quad      ! Open end pentagonal prism with quad mesh
    public Exam_DAEDALUS2_End_Cylinder_Quad         ! Open end cylinder with quad mesh

    ! Tri mesh
    public Exam_DAEDALUS2_End_Tri_Prism_Tri         ! Open end triangular prism with tri mesh
    public Exam_DAEDALUS2_End_Cube_Tri              ! Open end cube with tri mesh
    public Exam_DAEDALUS2_End_Penta_Prism_Tri       ! Open end pentagonal prism with tri mesh

contains

! -----------------------------------------------------------------------------

! Example of Tetrahedron
subroutine Exam_DAEDALUS2_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "01_Tetrahedron"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 4
    geom.n_face = 4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.00000d0,  0.0d0,  6.12374d0 ]
    geom.iniP(2).pos(1:3) = [  5.77351d0,  0.0d0, -2.04125d0 ]
    geom.iniP(3).pos(1:3) = [ -2.88676d0,  5.0d0, -2.04125d0 ]
    geom.iniP(4).pos(1:3) = [ -2.88676d0, -5.0d0, -2.04125d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 2, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 1, 4, 2 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 4, 3 ]

    ! Rotate the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 90.0d0)
end subroutine Exam_DAEDALUS2_Tetrahedron

! -----------------------------------------------------------------------------

! Example of Cube
subroutine Exam_DAEDALUS2_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "02_Cube"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 8
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -1.0d0, -1.0d0, -1.0d0 ]; geom.iniP(2).pos(1:3) = [  1.0d0, -1.0d0, -1.0d0 ]
    geom.iniP(3).pos(1:3) = [  1.0d0,  1.0d0, -1.0d0 ]; geom.iniP(4).pos(1:3) = [ -1.0d0,  1.0d0, -1.0d0 ]
    geom.iniP(5).pos(1:3) = [ -1.0d0, -1.0d0,  1.0d0 ]; geom.iniP(6).pos(1:3) = [  1.0d0, -1.0d0,  1.0d0 ]
    geom.iniP(7).pos(1:3) = [  1.0d0,  1.0d0,  1.0d0 ]; geom.iniP(8).pos(1:3) = [ -1.0d0,  1.0d0,  1.0d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1, 4, 3, 2 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 5, 6, 7, 8 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2, 3, 7, 6 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1, 5, 8, 4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 1, 2, 6, 5 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 3, 4, 8, 7 ]

    ! Rotate the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0], -30.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  30.0d0)
end subroutine Exam_DAEDALUS2_Cube

! -----------------------------------------------------------------------------

! Example of Octahedron
subroutine Exam_DAEDALUS2_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "03_Octahedron"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 6
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.00000d0,  0.00000d0,  7.07107d0 ]
    geom.iniP(2).pos(1:3) = [  7.07107d0,  0.00000d0,  0.00000d0 ]
    geom.iniP(3).pos(1:3) = [  0.00000d0,  7.07107d0,  0.00000d0 ]
    geom.iniP(4).pos(1:3) = [ -7.07107d0,  0.00000d0,  0.00000d0 ]
    geom.iniP(5).pos(1:3) = [  0.00000d0, -7.07107d0,  0.00000d0 ]
    geom.iniP(6).pos(1:3) = [  0.00000d0,  0.00000d0, -7.07107d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 2, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 1, 4, 5 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 1, 5, 2 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 2, 5, 6 ]
    geom.face(6).n_poi = 3; allocate(geom.face(6).poi(3)); geom.face(6).poi(1:3) = [ 2, 6, 3 ]
    geom.face(7).n_poi = 3; allocate(geom.face(7).poi(3)); geom.face(7).poi(1:3) = [ 3, 6, 4 ]
    geom.face(8).n_poi = 3; allocate(geom.face(8).poi(3)); geom.face(8).poi(1:3) = [ 4, 6, 5 ]

    ! Rotate the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0], -15.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  15.0d0)
end subroutine Exam_DAEDALUS2_Octahedron

! -----------------------------------------------------------------------------

! Example of Dodecahedron
subroutine Exam_DAEDALUS2_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "04_Dodecahedron"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  14.01264d0 ]; geom.iniP( 2).pos(1:3) = [   9.34173d0,   0.00000d0,  10.44437d0 ]
    geom.iniP( 3).pos(1:3) = [  -4.67086d0,   8.09018d0,  10.44437d0 ]; geom.iniP( 4).pos(1:3) = [  -4.67086d0,  -8.09018d0,  10.44437d0 ]
    geom.iniP( 5).pos(1:3) = [  10.44437d0,   8.09018d0,   4.67086d0 ]; geom.iniP( 6).pos(1:3) = [  10.44437d0,  -8.09018d0,   4.67086d0 ]
    geom.iniP( 7).pos(1:3) = [ -12.22848d0,   5.00000d0,   4.67086d0 ]; geom.iniP( 8).pos(1:3) = [   1.78411d0,  13.09018d0,   4.67086d0 ]
    geom.iniP( 9).pos(1:3) = [   1.78411d0, -13.09018d0,   4.67086d0 ]; geom.iniP(10).pos(1:3) = [ -12.22848d0,  -5.00000d0,   4.67086d0 ]
    geom.iniP(11).pos(1:3) = [  12.22848d0,   5.00000d0,  -4.67086d0 ]; geom.iniP(12).pos(1:3) = [  12.22848d0,  -5.00000d0,  -4.67086d0 ]
    geom.iniP(13).pos(1:3) = [ -10.44437d0,   8.09018d0,  -4.67086d0 ]; geom.iniP(14).pos(1:3) = [  -1.78411d0,  13.09018d0,  -4.67086d0 ]
    geom.iniP(15).pos(1:3) = [  -1.78411d0, -13.09018d0,  -4.67086d0 ]; geom.iniP(16).pos(1:3) = [ -10.44437d0,  -8.09018d0,  -4.67086d0 ]
    geom.iniP(17).pos(1:3) = [   4.67086d0,   8.09018d0, -10.44437d0 ]; geom.iniP(18).pos(1:3) = [   4.67086d0,  -8.09018d0, -10.44437d0 ]
    geom.iniP(19).pos(1:3) = [  -9.34173d0,   0.00000d0, -10.44437d0 ]; geom.iniP(20).pos(1:3) = [   0.00000d0,   0.00000d0, -14.01264d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi  = 5; allocate(geom.face(1).poi(5));  geom.face( 1).poi(1:5) = [  1,  2,  5,  8,  3 ]
    geom.face(2).n_poi  = 5; allocate(geom.face(2).poi(5));  geom.face( 2).poi(1:5) = [  1,  3,  7, 10,  4 ]
    geom.face(3).n_poi  = 5; allocate(geom.face(3).poi(5));  geom.face( 3).poi(1:5) = [  1,  4,  9,  6,  2 ]
    geom.face(4).n_poi  = 5; allocate(geom.face(4).poi(5));  geom.face( 4).poi(1:5) = [  2,  6, 12, 11,  5 ]
    geom.face(5).n_poi  = 5; allocate(geom.face(5).poi(5));  geom.face( 5).poi(1:5) = [  3,  8, 14, 13,  7 ]
    geom.face(6).n_poi  = 5; allocate(geom.face(6).poi(5));  geom.face( 6).poi(1:5) = [  4, 10, 16, 15,  9 ]
    geom.face(7).n_poi  = 5; allocate(geom.face(7).poi(5));  geom.face( 7).poi(1:5) = [  5, 11, 17, 14,  8 ]
    geom.face(8).n_poi  = 5; allocate(geom.face(8).poi(5));  geom.face( 8).poi(1:5) = [  6,  9, 15, 18, 12 ]
    geom.face(9).n_poi  = 5; allocate(geom.face(9).poi(5));  geom.face( 9).poi(1:5) = [  7, 13, 19, 16, 10 ]
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 11, 12, 18, 20, 17 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 13, 14, 17, 20, 19 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 15, 16, 19, 20, 18 ]
end subroutine Exam_DAEDALUS2_Dodecahedron

! -----------------------------------------------------------------------------

! Example of Icosahedron
subroutine Exam_DAEDALUS2_Icosahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "05_Icosahedron"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.00000d0,  0.00000d0,  9.51058d0 ]; geom.iniP( 2).pos(1:3) = [  8.50650d0,  0.00000d0,  4.25326d0 ]
    geom.iniP( 3).pos(1:3) = [  2.62866d0,  8.09018d0,  4.25326d0 ]; geom.iniP( 4).pos(1:3) = [ -6.88192d0,  5.00001d0,  4.25326d0 ]
    geom.iniP( 5).pos(1:3) = [ -6.88192d0, -5.00001d0,  4.25326d0 ]; geom.iniP( 6).pos(1:3) = [  2.62866d0, -8.09018d0,  4.25326d0 ]
    geom.iniP( 7).pos(1:3) = [  6.88192d0,  5.00001d0, -4.25326d0 ]; geom.iniP( 8).pos(1:3) = [  6.88192d0, -5.00001d0, -4.25326d0 ]
    geom.iniP( 9).pos(1:3) = [ -2.62866d0,  8.09018d0, -4.25326d0 ]; geom.iniP(10).pos(1:3) = [ -8.50650d0,  0.00000d0, -4.25326d0 ]
    geom.iniP(11).pos(1:3) = [ -2.62866d0, -8.09018d0, -4.25326d0 ]; geom.iniP(12).pos(1:3) = [  0.00000d0,  0.00000d0, -9.51058d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi  = 3; allocate(geom.face(1).poi(3));  geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face(2).n_poi  = 3; allocate(geom.face(2).poi(3));  geom.face( 2).poi(1:3) = [  1,  3,  4 ]
    geom.face(3).n_poi  = 3; allocate(geom.face(3).poi(3));  geom.face( 3).poi(1:3) = [  1,  4,  5 ]
    geom.face(4).n_poi  = 3; allocate(geom.face(4).poi(3));  geom.face( 4).poi(1:3) = [  1,  5,  6 ]
    geom.face(5).n_poi  = 3; allocate(geom.face(5).poi(3));  geom.face( 5).poi(1:3) = [  1,  6,  2 ]
    geom.face(6).n_poi  = 3; allocate(geom.face(6).poi(3));  geom.face( 6).poi(1:3) = [  2,  6,  8 ]
    geom.face(7).n_poi  = 3; allocate(geom.face(7).poi(3));  geom.face( 7).poi(1:3) = [  2,  8,  7 ]
    geom.face(8).n_poi  = 3; allocate(geom.face(8).poi(3));  geom.face( 8).poi(1:3) = [  2,  7,  3 ]
    geom.face(9).n_poi  = 3; allocate(geom.face(9).poi(3));  geom.face( 9).poi(1:3) = [  3,  7,  9 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  3,  9,  4 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  4,  9, 10 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  4, 10,  5 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  5, 10, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  5, 11,  6 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  6, 11,  8 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  7,  8, 12 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  7, 12,  9 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  8, 11, 12 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  9, 12, 10 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 10, 12, 11 ]
end subroutine Exam_DAEDALUS2_Icosahedron

! -----------------------------------------------------------------------------

! Example of pentagonal bipyramid(J13)
subroutine Exam_DAEDALUS2_Penta_Bipyramid(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "Penta_Bipyramid"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 7
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [-15.8832d0,   6.0527d0,  -0.7538d0 ]
    geom.iniP(2).pos(1:3) = [ -9.8887d0,  -9.9765d0,   9.6010d0 ]
    geom.iniP(3).pos(1:3) = [ -1.9321d0,  -6.1082d0,  -8.3383d0 ]
    geom.iniP(4).pos(1:3) = [  0.0720d0,  13.7173d0, -10.0665d0 ]
    geom.iniP(5).pos(1:3) = [  1.9321d0,   6.1087d0,   8.3389d0 ]
    geom.iniP(6).pos(1:3) = [  9.7727d0, -12.2186d0,   6.6868d0 ]
    geom.iniP(7).pos(1:3) = [ 15.9272d0,   2.4245d0,  -5.4681d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 4, 3, 1 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 3, 2, 1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 3, 6, 2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 1, 5, 4 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 1, 2, 5 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 5, 2, 6 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 3, 4, 7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 4, 5, 7 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 3, 7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 5, 6, 7 ]
end subroutine Exam_DAEDALUS2_Penta_Bipyramid

! -----------------------------------------------------------------------------

! Example of pentagonal pyramid(J2)
subroutine Exam_DAEDALUS2_Penta_Pyramid(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: a, angle
    integer :: i

    prob.name_prob = "06_Penta_Pyramid"
    call Mani_Set_Prob(prob, 'red')

    ! Allocate point and face structure
    geom.n_iniP = 6
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    do i = 1, 6
        geom.iniP(i).pos(:) = 0.0d0
    end do

    geom.iniP(2).pos(1) = 0.0d0
    geom.iniP(2).pos(2) = 1.0d0

    do i = 3, 6
        angle = dble(i-2) * 72.0d0 * pi / 180.0d0
        geom.iniP(i).pos(1) = cos(angle) * geom.iniP(2).pos(1) - sin(angle) * geom.iniP(2).pos(2)
        geom.iniP(i).pos(2) = sin(angle) * geom.iniP(2).pos(1) + cos(angle) * geom.iniP(2).pos(2)
    end do
    geom.iniP(1).pos(3) = sqrt( ( 5.0d0-sqrt(5.0d0) ) / 10.0d0   ) * Norm(geom.iniP(2).pos(:) - geom.iniP(3).pos(:))

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 2, 3, 1 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 3, 4, 1 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 4, 5, 1 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 5, 6, 1 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 6, 2, 1 ]
    geom.face(6).n_poi = 5; allocate(geom.face(6).poi(5)); geom.face(6).poi(1:5) = [ 6, 5, 4, 3, 2 ]

    ! Rotate the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 180.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -60.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  30.0d0)
end subroutine Exam_DAEDALUS2_Penta_Pyramid

! -----------------------------------------------------------------------------

! Example of pentagonal cupola
subroutine Exam_DAEDALUS2_Penta_Cupola(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "07_Penta_Cupola"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 15
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [ -0.973114d0,  0.120196d0, -0.57615d0  ]
    geom.iniP( 2).pos(1:3) = [ -0.844191d0, -0.563656d0, -0.512814d0 ]
    geom.iniP( 3).pos(1:3) = [ -0.711039d0,  0.75783d0,  -0.46202d0  ]
    geom.iniP( 4).pos(1:3) = [ -0.594483d0,  0.244733d0, -0.002202d0 ]
    geom.iniP( 5).pos(1:3) = [ -0.46556d0,  -0.439119d0,  0.061133d0 ]
    geom.iniP( 6).pos(1:3) = [ -0.373515d0, -1.03252d0,  -0.296206d0 ]
    geom.iniP( 7).pos(1:3) = [ -0.15807d0,   1.10569d0,  -0.21402d0  ]
    geom.iniP( 8).pos(1:3) = [ -0.041514d0,  0.592595d0,  0.245798d0 ]
    geom.iniP( 9).pos(1:3) = [  0.167087d0, -0.513901d0,  0.348277d0 ]
    geom.iniP(10).pos(1:3) = [  0.259132d0, -1.1073d0,   -0.009062d0 ]
    geom.iniP(11).pos(1:3) = [  0.429162d0,  0.123733d0,  0.462406d0 ]
    geom.iniP(12).pos(1:3) = [  0.474577d0,  1.03091d0,   0.073124d0 ]
    geom.iniP(13).pos(1:3) = [  0.812101d0, -0.759438d0,  0.238938d0 ]
    geom.iniP(14).pos(1:3) = [  0.945253d0,  0.562048d0,  0.289732d0 ]
    geom.iniP(15).pos(1:3) = [  1.07418d0,  -0.121804d0,  0.353067d0 ]

    ! Face connnectivity
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [  4,  1,  5 ] + 1
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [  8,  9, 12 ] + 1
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 10, 14, 13 ] + 1
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [  7, 11,  6 ] + 1
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [  3,  2,  0 ] + 1
    geom.face( 6).n_poi =  4; allocate(geom.face( 6).poi( 4)); geom.face( 6).poi(1: 4) = [  4,  3,  0,  1 ] + 1
    geom.face( 7).n_poi =  4; allocate(geom.face( 7).poi( 4)); geom.face( 7).poi(1: 4) = [  8,  4,  5,  9 ] + 1
    geom.face( 8).n_poi =  4; allocate(geom.face( 8).poi( 4)); geom.face( 8).poi(1: 4) = [ 10,  8, 12, 14 ] + 1
    geom.face( 9).n_poi =  4; allocate(geom.face( 9).poi( 4)); geom.face( 9).poi(1: 4) = [  7, 10, 13, 11 ] + 1
    geom.face(10).n_poi =  4; allocate(geom.face(10).poi( 4)); geom.face(10).poi(1: 4) = [  3,  7,  6,  2 ] + 1
    geom.face(11).n_poi =  5; allocate(geom.face(11).poi( 5)); geom.face(11).poi(1: 5) = [  3,  4,  8, 10,  7 ] + 1
    geom.face(12).n_poi = 10; allocate(geom.face(12).poi(10)); geom.face(12).poi(1:10) = [  2,  6, 11, 13, 14, 12, 9, 5, 1, 0 ] + 1

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],  25.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0],   5.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -45.0d0)
end subroutine Exam_DAEDALUS2_Penta_Cupola

! -----------------------------------------------------------------------------

! Example of hexagonal prism
subroutine Exam_DAEDALUS2_Prism_Hexa(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "08_Prism_Hexa"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  22.3606d0 ]
    geom.iniP( 2).pos(1:3) = [  17.8886d0,  -0.0000d0,  13.4164d0 ]
    geom.iniP( 3).pos(1:3) = [  -4.4721d0,  17.3205d0,  13.4164d0 ]
    geom.iniP( 4).pos(1:3) = [ -15.6525d0,  -8.6603d0,  13.4164d0 ]
    geom.iniP( 5).pos(1:3) = [  13.4164d0,  17.3205d0,   4.4721d0 ]
    geom.iniP( 6).pos(1:3) = [  20.1246d0,  -8.6603d0,  -4.4721d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.1246d0,   8.6603d0,   4.4721d0 ]
    geom.iniP( 8).pos(1:3) = [ -13.4164d0, -17.3205d0,  -4.4721d0 ]
    geom.iniP( 9).pos(1:3) = [  15.6525d0,   8.6603d0, -13.4164d0 ]
    geom.iniP(10).pos(1:3) = [   4.4721d0, -17.3205d0, -13.4164d0 ]
    geom.iniP(11).pos(1:3) = [ -17.8886d0,  -0.0000d0, -13.4164d0 ]
    geom.iniP(12).pos(1:3) = [   0.0000d0,  -0.0000d0, -22.3606d0 ]

    ! Position vectors
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 10, 12,  9 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 8, 11, 12, 10 ]
    geom.face(7).n_poi = 6; allocate(geom.face(7).poi(6)); geom.face(7).poi(1:6) = [ 1,  4,  8, 10,  6, 2 ]
    geom.face(8).n_poi = 6; allocate(geom.face(8).poi(6)); geom.face(8).poi(1:6) = [ 3,  5,  9, 12, 11, 7 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  27.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], -13.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],  36.5d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0],  -0.5d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  40.0d0)
end subroutine Exam_DAEDALUS2_Prism_Hexa

! -----------------------------------------------------------------------------

! Example of pentagonal antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Penta(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "09_Antiprism_Penta"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 10
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,   0.0000d0,  19.0212d0 ]
    geom.iniP( 2).pos(1:3) = [  17.0130d0,   0.0000d0,   8.5065d0 ]
    geom.iniP( 3).pos(1:3) = [   5.2573d0,  16.1804d0,   8.5065d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.7638d0,  10.0000d0,   8.5065d0 ]
    geom.iniP( 5).pos(1:3) = [ -13.7638d0, -10.0000d0,   8.5065d0 ]
    geom.iniP( 6).pos(1:3) = [  13.7638d0, -10.0000d0,  -8.5065d0 ]
    geom.iniP( 7).pos(1:3) = [  13.7638d0,  10.0000d0,  -8.5065d0 ]
    geom.iniP( 8).pos(1:3) = [ -17.0130d0,   0.0000d0,  -8.5065d0 ]
    geom.iniP( 9).pos(1:3) = [  -5.2573d0, -16.1804d0,  -8.5065d0 ]
    geom.iniP(10).pos(1:3) = [  -0.0000d0,   0.0000d0, -19.0212d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6,  9, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 10,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 10,  9 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 1,  5,  9, 6, 2 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 3,  7, 10, 8, 4 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  27.5d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], -16.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],   8.5d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  25.0d0)
end subroutine Exam_DAEDALUS2_Anti_Prism_Penta

! -----------------------------------------------------------------------------

! Example of asymmetric triangles
subroutine Exam_DAEDALUS2_Chiral_Object(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "10_Chiral_Object"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 9
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP(1).pos(1:3) = [  0.000000d0,  0.000000d0,  1.000000d0 ]
    geom.iniP(2).pos(1:3) = [  0.000000d0,  0.000000d0, -1.000000d0 ]
    geom.iniP(3).pos(1:3) = [  0.000000d0,  1.000000d0,  0.000000d0 ]
    geom.iniP(4).pos(1:3) = [  0.000000d0, -1.000000d0,  0.000000d0 ]
    geom.iniP(5).pos(1:3) = [  1.000000d0,  0.000000d0,  0.000000d0 ]
    geom.iniP(6).pos(1:3) = [ -1.000000d0,  0.000000d0,  0.000000d0 ]
    geom.iniP(7).pos(1:3) = [ -1.000000d0,  1.000000d0, -1.000000d0 ]
    geom.iniP(8).pos(1:3) = [  0.333333d0,  1.333333d0, -1.333333d0 ]
    geom.iniP(9).pos(1:3) = [  1.000000d0,  1.000000d0,  1.000000d0 ]

    ! Face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1, 5, 9 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1, 9, 3 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 9, 5, 3 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 1, 4, 5 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 1, 6, 4 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 1, 3, 6 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 2, 3, 5 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 2, 5, 4 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 2, 4, 6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 2, 6, 7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 3, 7, 6 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 3, 8, 7 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 3, 2, 8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 2, 7, 8 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0], -15.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  15.0d0)
end subroutine Exam_DAEDALUS2_Chiral_Object

! -----------------------------------------------------------------------------

! Example of Cubeoctahedron
subroutine Exam_DAEDALUS2_Cubeocta(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    prob.name_prob = "Cubeocta"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.00000d0,  0.00000d0, 10.00000d0 ]; geom.iniP( 2).pos(1:3) = [  8.66026d0,  0.00000d0,   5.00000d0 ]
    geom.iniP( 3).pos(1:3) = [  2.88675d0,  8.16497d0,  5.00000d0 ]; geom.iniP( 4).pos(1:3) = [ -8.66026d0,  0.00000d0,   5.00000d0 ]
    geom.iniP( 5).pos(1:3) = [ -2.88675d0, -8.16497d0,  5.00000d0 ]; geom.iniP( 6).pos(1:3) = [  8.66026d0,  0.00000d0,  -5.00000d0 ]
    geom.iniP( 7).pos(1:3) = [  5.77351d0, -8.16497d0,  0.00000d0 ]; geom.iniP( 8).pos(1:3) = [ -5.77351d0,  8.16497d0,   0.00000d0 ]
    geom.iniP( 9).pos(1:3) = [  2.88675d0,  8.16497d0, -5.00000d0 ]; geom.iniP(10).pos(1:3) = [ -8.66026d0,  0.00000d0,  -5.00000d0 ]
    geom.iniP(11).pos(1:3) = [ -2.88675d0, -8.16497d0, -5.00000d0 ]; geom.iniP(12).pos(1:3) = [  0.00000d0,  0.00000d0, -10.00000d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  1,  4,  5 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  7,  6 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  3,  9,  8 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  4,  8, 10 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  5, 11,  7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  6, 12,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10, 12, 11 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  1,  3,  8,  4 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  1,  5,  7,  2 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  2,  6,  9,  3 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  4, 10, 11,  5 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  6,  7, 11, 12 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  8,  9, 12, 10 ]

    ! Rotate the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],  -30.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -100.0d0)
end subroutine Exam_DAEDALUS2_Cubeocta

! -----------------------------------------------------------------------------

! Example of asymmetric tetrahedron
subroutine Exam_DAEDALUS2_Asym_Tetra(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "11_Asym_Tetra"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 4
    geom.n_face = 4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    geom.iniP(1).pos(1:3) = [ -5.9812d0,  12.5978d0, -4.5636d0 ]
    geom.iniP(2).pos(1:3) = [  0.8677d0,  -0.9451d0, 12.6382d0 ]
    geom.iniP(3).pos(1:3) = [ 11.6103d0,   1.0457d0, -4.1139d0 ]
    geom.iniP(4).pos(1:3) = [ -6.4967d0, -12.6984d0, -3.9607d0 ]

    ! Set connectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 2, 1, 4 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 2, 4, 3 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 2, 3, 1 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 4, 1, 3 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 90.0d0)
end subroutine Exam_DAEDALUS2_Asym_Tetra

! -----------------------------------------------------------------------------

! Example of the twisted cube
subroutine Exam_DAEDALUS2_Twisted_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "12_Twisted_Cube"
    call Mani_Set_Prob(prob, 'green')

    ! Allocate point and face structure
    geom.n_iniP =  9
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -11.9283d0,  -8.4001d0, -24.8515d0 ]
    geom.iniP(2).pos(1:3) = [  -2.9012d0,  13.0663d0,  21.7103d0 ]
    geom.iniP(3).pos(1:3) = [  -5.3647d0,   1.7208d0,  -1.5910d0 ]
    geom.iniP(4).pos(1:3) = [   6.8247d0, -10.0942d0, -11.4217d0 ]
    geom.iniP(5).pos(1:3) = [ -15.7149d0, -14.3456d0,  18.3808d0 ]
    geom.iniP(6).pos(1:3) = [  -6.1190d0,  18.1620d0, -12.9539d0 ]
    geom.iniP(7).pos(1:3) = [  15.0975d0,  11.4715d0, -12.2916d0 ]
    geom.iniP(8).pos(1:3) = [   5.4993d0, -16.5693d0,  12.0977d0 ]
    geom.iniP(9).pos(1:3) = [  14.6067d0,   4.9886d0,  10.9209d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 9, 2, 8 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 8, 2, 5 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 2, 3, 5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 3, 1, 5 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 6, 1, 3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 2, 6, 3 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 8, 4, 9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 4, 7, 9 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 1, 4, 5 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 5, 4, 8 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 7, 6, 9 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 9, 6, 2 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 7, 4, 1 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 6, 7, 1 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 110.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  22.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0], -25.0d0)
end subroutine Exam_DAEDALUS2_Twisted_Cube

! -----------------------------------------------------------------------------

! Example of twisted triangular prism
subroutine Exam_DAEDALUS2_Twisted_Tri_Prism(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: r1, r2, height, t_angle

    prob.name_prob = "13_Twisted_Tri_Prism"
    call Mani_Set_Prob(prob, 'green')

    ! Allocate point and face structure
    geom.n_iniP = 6
    geom.n_face = 5

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    r1      = 2.000d0
    r2      = 3.315d0
    r2      = 3.000d0
    height  = 3.500d0
    t_angle = 30.0d0

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.0d0,                       r1,                           height/2.0d0 ]
    geom.iniP(2).pos(1:3) = [ -r1*dsin( 120.d0*pi/180.0d0), r1*dcos( 120.d0*pi/180.0d0),  height/2.0d0 ]
    geom.iniP(3).pos(1:3) = [ -r1*dsin(-120.d0*pi/180.0d0), r1*dcos(-120.d0*pi/180.0d0),  height/2.0d0 ]
    geom.iniP(4).pos(1:3) = [  0.0d0,                       r2,                          -height/2.0d0 ]
    geom.iniP(5).pos(1:3) = [ -r2*dsin( 120.d0*pi/180.0d0), r2*dcos( 120.d0*pi/180.0d0), -height/2.0d0 ]
    geom.iniP(6).pos(1:3) = [ -r2*dsin(-120.d0*pi/180.0d0), r2*dcos(-120.d0*pi/180.0d0), -height/2.0d0 ]

    call Rotate_Vector(geom.iniP(4).pos, [0.0d0, 0.0d0, 1.0d0], t_angle*pi/180.0d0)
    call Rotate_Vector(geom.iniP(5).pos, [0.0d0, 0.0d0, 1.0d0], t_angle*pi/180.0d0)
    call Rotate_Vector(geom.iniP(6).pos, [0.0d0, 0.0d0, 1.0d0], t_angle*pi/180.0d0)

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 2, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 6, 5, 4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2, 5, 6, 3 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1, 3, 6, 4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 1, 4, 5, 2 ]

    !r1 = 2.0d0
    !r2 = 4.0d0
    !h  = 3.15d0
    !twist_angle = 45.0d0
    !
    !! Set point position vectors
    !geom.iniP(1).pos(1:3) = [  0.0d0,                        r1,                           h/2.0d0 ]
    !geom.iniP(2).pos(1:3) = [ -r1*dsin( 120.d0*pi/180.0d0),  r1*dcos( 120.d0*pi/180.0d0),  h/2.0d0 ]
    !geom.iniP(3).pos(1:3) = [ -r1*dsin(-120.d0*pi/180.0d0),  r1*dcos(-120.d0*pi/180.0d0),  h/2.0d0 ]
    !geom.iniP(4).pos(1:3) = [  0.0d0,                       -r2,                          -h/2.0d0 ]
    !geom.iniP(5).pos(1:3) = [  r2*dsin( 120.d0*pi/180.0d0), -r2*dcos( 120.d0*pi/180.0d0), -h/2.0d0 ]
    !geom.iniP(6).pos(1:3) = [  r2*dsin(-120.d0*pi/180.0d0), -r2*dcos(-120.d0*pi/180.0d0), -h/2.0d0 ]
    !
    !! Set point position vectors
    !geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 2, 3 ]
    !geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 6, 5, 4 ]
    !geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2, 4, 5, 3 ]
    !geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1, 3, 5, 6 ]
    !geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 2, 1, 6, 4 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -90.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],  30.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  15.0d0)
end subroutine Exam_DAEDALUS2_Twisted_Tri_Prism

! -----------------------------------------------------------------------------

! Example of chiral biscribed propello tetrahedron
subroutine Exam_DAEDALUS2_Bi_Propello_Tetra(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2, c3

    ! Set problem
    prob.name_prob = "14_Bi_Pro_Tetra"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    c0 = 0.041226573367477290238151003637d0
    c1 = 0.447594291856692697052238850001d0
    c2 = 0.577350269189625764509148780502d0
    c3 = 0.893285911422363030309232997113d0

    geom.iniP( 1).pos(1:3) = [  c1,  c0,  c3 ]
    geom.iniP( 2).pos(1:3) = [  c1, -c0, -c3 ]
    geom.iniP( 3).pos(1:3) = [ -c1, -c0,  c3 ]
    geom.iniP( 4).pos(1:3) = [ -c1,  c0, -c3 ]
    geom.iniP( 5).pos(1:3) = [  c3,  c1,  c0 ]
    geom.iniP( 6).pos(1:3) = [  c3, -c1, -c0 ]
    geom.iniP( 7).pos(1:3) = [ -c3, -c1,  c0 ]
    geom.iniP( 8).pos(1:3) = [ -c3,  c1, -c0 ]
    geom.iniP( 9).pos(1:3) = [  c0,  c3,  c1 ]
    geom.iniP(10).pos(1:3) = [  c0, -c3, -c1 ]
    geom.iniP(11).pos(1:3) = [ -c0, -c3,  c1 ]
    geom.iniP(12).pos(1:3) = [ -c0,  c3, -c1 ]
    geom.iniP(13).pos(1:3) = [  c2, -c2,  c2 ]
    geom.iniP(14).pos(1:3) = [  c2,  c2, -c2 ]
    geom.iniP(15).pos(1:3) = [ -c2,  c2,  c2 ]
    geom.iniP(16).pos(1:3) = [ -c2, -c2, -c2 ]

    ! Face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  13,  1,  3, 11 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  13, 11, 10,  6 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  13,  6,  5,  1 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  14,  2,  4, 12 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  14, 12,  9,  5 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  14,  5,  6,  2 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  15,  3,  1,  9 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  15,  9, 12,  8 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  15,  8,  7,  3 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  16,  4,  2, 10 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  16, 10, 11,  7 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  16,  7,  8,  4 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [   1,  5,  9 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [   2,  6, 10 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [   3,  7, 11 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [   4,  8, 12 ]
end subroutine Exam_DAEDALUS2_Bi_Propello_Tetra

! -----------------------------------------------------------------------------

! Example of chiral biscribed propello cube
subroutine Exam_DAEDALUS2_Bi_Propello_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2

    ! Set problem
    prob.name_prob = "15_Bi_Pro_Cube"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 32
    geom.n_face = 30

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    c0 = 0.268318503889892044924746743411d0
    c1 = 0.437902382065032587008255021493d0
    c2 = 0.649038600739057169667625426214d0

    geom.iniP( 1).pos(1:3) = [     c1,     c0,  1.0d0 ]
    geom.iniP( 2).pos(1:3) = [     c1,    -c0, -1.0d0 ]
    geom.iniP( 3).pos(1:3) = [    -c1,    -c0,  1.0d0 ]
    geom.iniP( 4).pos(1:3) = [    -c1,     c0, -1.0d0 ]
    geom.iniP( 5).pos(1:3) = [  1.0d0,     c1,     c0 ]
    geom.iniP( 6).pos(1:3) = [  1.0d0,    -c1,    -c0 ]
    geom.iniP( 7).pos(1:3) = [ -1.0d0,    -c1,     c0 ]
    geom.iniP( 8).pos(1:3) = [ -1.0d0,     c1,    -c0 ]
    geom.iniP( 9).pos(1:3) = [     c0,  1.0d0,     c1 ]
    geom.iniP(10).pos(1:3) = [     c0, -1.0d0,    -c1 ]
    geom.iniP(11).pos(1:3) = [    -c0, -1.0d0,     c1 ]
    geom.iniP(12).pos(1:3) = [    -c0,  1.0d0,    -c1 ]
    geom.iniP(13).pos(1:3) = [     c0,    -c1,  1.0d0 ]
    geom.iniP(14).pos(1:3) = [     c0,     c1, -1.0d0 ]
    geom.iniP(15).pos(1:3) = [    -c0,     c1,  1.0d0 ]
    geom.iniP(16).pos(1:3) = [    -c0,    -c1, -1.0d0 ]
    geom.iniP(17).pos(1:3) = [  1.0d0,    -c0,     c1 ]
    geom.iniP(18).pos(1:3) = [  1.0d0,     c0,    -c1 ]
    geom.iniP(19).pos(1:3) = [ -1.0d0,     c0,     c1 ]
    geom.iniP(20).pos(1:3) = [ -1.0d0,    -c0,    -c1 ]
    geom.iniP(21).pos(1:3) = [     c1, -1.0d0,     c0 ]
    geom.iniP(22).pos(1:3) = [     c1,  1.0d0,    -c0 ]
    geom.iniP(23).pos(1:3) = [    -c1,  1.0d0,     c0 ]
    geom.iniP(24).pos(1:3) = [    -c1, -1.0d0,    -c0 ]
    geom.iniP(25).pos(1:3) = [     c2,     c2,     c2 ]
    geom.iniP(26).pos(1:3) = [     c2,     c2,    -c2 ]
    geom.iniP(27).pos(1:3) = [     c2,    -c2,     c2 ]
    geom.iniP(28).pos(1:3) = [     c2,    -c2,    -c2 ]
    geom.iniP(29).pos(1:3) = [    -c2,     c2,     c2 ]
    geom.iniP(30).pos(1:3) = [    -c2,     c2,    -c2 ]
    geom.iniP(31).pos(1:3) = [    -c2,    -c2,     c2 ]
    geom.iniP(32).pos(1:3) = [    -c2,    -c2,    -c2 ]

    ! Face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  2, 12,  0, 14 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  3, 13,  1, 15 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  4, 16,  5, 17 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  7, 19,  6, 18 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  8, 21, 11, 22 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  9, 20, 10, 23 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [ 24,  0, 16,  4 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 24,  4, 21,  8 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 24,  8, 14,  0 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 25, 13, 11, 21 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 25, 21,  4, 17 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 25, 17,  1, 13 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 26, 12, 10, 20 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 26, 20,  5, 16 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 26, 16,  0, 12 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 27,  1, 17,  5 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 27,  5, 20,  9 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 27,  9, 15,  1 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 28, 14,  8, 22 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 28, 22,  7, 18 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 28, 18,  2, 14 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 29,  3, 19,  7 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 29,  7, 22, 11 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 29, 11, 13,  3 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 30,  2, 18,  6 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 30,  6, 23, 10 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 30, 10, 12,  2 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 31, 15,  9, 23 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 31, 23,  6, 19 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 31, 19,  3, 15 ]
end subroutine Exam_DAEDALUS2_Bi_Propello_Cube

! -----------------------------------------------------------------------------

! Example of open end cylinder with tri mesh
subroutine Exam_DAEDALUS2_End_Cylinder_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "16_End_Cylinder_Tri"
    call Mani_Set_Prob(prob, 'violet')

    ! Set mesh
    n  = 5
    nz = 2
    nr = 7

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr * 2 + 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face - 2
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(3) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(1) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(3) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(1) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(3) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(1) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(3) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(1) = i * nr + j
    end do

    n_face = n_face + 1
    geom.face(n_face).n_poi = 7
    allocate(geom.face(n_face).poi(7))
    geom.face(n_face).poi(1) = 1
    geom.face(n_face).poi(2) = 2
    geom.face(n_face).poi(3) = 3
    geom.face(n_face).poi(4) = 4
    geom.face(n_face).poi(5) = 5
    geom.face(n_face).poi(6) = 6
    geom.face(n_face).poi(7) = 7

    n_face = n_face + 1
    geom.face(n_face).n_poi = 7
    allocate(geom.face(n_face).poi(7))
    geom.face(n_face).poi(1) = 21
    geom.face(n_face).poi(2) = 20
    geom.face(n_face).poi(3) = 19
    geom.face(n_face).poi(4) = 18
    geom.face(n_face).poi(5) = 17
    geom.face(n_face).poi(6) = 16
    geom.face(n_face).poi(7) = 15

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0], 12.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], 30.0d0)
end subroutine Exam_DAEDALUS2_End_Cylinder_Tri

! -----------------------------------------------------------------------------

! Example of hemisphere with quad mesh
subroutine Exam_DAEDALUS2_Hemisphere_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: radius, ang1, ang2, length, ang_hole
    integer :: i, j, nz, nr, index

    ! Set problem
    prob.name_prob = "17_Hemisphere_Quad"
    call Mani_Set_Prob(prob, 'violet')

    ! Set mesh
    nz = 2
    nr = 8

    ! Set dimension
    radius   = 12.0d0
    ang_hole = 30.0d0

    ! The number of points and faces
    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    index = 0
    do i = 1, nz + 1
        ang2 = ((90.0d0 - ang_hole)-dble(i-1)*((90.0d0 - ang_hole)/dble(nz)))*(pi/180.0d0)

        do j = 1, nr
            index = index + 1
            ang1 = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(index).pos(1) = radius*cos(ang2)*cos(ang1)
            geom.iniP(index).pos(2) = radius*cos(ang2)*sin(ang1)
            geom.iniP(index).pos(3) = radius*sin(ang2)
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nz
        do j = 1, nr - 1
            index = index + 1
            geom.face(index).poi(1) = (i - 1) * nr + j
            geom.face(index).poi(2) = (i - 1) * nr + j + 1
            geom.face(index).poi(3) = i * nr + j + 1
            geom.face(index).poi(4) = i * nr + j
        end do

        index = index + 1
        geom.face(index).poi(1) = (i - 1) * nr + j
        geom.face(index).poi(2) = (i - 1) * nr + 1
        geom.face(index).poi(3) = i * nr + 1
        geom.face(index).poi(4) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -60.0d0)
end subroutine Exam_DAEDALUS2_Hemisphere_Quad

! -----------------------------------------------------------------------------

! Example of hemisphere with tri mesh
subroutine Exam_DAEDALUS2_Hemisphere_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: radius, ang1, ang2, length, ang_hole
    integer :: i, j, nz, nr, index

    ! Set problem
    prob.name_prob = "18_Hemisphere_Tri"
    call Mani_Set_Prob(prob, 'violet')

    ! Set mesh
    nz = 2
    nr = 6

    ! Set dimension
    radius   = 12.0d0
    ang_hole = 28.0d0

    ! The number of points and faces
    geom.n_face = nz * nr * 2 + 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face - 2
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    index = 0
    do i = 1, nz + 1
        ang2 = ((90.0d0 - ang_hole)-dble(i-1)*((90.0d0 - ang_hole)/dble(nz)))*(pi/180.0d0)

        do j = 1, nr
            index = index + 1
            ang1 = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(index).pos(1) = radius*cos(ang2)*cos(ang1)
            geom.iniP(index).pos(2) = radius*cos(ang2)*sin(ang1)
            geom.iniP(index).pos(3) = radius*sin(ang2)
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nz
        do j = 1, nr - 1
            index = index + 1
            geom.face(index).poi(1) = (i - 1) * nr + j
            geom.face(index).poi(2) = (i - 1) * nr + j + 1
            geom.face(index).poi(3) = i * nr + j

            index = index + 1
            geom.face(index).poi(1) = (i - 1) * nr + j + 1
            geom.face(index).poi(2) = i * nr + j + 1
            geom.face(index).poi(3) = i * nr + j
        end do

        index = index + 1
        geom.face(index).poi(1) = (i - 1) * nr + j
        geom.face(index).poi(2) = (i - 1) * nr + 1
        geom.face(index).poi(3) = i * nr + j

        index = index + 1
        geom.face(index).poi(1) = (i - 1) * nr + 1
        geom.face(index).poi(2) = i * nr + 1
        geom.face(index).poi(3) = i * nr + j
    end do

    index = index + 1
    geom.face(index).n_poi = 6
    allocate(geom.face(index).poi(6))
    geom.face(index).poi(6) = 1
    geom.face(index).poi(5) = 2
    geom.face(index).poi(4) = 3
    geom.face(index).poi(3) = 4
    geom.face(index).poi(2) = 5
    geom.face(index).poi(1) = 6

    index = index + 1
    geom.face(index).n_poi = 6
    allocate(geom.face(index).poi(6))
    geom.face(index).poi(6) = 18
    geom.face(index).poi(5) = 17
    geom.face(index).poi(4) = 16
    geom.face(index).poi(3) = 15
    geom.face(index).poi(2) = 14
    geom.face(index).poi(1) = 13

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -60.0d0)
end subroutine Exam_DAEDALUS2_Hemisphere_Tri

! -----------------------------------------------------------------------------

! Example of hourglass with quad mesh
subroutine Exam_DAEDALUS2_Hourglass_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: length, angle, radius, x_width, y_width, del_x, del_y
    integer :: i, j, jj(6), n_jnt, n_con, index, n_i_poi, n_j_poi, n, nr, nx, ny

    double precision, allocatable :: jnt(:,:)
    integer, allocatable :: con(:,:)

    !jj(:) = [1, 2, 3, 4, 3, 2]      ! "/" mesh
    jj(:) = [2, 4, 1, 3, 1, 4]      ! "\" mesh

    ! Set problem
    prob.name_prob = "19_Hourglass_Quad"
    call Mani_Set_Prob(prob, 'violet')

    nr = 5
    n  = 5

    ! Allocate memories
    allocate(con(n*nr, 4), jnt((n+1)*nr, 3))

    ! Generation of joints
    n_jnt = 0
    do i = 1, n + 1
        do j = 1, nr

            n_jnt        = n_jnt + 1
            angle        = (360.0d0 - dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)
            jnt(n_jnt,2) = 4.0d0 - dble(i-1)*(8.0d0/dble(n))
            radius       = dsqrt(1.0d0 + (jnt(n_jnt, 2) * jnt(n_jnt, 2))**(0.55d0))
            jnt(n_jnt,1) = radius * dcos(angle)
            jnt(n_jnt,3) = radius * dsin(angle)
        end do
    end do

    ! Connectivities
    con(:,:) = 0
    n_con = 0
    do i = 1, n
        do j = 1, nr - 1

            n_con = n_con + 1

            con(n_con, 1) = (i-1)*(nr) + j
            con(n_con, 2) = (i-1)*(nr) + j + 1
            con(n_con, 3) = (i-0)*(nr) + j
            con(n_con, 4) = (i-0)*(nr) + j + 1
        end do

        n_con = n_con + 1

        con(n_con, 1) = (i-1)*(nr) + nr
        con(n_con, 2) = (i-1)*(nr) + 1
        con(n_con, 3) = (i-0)*(nr) + nr
        con(n_con, 4) = (i-0)*(nr) + 1
    end do

    ! The number of points and faces
    geom.n_iniP = n_jnt
    geom.n_face = n_con
    
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    do i = 1, n_jnt
        geom.iniP(i).pos(1:3) = jnt(i,:)
    end do

    ! Set connectivity
    do i = 1, n_con
        geom.face(i).poi(1) = con(i, 1)
        geom.face(i).poi(2) = con(i, 2)
        geom.face(i).poi(3) = con(i, 4)
        geom.face(i).poi(4) = con(i, 3)
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0], 17.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], 30.0d0)
end subroutine Exam_DAEDALUS2_Hourglass_Quad

! -----------------------------------------------------------------------------

! Example of hourglass with tri mesh
subroutine Exam_DAEDALUS2_Hourglass_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: length, angle, radius, x_width, y_width, del_x, del_y
    integer :: i, j, jj(6), n_jnt, n_con, index, n_i_poi, n_j_poi, n, nr, nx, ny

    double precision, allocatable :: jnt(:,:)
    integer, allocatable :: con(:,:), out_con(:,:)

    !jj(:) = [1, 2, 3, 4, 3, 2]      ! "/" mesh
    jj(:) = [2, 4, 1, 3, 1, 4]      ! "\" mesh

    ! Set problem
    prob.name_prob = "20_Hourglass_Tri"
    call Mani_Set_Prob(prob, 'violet')

    nr = 6
    n  = 3

    ! Allocate memories
    allocate(con(n*nr, 4), jnt((n+1)*nr, 3))

    ! Generation of joints
    n_jnt = 0
    do i = 1, n + 1
        do j = 1, nr

            n_jnt        = n_jnt + 1
            angle        = (360.0d0 - dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)
            jnt(n_jnt,2) = 2.0d0 - dble(i-1)*(4.0d0/dble(n))
            radius       = dsqrt(1.0d0 + (jnt(n_jnt, 2) * jnt(n_jnt, 2))**(0.55d0))
            jnt(n_jnt,1) = radius * dcos(angle)
            jnt(n_jnt,3) = radius * dsin(angle)
        end do
    end do

    ! Connectivities
    con(:,:) = 0
    n_con = 0
    do i = 1, n
        do j = 1, nr - 1

            n_con = n_con + 1

            con(n_con, 1) = (i-1)*(nr) + j
            con(n_con, 2) = (i-1)*(nr) + j + 1
            con(n_con, 3) = (i-0)*(nr) + j
            con(n_con, 4) = (i-0)*(nr) + j + 1
        end do

        n_con = n_con + 1

        con(n_con, 1) = (i-1)*(nr) + nr
        con(n_con, 2) = (i-1)*(nr) + 1
        con(n_con, 3) = (i-0)*(nr) + nr
        con(n_con, 4) = (i-0)*(nr) + 1
    end do

    allocate(out_con(2*n_con, 3))
    call Prob_Elem_Mesh(jj, con, n_con, out_con, 3)

    ! The number of points and faces
    geom.n_iniP = n_jnt
    geom.n_face = n_con
    
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    do i = 1, n_jnt
        geom.iniP(i).pos(1:3) = jnt(i,:)
    end do

    ! Set connectivity
    do i = 1, n_con
        geom.face(i).poi(1:3) = out_con(i, 1:3)
    end do

    ! Deallocate
    deallocate(jnt, con, out_con)

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], 30.0d0)
end subroutine Exam_DAEDALUS2_Hourglass_Tri

! -----------------------------------------------------------------------------

! Example of pentagonal trapezohedron-antiprism toroid
subroutine Exam_DAEDALUS2_Anti_Prism_Toroid_Penta_Trapezo(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: C0, C1, C2, C3, C4, C5, C6

    ! Set problem
    prob.name_prob = "Anti_Toroid_Penta_Trapezo"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    C0 = 0.236067977499789696409173668731d0
    C1 = 0.618033988749894848204586834366d0
    C2 = 1.23606797749978969640917366873d0
    C3 = 1.38196601125010515179541316563d0
    C4 = 1.61803398874989484820458683437d0
    C5 = 2.23606797749978969640917366873d0
    C6 = 2.61803398874989484820458683437d0

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [     C1,    -C2,    C5 ]
    geom.iniP( 2).pos(1:3) = [     C1,    -C2,   -C5 ]
    geom.iniP( 3).pos(1:3) = [    -C1,     C2,    C5 ]
    geom.iniP( 4).pos(1:3) = [    -C1,     C2,   -C5 ]
    geom.iniP( 5).pos(1:3) = [ -2.0d0, -1.0d0,    C3 ]
    geom.iniP( 6).pos(1:3) = [ -2.0d0, -1.0d0,   -C3 ]
    geom.iniP( 7).pos(1:3) = [  2.0d0,  1.0d0,    C3 ]
    geom.iniP( 8).pos(1:3) = [  2.0d0,  1.0d0,   -C3 ]
    geom.iniP( 9).pos(1:3) = [    -C0,    -C6, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [     C0,     C6, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [  0.0d0, -1.0d0,    C6 ]
    geom.iniP(12).pos(1:3) = [  0.0d0, -1.0d0,   -C6 ]
    geom.iniP(13).pos(1:3) = [  0.0d0,  1.0d0,    C6 ]
    geom.iniP(14).pos(1:3) = [  0.0d0,  1.0d0,   -C6 ]
    geom.iniP(15).pos(1:3) = [    -C4,    -C4,    C4 ]
    geom.iniP(16).pos(1:3) = [    -C4,    -C4,   -C4 ]
    geom.iniP(17).pos(1:3) = [     C4,     C4,    C4 ]
    geom.iniP(18).pos(1:3) = [     C4,     C4,   -C4 ]
    geom.iniP(19).pos(1:3) = [ -1.0d0,    -C6, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  1.0d0,     C6, 0.0d0 ]

    ! Face connnectivity
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [ 0,  6, 16, 12, 10 ] + 1
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [ 1,  8, 18, 15, 11 ] + 1
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [ 2,  4, 14, 10, 12 ] + 1
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [ 3,  9, 19, 17, 13 ] + 1
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [ 4,  5, 15, 18, 14 ] + 1
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [ 5,  3, 13, 11, 15 ] + 1
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [ 6,  7, 17, 19, 16 ] + 1
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [ 7,  1, 11, 13, 17 ] + 1
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [ 8,  0, 10, 14, 18 ] + 1
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 9,  2, 12, 16, 19 ] + 1
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 2,  6,  0 ] + 1
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 2,  0,  4 ] + 1
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 3,  1,  7 ] + 1
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 3,  7,  9 ] + 1
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 4,  0,  8 ] + 1
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 4,  8,  5 ] + 1
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 5,  8,  1 ] + 1
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 5,  1,  3 ] + 1
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 6,  2,  9 ] + 1
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 6,  9,  7 ] + 1

    ! Set the orientation of the geometry
    !call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],  25.0d0)
    !call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0],   5.0d0)
    !call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], -45.0d0)
end subroutine Exam_DAEDALUS2_Anti_Prism_Toroid_Penta_Trapezo

! -----------------------------------------------------------------------------

! Example of isohedral toroid with 24 faces
subroutine Exam_DAEDALUS2_Isohedral_Toroid(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: C0, C1, C2, C3, C4, C5

    ! Set problem
    prob.name_prob = "Isohedral_Toroid"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    C0 = 0.103553390593273762200422181052d0
    C1 = 0.866025403784438646763723170753d0
    C2 = 1.103553390593273762200422181052d0
    C3 = 1.20710678118654752440084436210d0
    C4 = 2.09077027517602769586236520811d0
    C5 = 2.41421356237309504880168872421d0

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [  1.0d0, 0.0d0,  C2 ]
    geom.iniP( 2).pos(1:3) = [ -1.0d0, 0.0d0, -C2 ]
    geom.iniP( 3).pos(1:3) = [ -0.5d0,    C1,  C2 ]
    geom.iniP( 4).pos(1:3) = [  0.5d0,    C1, -C2 ]
    geom.iniP( 5).pos(1:3) = [ -0.5d0,   -C1,  C2 ]
    geom.iniP( 6).pos(1:3) = [  0.5d0,   -C1, -C2 ]
    geom.iniP( 7).pos(1:3) = [     C5, 0.0d0,  C0 ]
    geom.iniP( 8).pos(1:3) = [    -C5, 0.0d0, -C0 ]
    geom.iniP( 9).pos(1:3) = [    -C3,    C4,  C0 ]
    geom.iniP(10).pos(1:3) = [     C3,    C4, -C0 ]
    geom.iniP(11).pos(1:3) = [    -C3,   -C4,  C0 ]
    geom.iniP(12).pos(1:3) = [     C3,   -C4, -C0 ]

    ! Face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 0,  2,  3 ] + 1
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 0,  3,  5 ] + 1
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 0,  5,  4 ] + 1
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 0,  4, 11 ] + 1
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 0, 11,  6 ] + 1
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 0,  6,  9 ] + 1
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 0,  9,  2 ] + 1
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 1,  2,  4 ] + 1
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 1,  4,  5 ] + 1
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 1,  5, 10 ] + 1
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 1, 10,  7 ] + 1
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 1,  7,  8 ] + 1
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 1,  8,  3 ] + 1
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 1,  3,  2 ] + 1
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 2,  9,  8 ] + 1
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 2,  8,  7 ] + 1
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 2,  7,  4 ] + 1
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 3,  8,  9 ] + 1
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 3,  9,  6 ] + 1
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 3,  6,  5 ] + 1
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 4,  7, 10 ] + 1
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 4, 10, 11 ] + 1
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 5,  6, 11 ] + 1
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 5, 11, 10 ] + 1

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_Isohedral_Toroid

! -----------------------------------------------------------------------------

! Example of chiral biscribed propello octahedron
subroutine Exam_DAEDALUS2_Bi_Propello_Octa(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2

    ! Set problem
    prob.name_prob = "Bi_Pro_Octa"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 30
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    c0 = 0.149825404557484083971488566549d0
    c1 = 0.623938024555993280267936683079d0
    c2 = 0.766976981181541688825873625866d0

    geom.iniP( 1).pos(1:3) = [  0.0d0,  0.0d0,  1.0d0 ]
    geom.iniP( 2).pos(1:3) = [  0.0d0,  0.0d0, -1.0d0 ]
    geom.iniP( 3).pos(1:3) = [  1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 4).pos(1:3) = [ -1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  0.0d0,  1.0d0,  0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  0.0d0, -1.0d0,  0.0d0 ]
    geom.iniP( 7).pos(1:3) = [   c1,  -c0,   c2 ]
    geom.iniP( 8).pos(1:3) = [   c1,   c0,  -c2 ]
    geom.iniP( 9).pos(1:3) = [  -c1,   c0,   c2 ]
    geom.iniP(10).pos(1:3) = [  -c1,  -c0,  -c2 ]
    geom.iniP(11).pos(1:3) = [   c2,  -c1,   c0 ]
    geom.iniP(12).pos(1:3) = [   c2,   c1,  -c0 ]
    geom.iniP(13).pos(1:3) = [  -c2,   c1,   c0 ]
    geom.iniP(14).pos(1:3) = [  -c2,  -c1,  -c0 ]
    geom.iniP(15).pos(1:3) = [   c0,  -c2,   c1 ]
    geom.iniP(16).pos(1:3) = [   c0,   c2,  -c1 ]
    geom.iniP(17).pos(1:3) = [  -c0,   c2,   c1 ]
    geom.iniP(18).pos(1:3) = [  -c0,  -c2,  -c1 ]
    geom.iniP(19).pos(1:3) = [   c0,   c1,   c2 ]
    geom.iniP(20).pos(1:3) = [   c0,  -c1,  -c2 ]
    geom.iniP(21).pos(1:3) = [  -c0,  -c1,   c2 ]
    geom.iniP(22).pos(1:3) = [  -c0,   c1,  -c2 ]
    geom.iniP(23).pos(1:3) = [   c2,   c0,   c1 ]
    geom.iniP(24).pos(1:3) = [   c2,  -c0,  -c1 ]
    geom.iniP(25).pos(1:3) = [  -c2,  -c0,   c1 ]
    geom.iniP(26).pos(1:3) = [  -c2,   c0,  -c1 ]
    geom.iniP(27).pos(1:3) = [   c1,   c2,   c0 ]
    geom.iniP(28).pos(1:3) = [   c1,  -c2,  -c0 ]
    geom.iniP(29).pos(1:3) = [  -c1,  -c2,   c0 ]
    geom.iniP(30).pos(1:3) = [  -c1,   c2,  -c0 ]

    ! Set connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  0,  6, 22, 18 ] + 1
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  0, 18, 16,  8 ] + 1
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  0,  8, 24, 20 ] + 1
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  0, 20, 14,  6 ] + 1
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  1,  7, 23, 19 ] + 1
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  1, 19, 17,  9 ] + 1
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  1,  9, 25, 21 ] + 1
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  1, 21, 15,  7 ] + 1
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  2, 10, 27, 23 ] + 1
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  2, 23,  7, 11 ] + 1
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  2, 11, 26, 22 ] + 1
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  2, 22,  6, 10 ] + 1
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  3, 12, 29, 25 ] + 1
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  3, 25,  9, 13 ] + 1
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  3, 13, 28, 24 ] + 1
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  3, 24,  8, 12 ] + 1
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [  4, 15, 21, 29 ] + 1
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  4, 29, 12, 16 ] + 1
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [  4, 16, 18, 26 ] + 1
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [  4, 26, 11, 15 ] + 1
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [  5, 14, 20, 28 ] + 1
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [  5, 28, 13, 17 ] + 1
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [  5, 17, 19, 27 ] + 1
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [  5, 27, 10, 14 ] + 1
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 14, 10,  6 ] + 1
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 15, 11,  7 ] + 1
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 16, 12,  8 ] + 1
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 17, 13,  9 ] + 1
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 18, 22, 26 ] + 1
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 19, 23, 27 ] + 1
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 20, 24, 28 ] + 1
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 21, 25, 29 ] + 1
end subroutine Exam_DAEDALUS2_Bi_Propello_Octa

! -----------------------------------------------------------------------------

! Example of chiral biscribed snub cube
subroutine Exam_DAEDALUS2_Bi_Snub_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2

    ! Set problem
    prob.name_prob = "Bi_Snub_Cube"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 24
    geom.n_face = 38

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    c0 = 0.33775397381375235213753224516503d0
    c1 = 0.62122641055658531169250095449000d0
    c2 = 1.14261350892596209347948408672000d0

    geom.iniP( 1).pos(1:3) = [  c1,  c0,  c2 ]
    geom.iniP( 2).pos(1:3) = [  c1, -c0, -c2 ]
    geom.iniP( 3).pos(1:3) = [ -c1, -c0,  c2 ]
    geom.iniP( 4).pos(1:3) = [ -c1,  c0, -c2 ]
    geom.iniP( 5).pos(1:3) = [  c2,  c1,  c0 ]
    geom.iniP( 6).pos(1:3) = [  c2, -c1, -c0 ]
    geom.iniP( 7).pos(1:3) = [ -c2, -c1,  c0 ]
    geom.iniP( 8).pos(1:3) = [ -c2,  c1, -c0 ]
    geom.iniP( 9).pos(1:3) = [  c0,  c2,  c1 ]
    geom.iniP(10).pos(1:3) = [  c0, -c2, -c1 ]
    geom.iniP(11).pos(1:3) = [ -c0, -c2,  c1 ]
    geom.iniP(12).pos(1:3) = [ -c0,  c2, -c1 ]
    geom.iniP(13).pos(1:3) = [  c0, -c1,  c2 ]
    geom.iniP(14).pos(1:3) = [  c0,  c1, -c2 ]
    geom.iniP(15).pos(1:3) = [ -c0,  c1,  c2 ]
    geom.iniP(16).pos(1:3) = [ -c0, -c1, -c2 ]
    geom.iniP(17).pos(1:3) = [  c2, -c0,  c1 ]
    geom.iniP(18).pos(1:3) = [  c2,  c0, -c1 ]
    geom.iniP(19).pos(1:3) = [ -c2,  c0,  c1 ]
    geom.iniP(20).pos(1:3) = [ -c2, -c0, -c1 ]
    geom.iniP(21).pos(1:3) = [  c1, -c2,  c0 ]
    geom.iniP(22).pos(1:3) = [  c1,  c2, -c0 ]
    geom.iniP(23).pos(1:3) = [ -c1,  c2,  c0 ]
    geom.iniP(24).pos(1:3) = [ -c1, -c2, -c0 ]

    ! Set connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  2, 12,  0, 14 ] + 1
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  3, 13,  1, 15 ] + 1
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  4, 16,  5, 17 ] + 1
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  7, 19,  6, 18 ] + 1
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  8, 21, 11, 22 ] + 1
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  9, 20, 10, 23 ] + 1
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  0,  8, 14 ] + 1
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  1,  9, 15 ] + 1
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  2, 10, 12 ] + 1
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  3, 11, 13 ] + 1
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  4,  0, 16 ] + 1
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  5,  1, 17 ] + 1
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  6,  2, 18 ] + 1
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  7,  3, 19 ] + 1
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  8,  4, 21 ] + 1
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  9,  5, 20 ] + 1
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 10,  6, 23 ] + 1
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 11,  7, 22 ] + 1
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 12, 16,  0 ] + 1
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 13, 17,  1 ] + 1
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 14, 18,  2 ] + 1
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 15, 19,  3 ] + 1
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 16, 20,  5 ] + 1
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 17, 21,  4 ] + 1
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 18, 22,  7 ] + 1
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 19, 23,  6 ] + 1
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 20, 12, 10 ] + 1
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 21, 13, 11 ] + 1
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 22, 14,  8 ] + 1
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 23, 15,  9 ] + 1
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [  8,  0,  4 ] + 1
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [  9,  1,  5 ] + 1
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 10,  2,  6 ] + 1
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 11,  3,  7 ] + 1
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 12, 20, 16 ] + 1
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 13, 21, 17 ] + 1
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 14, 22, 18 ] + 1
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 15, 23, 19 ] + 1
end subroutine Exam_DAEDALUS2_Bi_Snub_Cube

! -----------------------------------------------------------------------------

! Example of chiral biscribed pentagonal icositetrahedron
subroutine Exam_DAEDALUS2_Bi_Penta_Icositetra(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2, c3

    ! Set problem
    prob.name_prob = "Bi_Penta_Icositetra"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 38
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    c0 = 0.107540025570073586711065317103d0
    c1 = 0.577350269189625764509148780502d0
    c2 = 0.643438410432318999656172634091d0
    c3 = 0.757906428842451843370418995303d0

    geom.iniP( 1).pos(1:3) = [  0.0d0,  0.0d0,  1.0d0 ]
    geom.iniP( 2).pos(1:3) = [  0.0d0,  0.0d0, -1.0d0 ]
    geom.iniP( 3).pos(1:3) = [  1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 4).pos(1:3) = [ -1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  0.0d0,  1.0d0,  0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  0.0d0, -1.0d0,  0.0d0 ]
    geom.iniP( 7).pos(1:3) = [   c2,  -c0,   c3 ]
    geom.iniP( 8).pos(1:3) = [   c2,   c0,  -c3 ]
    geom.iniP( 9).pos(1:3) = [  -c2,   c0,   c3 ]
    geom.iniP(10).pos(1:3) = [  -c2,  -c0,  -c3 ]
    geom.iniP(11).pos(1:3) = [   c3,  -c2,   c0 ]
    geom.iniP(12).pos(1:3) = [   c3,   c2,  -c0 ]
    geom.iniP(13).pos(1:3) = [  -c3,   c2,   c0 ]
    geom.iniP(14).pos(1:3) = [  -c3,  -c2,  -c0 ]
    geom.iniP(15).pos(1:3) = [   c0,  -c3,   c2 ]
    geom.iniP(16).pos(1:3) = [   c0,   c3,  -c2 ]
    geom.iniP(17).pos(1:3) = [  -c0,   c3,   c2 ]
    geom.iniP(18).pos(1:3) = [  -c0,  -c3,  -c2 ]
    geom.iniP(19).pos(1:3) = [   c0,   c2,   c3 ]
    geom.iniP(20).pos(1:3) = [   c0,  -c2,  -c3 ]
    geom.iniP(21).pos(1:3) = [  -c0,  -c2,   c3 ]
    geom.iniP(22).pos(1:3) = [  -c0,   c2,  -c3 ]
    geom.iniP(23).pos(1:3) = [   c3,   c0,   c2 ]
    geom.iniP(24).pos(1:3) = [   c3,  -c0,  -c2 ]
    geom.iniP(25).pos(1:3) = [  -c3,  -c0,   c2 ]
    geom.iniP(26).pos(1:3) = [  -c3,   c0,  -c2 ]
    geom.iniP(27).pos(1:3) = [   c2,   c3,   c0 ]
    geom.iniP(28).pos(1:3) = [   c2,  -c3,  -c0 ]
    geom.iniP(29).pos(1:3) = [  -c2,  -c3,   c0 ]
    geom.iniP(30).pos(1:3) = [  -c2,   c3,  -c0 ]
    geom.iniP(31).pos(1:3) = [   c1,   c1,   c1 ]
    geom.iniP(32).pos(1:3) = [   c1,   c1,  -c1 ]
    geom.iniP(33).pos(1:3) = [   c1,  -c1,   c1 ]
    geom.iniP(34).pos(1:3) = [   c1,  -c1,  -c1 ]
    geom.iniP(35).pos(1:3) = [  -c1,   c1,   c1 ]
    geom.iniP(36).pos(1:3) = [  -c1,   c1,  -c1 ]
    geom.iniP(37).pos(1:3) = [  -c1,  -c1,   c1 ]
    geom.iniP(38).pos(1:3) = [  -c1,  -c1,  -c1 ]

    ! Set connectivity
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [ 0,  6, 22, 30, 18 ] + 1
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [ 0, 18, 16, 34,  8 ] + 1
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [ 0,  8, 24, 36, 20 ] + 1
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [ 0, 20, 14, 32,  6 ] + 1
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [ 1,  7, 23, 33, 19 ] + 1
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [ 1, 19, 17, 37,  9 ] + 1
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [ 1,  9, 25, 35, 21 ] + 1
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [ 1, 21, 15, 31,  7 ] + 1
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [ 2, 10, 27, 33, 23 ] + 1
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 2, 23,  7, 31, 11 ] + 1
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 2, 11, 26, 30, 22 ] + 1
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 2, 22,  6, 32, 10 ] + 1
    geom.face(13).n_poi = 5; allocate(geom.face(13).poi(5)); geom.face(13).poi(1:5) = [ 3, 12, 29, 35, 25 ] + 1
    geom.face(14).n_poi = 5; allocate(geom.face(14).poi(5)); geom.face(14).poi(1:5) = [ 3, 25,  9, 37, 13 ] + 1
    geom.face(15).n_poi = 5; allocate(geom.face(15).poi(5)); geom.face(15).poi(1:5) = [ 3, 13, 28, 36, 24 ] + 1
    geom.face(16).n_poi = 5; allocate(geom.face(16).poi(5)); geom.face(16).poi(1:5) = [ 3, 24,  8, 34, 12 ] + 1
    geom.face(17).n_poi = 5; allocate(geom.face(17).poi(5)); geom.face(17).poi(1:5) = [ 4, 15, 21, 35, 29 ] + 1
    geom.face(18).n_poi = 5; allocate(geom.face(18).poi(5)); geom.face(18).poi(1:5) = [ 4, 29, 12, 34, 16 ] + 1
    geom.face(19).n_poi = 5; allocate(geom.face(19).poi(5)); geom.face(19).poi(1:5) = [ 4, 16, 18, 30, 26 ] + 1
    geom.face(20).n_poi = 5; allocate(geom.face(20).poi(5)); geom.face(20).poi(1:5) = [ 4, 26, 11, 31, 15 ] + 1
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [ 5, 14, 20, 36, 28 ] + 1
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [ 5, 28, 13, 37, 17 ] + 1
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [ 5, 17, 19, 33, 27 ] + 1
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [ 5, 27, 10, 32, 14 ] + 1
end subroutine Exam_DAEDALUS2_Bi_Penta_Icositetra

! -----------------------------------------------------------------------------

! Example of asymmetric tetrahedron 63-63-63-73-52-42bp
subroutine Exam_DAEDALUS2_Tetra_63_63_63_73_52_42(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Asym_Tetra_63_63_63_73_52_42"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 4
    geom.n_face = 4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP(1).pos(1:3) = [  0.000000d0,  0.000000d0,  0.000000d0 ]
    geom.iniP(2).pos(1:3) = [  0.000000d0, 63.000000d0,  0.000000d0 ]
    geom.iniP(3).pos(1:3) = [ 59.499976d0, 20.706349d0,  0.000000d0 ]
    geom.iniP(4).pos(1:3) = [ 37.426315d0, 41.539683d0, 29.029739d0 ]

    ! Face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 4, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 1, 2, 3 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 4, 3 ]
end subroutine Exam_DAEDALUS2_Tetra_63_63_63_73_52_42

! -----------------------------------------------------------------------------

! Example of asymmetric tetrahedron 105-63-84-84-63-42bp
subroutine Exam_DAEDALUS2_Tetra_105_63_84_84_63_42(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Asym_Tetra_105_63_84_84_63_42"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 4
    geom.n_face = 4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! 105-63-84-84-63-42
    geom.iniP(1).pos(1:3) = [  0.000000d0,   0.000000d0,  0.000000d0 ]
    geom.iniP(2).pos(1:3) = [  0.000000d0, 105.000000d0,  0.000000d0 ]
    geom.iniP(3).pos(1:3) = [ 50.400000d0,  37.800000d0,  0.000000d0 ]
    geom.iniP(4).pos(1:3) = [ 41.475000d0,  67.200000d0, 28.635369d0 ]

    ! Face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 4, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 1, 2, 3 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 4, 3 ]
end subroutine Exam_DAEDALUS2_Tetra_105_63_84_84_63_42

! -----------------------------------------------------------------------------

! Example of asymetric triangular bipyramid 42-63-84-105
subroutine Exam_DAEDALUS2_Tri_Bipyramid_42_63_84_105(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Asym_Tri_Bipyramid_42_63_84_105"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 5
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! 105-63-84-84-63-42
    geom.iniP(1).pos(1:3) = [  0.000000d0,   0.000000d0,  0.000000d0 ]
    geom.iniP(2).pos(1:3) = [  0.000000d0, 105.000000d0,  0.000000d0 ]
    geom.iniP(3).pos(1:3) = [ 50.400000d0,  37.800000d0,  0.000000d0 ]
    geom.iniP(4).pos(1:3) = [ 41.475000d0,  67.200000d0, 28.635369d0 ]
    !geom.iniP(5).pos(1:3) = [ -119.0d0/40.0d0, 273.0d0/10.0d0, -( 7.0d0*1319.0d0**0.5d0)/8.0d0 ]    ! 84-63-42
    geom.iniP(5).pos(1:3) = [  189.0d0/10.0d0, 273.0d0/10.0d0, -(21.0d0*   6.0d0**0.5d0)/2.0d0 ]   ! 84-42-42
    !geom.iniP(5).pos(1:3) = [  329.0d0/10.0d0, 189.0d0/ 5.0d0, -( 7.0d0* 119.0d0**0.5d0)/2.0d0 ]    ! 84-42-63

    ! Face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 4, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 2, 4, 3 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 5, 1 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 5, 3, 1 ]
    geom.face(6).n_poi = 3; allocate(geom.face(6).poi(3)); geom.face(6).poi(1:3) = [ 3, 5, 2 ]
end subroutine Exam_DAEDALUS2_Tri_Bipyramid_42_63_84_105

! -----------------------------------------------------------------------------

! Example of triangular prism
subroutine Exam_DAEDALUS2_Prism_Triangle(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Prism_Triangle"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 6
    geom.n_face = 5

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP(1).pos(1:3) = [  -0.0000d0,   0.0000d0,  15.2753d0 ]
    geom.iniP(2).pos(1:3) = [  15.1186d0,   0.0000d0,   2.1822d0 ]
    geom.iniP(3).pos(1:3) = [ -11.3389d0,  10.0000d0,   2.1822d0 ]
    geom.iniP(4).pos(1:3) = [   1.8898d0, -15.0000d0,   2.1822d0 ]
    geom.iniP(5).pos(1:3) = [   3.7796d0,  10.0000d0, -10.9109d0 ]
    geom.iniP(6).pos(1:3) = [  -9.4491d0,  -5.0000d0, -10.9109d0 ]

    ! Position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 4, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 3, 5, 6 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 1, 2, 5, 3 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1, 3, 6, 4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 2, 4, 6, 5 ]
end subroutine Exam_DAEDALUS2_Prism_Triangle

! -----------------------------------------------------------------------------

! Example of pentagonal prism
subroutine Exam_DAEDALUS2_Prism_Pentagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Prism_Pentagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 10
    geom.n_face = 7

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,  -0.0000d0,  19.7343d0 ]
    geom.iniP( 2).pos(1:3) = [  17.2421d0,  -0.0000d0,   9.5997d0 ]
    geom.iniP( 3).pos(1:3) = [  -5.9570d0,  16.1804d0,   9.5997d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.1260d0, -11.1804d0,   9.5997d0 ]
    geom.iniP( 5).pos(1:3) = [  11.2851d0,  16.1804d0,  -0.5350d0 ]
    geom.iniP( 6).pos(1:3) = [  14.7724d0, -11.1804d0,  -6.7985d0 ]
    geom.iniP( 7).pos(1:3) = [ -19.0829d0,   5.0000d0,  -0.5350d0 ]
    geom.iniP( 8).pos(1:3) = [  -3.9961d0, -18.0902d0,  -6.7985d0 ]
    geom.iniP( 9).pos(1:3) = [   8.8154d0,   5.0000d0, -16.9332d0 ]
    geom.iniP(10).pos(1:3) = [  -9.9531d0,  -1.9098d0, -16.9332d0 ]

    ! Position vectors
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1, 2,  5,  3 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 1, 3,  7,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2, 6,  9,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4, 7, 10,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 8, 10,  9 ]
    geom.face(6).n_poi = 5; allocate(geom.face(6).poi(5)); geom.face(6).poi(1:5) = [ 1, 4,  8,  6, 2 ]
    geom.face(7).n_poi = 5; allocate(geom.face(7).poi(5)); geom.face(7).poi(1:5) = [ 3, 5,  9, 10, 7 ]
end subroutine Exam_DAEDALUS2_Prism_Pentagon

! -----------------------------------------------------------------------------

! Example of hetagonal prism
subroutine Exam_DAEDALUS2_Prism_Heptagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Prism_Heptagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 9

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,  -0.0000d0,  25.1236d0 ]
    geom.iniP( 2).pos(1:3) = [  18.3475d0,  -0.0000d0,  17.1630d0 ]
    geom.iniP( 3).pos(1:3) = [  -3.4540d0,  18.0194d0,  17.1630d0 ]
    geom.iniP( 4).pos(1:3) = [ -17.0470d0,  -6.7845d0,  17.1630d0 ]
    geom.iniP( 5).pos(1:3) = [  14.8935d0,  18.0194d0,   9.2023d0 ]
    geom.iniP( 6).pos(1:3) = [  24.1794d0,  -6.7845d0,  -0.7245d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.5010d0,  11.2349d0,   9.2023d0 ]
    geom.iniP( 8).pos(1:3) = [ -19.9568d0, -15.2446d0,  -0.7245d0 ]
    geom.iniP( 9).pos(1:3) = [  20.7254d0,  11.2349d0,  -8.6852d0 ]
    geom.iniP(10).pos(1:3) = [  13.1042d0, -15.2446d0, -15.0691d0 ]
    geom.iniP(11).pos(1:3) = [ -23.4107d0,   2.7748d0,  -8.6852d0 ]
    geom.iniP(12).pos(1:3) = [  -6.5383d0, -19.0097d0, -15.0691d0 ]
    geom.iniP(13).pos(1:3) = [   9.6502d0,   2.7748d0, -23.0298d0 ]
    geom.iniP(14).pos(1:3) = [  -9.9923d0,  -0.9903d0, -23.0298d0 ]

    ! Position vectors
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 10, 13,  9 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 8, 11, 14, 12 ]
    geom.face(7).n_poi = 4; allocate(geom.face(7).poi(4)); geom.face(7).poi(1:4) = [10, 12, 14, 13 ]
    geom.face(8).n_poi = 7; allocate(geom.face(8).poi(7)); geom.face(8).poi(1:7) = [ 1,  4,  8, 12, 10,  6, 2 ]
    geom.face(9).n_poi = 7; allocate(geom.face(9).poi(7)); geom.face(9).poi(1:7) = [ 3,  5,  9, 13, 14, 11, 7 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  24.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],   5.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], -10.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  30.0d0)
end subroutine Exam_DAEDALUS2_Prism_Heptagon

! -----------------------------------------------------------------------------

! Example of octagonal prism
subroutine Exam_DAEDALUS2_Prism_Octagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Prism_Octagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  27.9793d0 ]
    geom.iniP( 2).pos(1:3) = [  18.6790d0,  -0.0000d0,  20.8312d0 ]
    geom.iniP( 3).pos(1:3) = [  -2.7355d0,  18.4776d0,  20.8312d0 ]
    geom.iniP( 4).pos(1:3) = [ -17.8778d0,  -5.4120d0,  20.8312d0 ]
    geom.iniP( 5).pos(1:3) = [  15.9435d0,  18.4776d0,  13.6831d0 ]
    geom.iniP( 6).pos(1:3) = [  27.2173d0,  -5.4120d0,   3.5741d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.6133d0,  13.0657d0,  13.6831d0 ]
    geom.iniP( 8).pos(1:3) = [ -24.4818d0, -13.0657d0,   3.5741d0 ]
    geom.iniP( 9).pos(1:3) = [  24.4818d0,  13.0657d0,  -3.5741d0 ]
    geom.iniP(10).pos(1:3) = [  20.6133d0, -13.0657d0, -13.6831d0 ]
    geom.iniP(11).pos(1:3) = [ -27.2173d0,   5.4120d0,  -3.5741d0 ]
    geom.iniP(12).pos(1:3) = [ -15.9435d0, -18.4776d0, -13.6831d0 ]
    geom.iniP(13).pos(1:3) = [  17.8778d0,   5.4120d0, -20.8312d0 ]
    geom.iniP(14).pos(1:3) = [   2.7355d0, -18.4776d0, -20.8312d0 ]
    geom.iniP(15).pos(1:3) = [ -18.6790d0,  -0.0000d0, -20.8312d0 ]
    geom.iniP(16).pos(1:3) = [   0.0000d0,  -0.0000d0, -27.9793d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 6, 10, 13,  9 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 8, 11, 15, 12 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [10, 14, 16, 13 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [12, 15, 16, 14 ]
    geom.face( 9).n_poi = 8; allocate(geom.face( 9).poi(8)); geom.face( 9).poi(1:8) = [ 1,  4,  8, 12, 14, 10,  6, 2 ]
    geom.face(10).n_poi = 8; allocate(geom.face(10).poi(8)); geom.face(10).poi(1:8) = [ 3,  5,  9, 13, 16, 15, 11, 7 ]
end subroutine Exam_DAEDALUS2_Prism_Octagon

! -----------------------------------------------------------------------------

! Example of enneagonal prism
subroutine Exam_DAEDALUS2_Prism_Enneagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Prism_Enneagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 11

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  30.9008d0 ]
    geom.iniP( 2).pos(1:3) = [  18.9238d0,  -0.0000d0,  24.4286d0 ]
    geom.iniP( 3).pos(1:3) = [  -2.2137d0,  18.7939d0,  24.4286d0 ]
    geom.iniP( 4).pos(1:3) = [ -18.4059d0,  -4.3969d0,  24.4286d0 ]
    geom.iniP( 5).pos(1:3) = [  16.7101d0,  18.7939d0,  17.9562d0 ]
    geom.iniP( 6).pos(1:3) = [  29.5109d0,  -4.3969d0,   8.0401d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.6196d0,  14.3969d0,  17.9562d0 ]
    geom.iniP( 8).pos(1:3) = [ -27.6816d0, -11.1334d0,   8.0401d0 ]
    geom.iniP( 9).pos(1:3) = [  27.2972d0,  14.3969d0,   1.5678d0 ]
    geom.iniP(10).pos(1:3) = [  26.8073d0, -11.1334d0, -10.5962d0 ]
    geom.iniP(11).pos(1:3) = [ -29.8953d0,   7.6605d0,   1.5678d0 ]
    geom.iniP(12).pos(1:3) = [ -23.4868d0, -17.0574d0, -10.5962d0 ]
    geom.iniP(13).pos(1:3) = [  24.5937d0,   7.6605d0, -17.0685d0 ]
    geom.iniP(14).pos(1:3) = [  12.0783d0, -17.0574d0, -22.7602d0 ]
    geom.iniP(15).pos(1:3) = [ -25.7005d0,   1.7365d0, -17.0685d0 ]
    geom.iniP(16).pos(1:3) = [  -7.7844d0, -19.3970d0, -22.7602d0 ]
    geom.iniP(17).pos(1:3) = [   9.8646d0,   1.7365d0, -29.2325d0 ]
    geom.iniP(18).pos(1:3) = [  -9.9981d0,  -0.6031d0, -29.2325d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 6, 10, 13,  9 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 8, 11, 15, 12 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [10, 14, 17, 13 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [12, 15, 18, 16 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [14, 16, 18, 17 ]
    geom.face(10).n_poi = 9; allocate(geom.face(10).poi(9)); geom.face(10).poi(1:9) = [ 1,  4,  8, 12, 16, 14, 10,  6, 2 ]
    geom.face(11).n_poi = 9; allocate(geom.face(11).poi(9)); geom.face(11).poi(1:9) = [ 3,  5,  9, 13, 17, 18, 15, 11, 7 ]
end subroutine Exam_DAEDALUS2_Prism_Enneagon

! -----------------------------------------------------------------------------

! Example of decagonal prism
subroutine Exam_DAEDALUS2_Prism_Decagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Prism_Decagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,   0.0000d0,  33.8707d0 ]
    geom.iniP( 2).pos(1:3) = [  19.1085d0,   0.0000d0,  27.9657d0 ]
    geom.iniP( 3).pos(1:3) = [  -1.8247d0,  19.0211d0,  27.9657d0 ]
    geom.iniP( 4).pos(1:3) = [ -18.7600d0,  -3.6327d0,  27.9657d0 ]
    geom.iniP( 5).pos(1:3) = [  17.2838d0,  19.0211d0,  22.0609d0 ]
    geom.iniP( 6).pos(1:3) = [  31.2666d0,  -3.6327d0,  12.5066d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.5847d0,  15.3884d0,  22.0609d0 ]
    geom.iniP( 8).pos(1:3) = [ -30.0058d0,  -9.5106d0,  12.5066d0 ]
    geom.iniP( 9).pos(1:3) = [  29.4419d0,  15.3884d0,   6.6018d0 ]
    geom.iniP(10).pos(1:3) = [  31.8305d0,  -9.5106d0,  -6.6018d0 ]
    geom.iniP(11).pos(1:3) = [ -31.8305d0,   9.5106d0,   6.6018d0 ]
    geom.iniP(12).pos(1:3) = [ -29.4419d0, -15.3884d0,  -6.6018d0 ]
    geom.iniP(13).pos(1:3) = [  30.0058d0,   9.5106d0, -12.5066d0 ]
    geom.iniP(14).pos(1:3) = [  20.5847d0, -15.3884d0, -22.0609d0 ]
    geom.iniP(15).pos(1:3) = [ -31.2666d0,   3.6327d0, -12.5066d0 ]
    geom.iniP(16).pos(1:3) = [ -17.2838d0, -19.0211d0, -22.0609d0 ]
    geom.iniP(17).pos(1:3) = [  18.7600d0,   3.6327d0, -27.9657d0 ]
    geom.iniP(18).pos(1:3) = [   1.8247d0, -19.0211d0, -27.9657d0 ]
    geom.iniP(19).pos(1:3) = [ -19.1085d0,   0.0000d0, -27.9657d0 ]
    geom.iniP(20).pos(1:3) = [   0.0000d0,   0.0000d0, -33.8707d0 ]

    ! Position vectors
    geom.face( 1).n_poi =  4; allocate(geom.face( 1).poi( 4)); geom.face( 1).poi(1: 4) = [ 1,  2,  5,  3 ]
    geom.face( 2).n_poi =  4; allocate(geom.face( 2).poi( 4)); geom.face( 2).poi(1: 4) = [ 1,  3,  7,  4 ]
    geom.face( 3).n_poi =  4; allocate(geom.face( 3).poi( 4)); geom.face( 3).poi(1: 4) = [ 2,  6,  9,  5 ]
    geom.face( 4).n_poi =  4; allocate(geom.face( 4).poi( 4)); geom.face( 4).poi(1: 4) = [ 4,  7, 11,  8 ]
    geom.face( 5).n_poi =  4; allocate(geom.face( 5).poi( 4)); geom.face( 5).poi(1: 4) = [ 6, 10, 13,  9 ]
    geom.face( 6).n_poi =  4; allocate(geom.face( 6).poi( 4)); geom.face( 6).poi(1: 4) = [ 8, 11, 15, 12 ]
    geom.face( 7).n_poi =  4; allocate(geom.face( 7).poi( 4)); geom.face( 7).poi(1: 4) = [10, 14, 17, 13 ]
    geom.face( 8).n_poi =  4; allocate(geom.face( 8).poi( 4)); geom.face( 8).poi(1: 4) = [12, 15, 19, 16 ]
    geom.face( 9).n_poi =  4; allocate(geom.face( 9).poi( 4)); geom.face( 9).poi(1: 4) = [14, 18, 20, 17 ]
    geom.face(10).n_poi =  4; allocate(geom.face(10).poi( 4)); geom.face(10).poi(1: 4) = [16, 19, 20, 18 ]
    geom.face(11).n_poi = 10; allocate(geom.face(11).poi(10)); geom.face(11).poi(1:10) = [ 1,  4,  8, 12, 16, 18, 14, 10,  6, 2 ]
    geom.face(12).n_poi = 10; allocate(geom.face(12).poi(10)); geom.face(12).poi(1:10) = [ 3,  5,  9, 13, 17, 20, 19, 15, 11, 7 ]
end subroutine Exam_DAEDALUS2_Prism_Decagon

! -----------------------------------------------------------------------------

! Example of square antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Square(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Antiprism_Square"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 8
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  16.4533d0 ]
    geom.iniP( 2).pos(1:3) = [  15.8821d0,  -0.0000d0,   4.2977d0 ]
    geom.iniP( 3).pos(1:3) = [   3.2893d0,  15.5378d0,   4.2977d0 ]
    geom.iniP( 4).pos(1:3) = [ -14.5196d0,   6.4360d0,   4.2977d0 ]
    geom.iniP( 5).pos(1:3) = [  -9.3035d0, -12.8719d0,   4.2977d0 ]
    geom.iniP( 6).pos(1:3) = [   6.5786d0, -12.8719d0,  -7.8580d0 ]
    geom.iniP( 7).pos(1:3) = [   7.9411d0,   6.4360d0, -12.8930d0 ]
    geom.iniP( 8).pos(1:3) = [  -9.8679d0,  -2.6659d0, -12.8930d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1, 2, 3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1, 3, 4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1, 4, 5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2, 6, 7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2, 7, 3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4, 8, 5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5, 8, 6 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 8, 7 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 1, 5, 6, 2 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 3, 7, 8, 4 ]
end subroutine Exam_DAEDALUS2_Anti_Prism_Square

! -----------------------------------------------------------------------------

! Example of hexagonal antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Hexagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Antiprism_Hexagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  21.7533d0 ]
    geom.iniP( 2).pos(1:3) = [  17.7615d0,  -0.0000d0,  12.5593d0 ]
    geom.iniP( 3).pos(1:3) = [   6.5012d0,  16.5290d0,  12.5593d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.0023d0,  12.1000d0,  12.5593d0 ]
    geom.iniP( 5).pos(1:3) = [ -16.0195d0,  -7.6711d0,  12.5593d0 ]
    geom.iniP( 6).pos(1:3) = [  19.5035d0,  -7.6711d0,  -5.8288d0 ]
    geom.iniP( 7).pos(1:3) = [  17.7615d0,  12.1000d0,  -3.3652d0 ]
    geom.iniP( 8).pos(1:3) = [ -21.2456d0,   3.2422d0,  -3.3652d0 ]
    geom.iniP( 9).pos(1:3) = [ -14.2775d0, -15.3422d0,  -5.8288d0 ]
    geom.iniP(10).pos(1:3) = [   3.4840d0, -15.3422d0, -15.0228d0 ]
    geom.iniP(11).pos(1:3) = [   9.5184d0,   3.2422d0, -19.2898d0 ]
    geom.iniP(12).pos(1:3) = [  -9.9851d0,  -1.1867d0, -19.2898d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 10 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 12, 11 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [ 1,  5,  9, 10, 6, 2 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [ 3,  7, 11, 12, 8, 4 ]

    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  24.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 1.0d0, 0.0d0],   5.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], -11.0d0)
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0],  30.0d0)
end subroutine Exam_DAEDALUS2_Anti_Prism_Hexagon

! -----------------------------------------------------------------------------

! Example of heptagonal antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Heptagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Antiprism Heptagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,  -0.0000d0,  24.5946d0 ]
    geom.iniP( 2).pos(1:3) = [  18.2722d0,  -0.0000d0,  16.4627d0 ]
    geom.iniP( 3).pos(1:3) = [   7.3266d0,  16.7390d0,  16.4627d0 ]
    geom.iniP( 4).pos(1:3) = [ -12.3967d0,  13.4237d0,  16.4627d0 ]
    geom.iniP( 5).pos(1:3) = [ -17.2680d0,  -5.9741d0,  16.4627d0 ]
    geom.iniP( 6).pos(1:3) = [  23.7892d0,  -5.9741d0,  -1.8095d0 ]
    geom.iniP( 7).pos(1:3) = [  20.5286d0,  13.4237d0,   1.8095d0 ]
    geom.iniP( 8).pos(1:3) = [ -23.7892d0,   5.9741d0,   1.8095d0 ]
    geom.iniP( 9).pos(1:3) = [ -20.5286d0, -13.4237d0,  -1.8095d0 ]
    geom.iniP(10).pos(1:3) = [  12.3967d0, -13.4237d0, -16.4627d0 ]
    geom.iniP(11).pos(1:3) = [  17.2680d0,   5.9741d0, -16.4627d0 ]
    geom.iniP(12).pos(1:3) = [ -18.2722d0,  -0.0000d0, -16.4627d0 ]
    geom.iniP(13).pos(1:3) = [  -7.3266d0, -16.7390d0, -16.4627d0 ]
    geom.iniP(14).pos(1:3) = [  -0.0000d0,  -0.0000d0, -24.5946d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 13, 14 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [10, 14, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [12, 14, 13 ]
    geom.face(15).n_poi = 7; allocate(geom.face(15).poi(7)); geom.face(15).poi(1:7) = [ 1,  5,  9, 13, 10, 6, 2 ]
    geom.face(16).n_poi = 7; allocate(geom.face(16).poi(7)); geom.face(16).poi(1:7) = [ 3,  7, 11, 14, 12, 8, 4 ]
end subroutine Exam_DAEDALUS2_Anti_Prism_Heptagon

! -----------------------------------------------------------------------------

! Example of octagonal antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Octagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Antiprism_Octagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,   0.0000d0,  27.5111d0 ]
    geom.iniP( 2).pos(1:3) = [  18.6320d0,   0.0000d0,  20.2411d0 ]
    geom.iniP( 3).pos(1:3) = [   7.8977d0,  16.8753d0,  20.2411d0 ]
    geom.iniP( 4).pos(1:3) = [ -11.9366d0,  14.3062d0,  20.2411d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.0171d0,  -4.7471d0,  20.2411d0 ]
    geom.iniP( 6).pos(1:3) = [  26.9645d0,  -4.7471d0,   2.6902d0 ]
    geom.iniP( 7).pos(1:3) = [  22.4908d0,  14.3062d0,   6.8083d0 ]
    geom.iniP( 8).pos(1:3) = [ -25.3935d0,   8.1038d0,   6.8083d0 ]
    geom.iniP( 9).pos(1:3) = [ -24.8651d0, -11.4605d0,   2.6902d0 ]
    geom.iniP(10).pos(1:3) = [  20.1164d0, -11.4605d0, -14.8607d0 ]
    geom.iniP(11).pos(1:3) = [  23.2942d0,   8.1038d0, -12.1887d0 ]
    geom.iniP(12).pos(1:3) = [ -24.5901d0,   1.9014d0, -12.1887d0 ]
    geom.iniP(13).pos(1:3) = [ -16.5326d0, -16.2076d0, -14.8607d0 ]
    geom.iniP(14).pos(1:3) = [   2.0994d0, -16.2076d0, -22.1306d0 ]
    geom.iniP(15).pos(1:3) = [   9.8372d0,   1.9014d0, -25.6216d0 ]
    geom.iniP(16).pos(1:3) = [  -9.9971d0,  -0.6677d0, -25.6216d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 14, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [10, 15, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [12, 16, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [13, 16, 14 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [14, 16, 15 ]
    geom.face(17).n_poi = 8; allocate(geom.face(17).poi(8)); geom.face(17).poi(1:8) = [ 1,  5,  9, 13, 14, 10, 6, 2 ]
    geom.face(18).n_poi = 8; allocate(geom.face(18).poi(8)); geom.face(18).poi(1:8) = [ 3,  7, 11, 15, 16, 12, 8, 4 ]

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [1.0d0, 0.0d0, 0.0d0], 14.0d0)
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], -7.0d0)
end subroutine Exam_DAEDALUS2_Anti_Prism_Octagon

! -----------------------------------------------------------------------------

! Example of enneagonal antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Enneagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Antiprism_Enneagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  30.4809d0 ]
    geom.iniP( 2).pos(1:3) = [  18.8930d0,  -0.0000d0,  23.9195d0 ]
    geom.iniP( 3).pos(1:3) = [   8.3072d0,  16.9688d0,  23.9195d0 ]
    geom.iniP( 4).pos(1:3) = [ -11.5879d0,  14.9221d0,  23.9195d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.4974d0,  -3.8465d0,  23.9195d0 ]
    geom.iniP( 6).pos(1:3) = [  29.3415d0,  -3.8465d0,   7.3052d0 ]
    geom.iniP( 7).pos(1:3) = [  23.9195d0,  14.9221d0,  11.5879d0 ]
    geom.iniP( 8).pos(1:3) = [ -26.4565d0,   9.7397d0,  11.5879d0 ]
    geom.iniP( 9).pos(1:3) = [ -27.9439d0,  -9.7397d0,   7.3052d0 ]
    geom.iniP(10).pos(1:3) = [  26.4565d0,  -9.7397d0, -11.5879d0 ]
    geom.iniP(11).pos(1:3) = [  27.9439d0,   9.7397d0,  -7.3052d0 ]
    geom.iniP(12).pos(1:3) = [ -29.3415d0,   3.8465d0,  -7.3052d0 ]
    geom.iniP(13).pos(1:3) = [ -23.9195d0, -14.9221d0, -11.5879d0 ]
    geom.iniP(14).pos(1:3) = [  11.5879d0, -14.9221d0, -23.9195d0 ]
    geom.iniP(15).pos(1:3) = [  18.4974d0,   3.8465d0, -23.9195d0 ]
    geom.iniP(16).pos(1:3) = [ -18.8930d0,  -0.0000d0, -23.9195d0 ]
    geom.iniP(17).pos(1:3) = [  -8.3072d0, -16.9688d0, -23.9195d0 ]
    geom.iniP(18).pos(1:3) = [   0.0000d0,  -0.0000d0, -30.4809d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 14, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [10, 15, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [12, 16, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [13, 16, 17 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [14, 17, 18 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [14, 18, 15 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [16, 18, 17 ]
    geom.face(19).n_poi = 9; allocate(geom.face(19).poi(9)); geom.face(19).poi(1:9) = [ 1,  5,  9, 13, 17, 14, 10, 6, 2 ]
    geom.face(20).n_poi = 9; allocate(geom.face(20).poi(9)); geom.face(20).poi(1:9) = [ 3,  7, 11, 15, 18, 16, 12, 8, 4 ]
end subroutine Exam_DAEDALUS2_Anti_Prism_Enneagon

! -----------------------------------------------------------------------------

! Example of decagonal antiprism
subroutine Exam_DAEDALUS2_Anti_Prism_Decagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "Antiprism_Decagon"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 22

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  33.4901d0 ]
    geom.iniP( 2).pos(1:3) = [  19.0876d0,  -0.0000d0,  27.5182d0 ]
    geom.iniP( 3).pos(1:3) = [   8.6096d0,  17.0356d0,  27.5182d0 ]
    geom.iniP( 4).pos(1:3) = [ -11.3208d0,  15.3681d0,  27.5182d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.8222d0,  -3.1719d0,  27.5182d0 ]
    geom.iniP( 6).pos(1:3) = [  31.1498d0,  -3.1719d0,  11.8835d0 ]
    geom.iniP( 7).pos(1:3) = [  24.9860d0,  15.3681d0,  16.1589d0 ]
    geom.iniP( 8).pos(1:3) = [ -27.1924d0,  11.0023d0,  16.1589d0 ]
    geom.iniP( 9).pos(1:3) = [ -30.1896d0,  -8.3041d0,  11.8835d0 ]
    geom.iniP(10).pos(1:3) = [  31.5792d0,  -8.3041d0,  -7.4420d0 ]
    geom.iniP(11).pos(1:3) = [  31.5533d0,  11.0023d0,  -2.2207d0 ]
    geom.iniP(12).pos(1:3) = [ -32.9429d0,   5.6060d0,  -2.2207d0 ]
    geom.iniP(13).pos(1:3) = [ -29.7602d0, -13.4364d0,  -7.4420d0 ]
    geom.iniP(14).pos(1:3) = [  20.2118d0, -13.4364d0, -23.0767d0 ]
    geom.iniP(15).pos(1:3) = [  25.8028d0,   5.6060d0, -20.6004d0 ]
    geom.iniP(16).pos(1:3) = [ -26.3756d0,   1.2402d0, -20.6004d0 ]
    geom.iniP(17).pos(1:3) = [ -17.6980d0, -16.6083d0, -23.0767d0 ]
    geom.iniP(18).pos(1:3) = [   1.3896d0, -16.6083d0, -29.0487d0 ]
    geom.iniP(19).pos(1:3) = [   9.9312d0,   1.2402d0, -31.9597d0 ]
    geom.iniP(20).pos(1:3) = [  -9.9992d0,  -0.4274d0, -31.9597d0 ]

    ! Position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi( 3)); geom.face( 6).poi(1: 3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi( 3)); geom.face( 7).poi(1: 3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi( 3)); geom.face( 8).poi(1: 3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi( 3)); geom.face( 9).poi(1: 3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi( 3)); geom.face(10).poi(1: 3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi( 3)); geom.face(11).poi(1: 3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi( 3)); geom.face(12).poi(1: 3) = [10, 14, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi( 3)); geom.face(13).poi(1: 3) = [10, 15, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi( 3)); geom.face(14).poi(1: 3) = [12, 16, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi( 3)); geom.face(15).poi(1: 3) = [13, 16, 17 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi( 3)); geom.face(16).poi(1: 3) = [14, 18, 19 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi( 3)); geom.face(17).poi(1: 3) = [14, 19, 15 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi( 3)); geom.face(18).poi(1: 3) = [16, 20, 17 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi( 3)); geom.face(19).poi(1: 3) = [17, 20, 18 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi( 3)); geom.face(20).poi(1: 3) = [18, 20, 19 ]
    geom.face(21).n_poi =10; allocate(geom.face(21).poi(10)); geom.face(21).poi(1:10) = [ 1,  5,  9, 13, 17, 18, 14, 10, 6, 2 ]
    geom.face(22).n_poi =10; allocate(geom.face(22).poi(10)); geom.face(22).poi(1:10) = [ 3,  7, 11, 15, 19, 20, 16, 12, 8, 4 ]
end subroutine Exam_DAEDALUS2_Anti_Prism_Decagon

! -----------------------------------------------------------------------------

! Example of open end triangular prism with quad mesh
subroutine Exam_DAEDALUS2_End_Tri_Prism_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Tri_Prism_Quad"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 2
    nz = n
    nr = 3

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Tri_Prism_Quad

! -----------------------------------------------------------------------------

! Example of open end cube with quad mesh
subroutine Exam_DAEDALUS2_End_Cube_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Cube_Quad"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 3
    nz = n
    nr = 4

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Cube_Quad

! -----------------------------------------------------------------------------

! Example of open end pentagonal prism with quad mesh
subroutine Exam_DAEDALUS2_End_Penta_Prism_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Penta_Prism_Quad"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 3
    nz = n
    nr = 5

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Penta_Prism_Quad

! -----------------------------------------------------------------------------

! Example of open end cylinder with quad mesh
subroutine Exam_DAEDALUS2_End_Cylinder_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Cylinder_Quad"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 4
    nz = n
    nr = 8

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)  
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Cylinder_Quad

! -----------------------------------------------------------------------------

! Example of open end triangular prism with tri mesh
subroutine Exam_DAEDALUS2_End_Tri_Prism_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Tri_Prism_Tri"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 3
    nz = n
    nr = 3

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Tri_Prism_Tri

! -----------------------------------------------------------------------------

! Example of open end cube with tri mesh
subroutine Exam_DAEDALUS2_End_Cube_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Cube_Tri"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 3
    nz = n
    nr = 4

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Cube_Tri

! -----------------------------------------------------------------------------

! Example of open end pentagonal prism with tri mesh
subroutine Exam_DAEDALUS2_End_Penta_Prism_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face

    ! Set problem
    prob.name_prob = "End_Penta_Prism_Tri"
    call Mani_Set_Prob(prob, 'red')

    ! Set mesh
    n  = 3
    nz = n
    nr = 5

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    ! The number of points and faces
    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do

    ! Set the orientation of the geometry
    call Mani_Set_3DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 0.0d0)
end subroutine Exam_DAEDALUS2_End_Penta_Prism_Tri

! -----------------------------------------------------------------------------

subroutine Prob_Elem_Mesh(con, element, n_con, out_con, nNPE)
    integer, intent(in)    :: con(:)
    integer, intent(in)    :: element(:,:), nNPE
    integer, intent(inout) :: n_con
    integer, intent(out)   :: out_con(:,:)

    integer :: i, j, k

    j = 0
    do i = 1, n_con
        j = j + 1
        do k = 1, nNPE
            out_con(j,k) = element(i, con(k))
        end do
        j = j + 1
        do k = 1, nNPE
            out_con(j,k) = element(i, con(nNPE+k))
        end do
    end do
    n_con = j
end subroutine Prob_Elem_Mesh

! -----------------------------------------------------------------------------

end module Exam_DAEDALUS2