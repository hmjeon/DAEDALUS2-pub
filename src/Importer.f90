!
! =============================================================================
!
! Module - Importer
! Last Updated : 01/09/2019, by Hyungmin Jun (hyungminjun@outlook.com)
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
module Importer

    use Ifport

    use Para
    use Data_Prob
    use Data_Geom

    use Math

    implicit none

    public Importer_PLY
    public Importer_STL
    public Importer_WRL

contains

! -----------------------------------------------------------------------------

! Import PLY file
subroutine Importer_PLY(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(200) :: ctemp, path
    integer :: i, status, npoint, ioerr
    logical :: here

    ! Mesh data structure
    type :: MeshType
        integer :: cn(100)
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable, dimension (:) :: Basepair_con

    path = trim(prob.path_input)//trim(prob.name_file)//"."//trim(prob.type_file)

    inquire(file = trim(path), exist=here)
    if(here == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | No ply file.                                     |"
        write(p_redir, "(a)"), " +==================================================+"
        write(p_redir, "(a)")
        stop
    end if
    open(unit = 1001, file = trim(path), form = "formatted")

    do
        read(1001, "(a100)", iostat = status), ctemp

        ! Negative value, means the end-of-file (EOF) mark was read
        if(status < 0) exit

        if(index(ctemp, "format")) then

            ! Read the number of points
            read(1001, *, iostat = ioerr), ctemp, ctemp, geom.n_iniP
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | The format for # is not correct in GEO file.     |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if
            allocate(geom.iniP(geom.n_iniP))

        else if(index(ctemp, "property float32 z") .or. index(ctemp, "property float z")) then

            ! Read the number of faces
            read(1001, *, iostat = ioerr), ctemp, ctemp, geom.n_face
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | The format for faces is not correct in PLY file. |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

        else if(index(ctemp,"end_header")) then

            ! Read point
            do i = 1, geom.n_iniP
                read(1001, *, iostat = ioerr), geom.iniP(i).pos(1:3)
                if(ioerr /= 0) then
                    write(p_redir, "(a)")
                    write(p_redir, "(a)"), " +=== error ========================================+"
                    write(p_redir, "(a)"), " | The format for points is not correct in PLY file.|"
                    write(p_redir, "(a)"), " +==================================================+"
                    stop
                end if
            end do

            allocate(geom.face(geom.n_face))
            allocate(Basepair_con(geom.n_face))

            ! Read face connectivity
            do i = 1, geom.n_face
                read(1001, *, iostat = ioerr), npoint, Basepair_con(i).cn(1:npoint)
                if(ioerr /= 0) then
                    write(p_redir, "(a)")
                    write(p_redir, "(a)"), " +=== error ========================================+"
                    write(p_redir, "(a)"), " | The format for conns is not correct in PLY file. |"
                    write(p_redir, "(a)"), " +==================================================+"
                    stop
                end if
                geom.face(i).n_poi = npoint
                allocate(geom.face(i).poi(npoint))

                geom.face(i).poi(1:npoint) = Basepair_con(i).cn(1:npoint)
            end do

        end if
    end do

    deallocate(Basepair_con)
    close(unit = 1001)
end subroutine Importer_PLY

! -----------------------------------------------------------------------------

! Import STL format using meshconv
subroutine Importer_STL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! Run meshconv to generate *.PLY fileformat
    results = systemqq(trim("Library/meshconv -c ply Input/")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    if(results == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Failed to load the meshconv.                     |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_STL

! -----------------------------------------------------------------------------

! Import WRL format using meshconv
subroutine Importer_WRL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! Run meshconv to generate *.PLY fileformat
    results = systemqq(trim("Library/meshconv -c ply Input/")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    if(results == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Failed to load the meshconv.                     |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_WRL

! -----------------------------------------------------------------------------

end module Importer