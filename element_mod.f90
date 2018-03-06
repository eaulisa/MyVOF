    module element_mod

    type, public ::element
        sequence
        real umodal(6)
        real vertex(6,2)
        real node(9,2)
        real temp(5)
    end type

    type(element),allocatable, public :: Dij(:,:)
    ! allocatable DG solution space

    !---------------------------------------------------------------
    type, public :: segment1
        sequence
        real vl(2)
        real vr(2)
        integer :: id(2)

        real v11(2)
        real v22(2)
        real v33(2)
        real :: c_ab(6)

        integer :: icurve
    end type


    type, public ::element_star
        sequence

        real vertex_star(6,2)
        real node_star(9,2)
        integer :: nsub_outer
        integer :: nsub_inner
        type(segment1) segment_outer( 35 )
        type(segment1) segment_inner( 35 )

    end type

    type(element_star),allocatable, public :: Dij_star(:,:)

    end module element_mod