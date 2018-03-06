    !**************************************************************************
    ! code 2D SLDG pk
    ! accuracy test code for
    !                Swirling deformation flow.
    !
    ! u_t -( cos^2(x/2)*sin(y)*g(t)*u )_x + ( sin(x)*cos^2(y/2)*g(t)*u )_y = 0,
    !                              x\in [ -pi,pi ]\times[ -pi,pi ],
    !   where g(t) = cos( pi*t/T )\pi. The initial condition is set to be
    !   the following smooth cosine bells (with C^5 smoothness).
    !
    !**************************************************************************
    !  set parameters at subroutine parameters!
    !     isearch for QC or not
    !     n_moment for Pk
    !     limiter for WENO limiter
    !**************************************************************************
    !
    !                    code by
    !                            Jingmei Qiu (jingqiu@math.uh.edu)
    !                            Wei Guo (wguo@math.msu.edu)
    !                            Xiaofeng Cai (xfcai@math.uh.edu
    !                                      & xfcai@stu.xmu.edu.cn)
    !                            May 3, 2017
    !                            University of Houston
    !**************************************************************************
    program SLDG2D_transport

    use globals2d
    use element_mod
    use LU
    implicit none

    call parameters

    do kkkk = 1,5
        nx = 10*2**(kkkk-1)
        ny = nx

        ! Grid generator
        dx = (xright - xleft)/nx
        dy = (yright - yleft)/ny
        overdx = 1./dx
        overdy = 1./dy

        call allocate_variable

        call init
        if(limiter == 1)then
            call boundary
            call trouble_tvb
            if( n_moment == 3 )then
                call wenolimiter1
            elseif(n_moment == 6)then
                call wenolimiter2
            endif
        endif

        tnum = 0.

        nt =0

        !******************** BEGIN TIME EVOLUTION ***************************
        do while(tnum<tprint-1.d-11)
            call setdt
            call boundary
            call scheme_SLDG
            if(limiter == 1)then
                call boundary
                call trouble_tvb
                if( n_moment == 3 )then
                    call wenolimiter1
                elseif(n_moment == 6)then
                    call wenolimiter2
                endif
            endif
            nt = nt + 1
            if(nt/2*2==nt) print *,tnum,tnum/tprint*100,"%"
        enddo
        call order_DG
        call deallocate_variable
    enddo


    contains
    include "parameters.f90"
    include "init.f90"
    include "setdt.f90"
    include "boundary.f90"

    include "RK.f90"

    include "order.f90"
    include "get_matrix.f90"

    include "get_matrix2.f90"


    include "curve_get_line.f90"
    include "curve_coordinate_trans.f90"

    include "id_get.f90"

    include "curve_integral_function.f90"


    include "scheme_SLDG.f90"

    include "trouble_tvb.f90"
    include "weno_limiter1.f90"
    include "weno_limiter2.f90"

    include "allocate_variable.f90"

    include "get_matrix_vector.f90"


    include "green_p0_gauss_eulerian.f90"
    include "green_p1_gauss_eulerian.f90"

    include "green_p2_gauss3_eulerian.f90"

    include "green_p2QC_gauss3_eulerian.f90"


    include "search_quadrilateral_final.f90"

    include "search_QC_formula.f90"
    include "get_quadratic_root.f90"
    end program SLDG2D_transport