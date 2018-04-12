    subroutine parameters
    implicit none

    pi =4.*atan(1.)
    one_third = 1./3.
    two_third = 2./3.
    over12 = 1./12.
    over6 = 1./6.
    ! compatational domain

    xleft = -pi
    xright = pi
    yleft = -pi
    yright = pi

    xlength = xright - xleft
    ylength = yright - yleft

    tprint = 2*1.5
    !tprint =0.01
    !tprint = pi
    cfl = 1./2
    ! cfl = 1.
    irk = 5
    ighost = int(cfl) + 7
    ! if isearch = 1, call search_quadrilateral.
    ! if isearch = 2, call search_QC.
    isearch = 2
    n_moment = 6
    limiter = 1


    para_dx3 = 0.01
    para_limit = 1d-9
    !*******************************
    tvb_m = 1.
    r0_dgweno = 0.996
    r1_dgweno = 0.001
    r2_dgweno = 0.001
    r3_dgweno = 0.001
    r4_dgweno = 0.001
    eps = 1d-6
    !********************************************************please set
    !

    !
    xg(1)=-0.466234757101576013906150777246997304567d0
    xg(2)=-0.330604693233132256830699797509952673503d0
    xg(3)=-0.119309593041598454315250860840355967709d0
    xg(4)=-xg(3)
    xg(5)=-xg(2)
    xg(6)=-xg(1)
    wg(1)=1.71324492379170345040296142172733d-1/2d0
    wg(2)=3.60761573048138607569833513837716d-1/2d0
    wg(3)=4.67913934572691047389870343989551d-1/2d0
    wg(4)=wg(3)
    wg(5)=wg(2)
    wg(6)=wg(1)


    gauss2(1,1)=-sqrt(3.0)/6.0
    gauss2(2,1)=sqrt(3.0)/6.0
    gauss2(1,2)=0.5
    gauss2(2,2)=0.5

    gauss3(1,1) = -sqrt(15.)/10.
    gauss3(2,1) = 0.
    gauss3(3,1) = sqrt(15.)/10.
    gauss3(1,2) = 5./18.
    gauss3(2,2) =  4./9.
    gauss3(3,2) = 5./18.
    !
    ai(1) = 1.
    ai(2) = 12.
    ai(3) = 12.
    ai(4) = 180.
    ai(5) = 144.
    ai(6) = 180.

    end subroutine parameters