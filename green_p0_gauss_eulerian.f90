    subroutine green_p0_gauss_eulerian(sum16)
    ! The subroutine is the factorization version of
    ! subroutine green_p1_gauss17
    !
    implicit none
    integer :: kk
    real :: x1,y1,x2,y2
    real :: sum16(1)
    real :: sum

    real :: x11,y11,x22,y22
    integer :: idx,idy

    real :: slope
    real :: x2m1,y2m1

    !real :: x_base,y_base
    !
    !x_base = 0.25*( Dij_star(i,j)%vertex_star(1,1) + Dij_star(i,j)%vertex_star(2,1)  &
    !    + Dij_star(i,j)%vertex_star(3,1)  + Dij_star(i,j)%vertex_star(4,1)  )
    !y_base = 0.25*( Dij_star(i,j)%vertex_star(1,2) + Dij_star(i,j)%vertex_star(2,2)  &
    !    + Dij_star(i,j)%vertex_star(3,2)  + Dij_star(i,j)%vertex_star(4,2)  )

    sum = 0.
    !

    do kk = 1,Dij_star(i,j)%nsub_outer
        x11 = Dij_star(i,j)%segment_outer(kk)%vl(1)
        y11 = Dij_star(i,j)%segment_outer(kk)%vl(2)
        x22 = Dij_star(i,j)%segment_outer(kk)%vr(1)
        y22 = Dij_star(i,j)%segment_outer(kk)%vr(2)
        !
        idx = Dij_star(i,j)%segment_outer(kk)%id(1)
        idy = Dij_star(i,j)%segment_outer(kk)%id(2)
        x1 = (x11 -x(idx) )*overdx ;
        y1 = (y11 -y(idy) )*overdy ;
        x2 = (x22 -x(idx) )*overdx ;
        y2 = (y22 -y(idy) )*overdy ;


        x2m1 = x2-x1
        y2m1 = y2-y1
        !*******************
        ! c2_10 = c1_00
        ! c2_20 = c1_10
        ! c2_11 = c1_01
        !*******************
        ! c3_01 = c1_00
        ! c3_11 = c1_10
        ! c3_02 = c1_01
        !*******************
        if( abs(x2m1)< 1d-12 .and. abs(y2m1)<1d-12 )then

        else
            if( abs(x2m1) > abs(y2m1) )then
                slope = (y2m1)/(x2m1)
                sum = sum + Dij(idx,idy)%umodal(1)* ( x1+x2 )*0.5*slope *x2m1

            else
                sum = sum + Dij(idx,idy)%umodal(1) *( x1+x2 )*0.5*y2m1

            endif
        endif

    enddo

    do kk = 1,Dij_star(i,j)%nsub_inner
        x11 = Dij_star(i,j)%segment_inner(kk)%vl(1)
        y11 = Dij_star(i,j)%segment_inner(kk)%vl(2)
        x22 = Dij_star(i,j)%segment_inner(kk)%vr(1)
        y22 = Dij_star(i,j)%segment_inner(kk)%vr(2)
        !
        idx = Dij_star(i,j)%segment_inner(kk)%id(1)
        idy = Dij_star(i,j)%segment_inner(kk)%id(2)
        x1 = (x11 -x(idx) )*overdx ;
        y1 = (y11 -y(idy) )*overdy ;
        x2 = (x22 -x(idx) )*overdx ;
        y2 = (y22 -y(idy) )*overdy ;

        x2m1 = x2-x1
        y2m1 = y2-y1
        !*******************
        ! c2_10 = c1_00
        ! c2_20 = c1_10
        ! c2_11 = c1_01
        !*******************
        ! c3_01 = c1_00
        ! c3_11 = c1_10
        ! c3_02 = c1_01
        !*******************
        if( abs(x2m1)< 1d-12 .and. abs(y2m1)<1d-12 )then

        else
            if( abs(x2m1) > abs(y2m1) )then
                slope = (y2m1)/(x2m1)
                sum = sum + Dij(idx,idy)%umodal(1)* ( x1+x2 )*0.5*slope *x2m1
            else
                sum = sum + Dij(idx,idy)%umodal(1) *( x1+x2 )*0.5*y2m1
            endif
        endif
    enddo
    sum16(1) = sum

    end subroutine green_p0_gauss_eulerian
