    subroutine green_p2QC_gauss3_eulerian(sum16)
    implicit none
    !*******************************************************************************
    !
    !   Purpose   : line integral
    !
    !*******************************************************************************

    integer :: kk
    real :: x1,y1,x2,y2
    real :: sum16(6)
    real :: sum

    real :: x11,y11,x22,y22
    integer :: idx,idy
    !*********** least square ***********
    real :: A_temp(6,6)
    real :: temp_L(6,6),temp_U(6,6)
    real :: vert_temp(1:9,1:5)

    real :: b(6),aa(6),psi(9)
    !real,allocatable :: aaa(:,:,:)
    !integer,allocatable :: iswitch(:,:)
    !************************************
    integer :: nm

    real :: xxi1,eta1,xxi2,eta2,xxi3,eta3
    real :: e_y2,e_x2
    real :: a_trans,b_trans,c_trans,d_trans

    real :: ac1,bc1,cc1
    real :: ac2,bc2,cc2
    real :: exx1,exx2


    real :: distance,x33,y33,a16,b16,c16
    !
    real :: Ac ,para_dx3
    !****************************************************************17
    real :: uhat(6)

    real :: xshift, yshift
    real :: xgs(3),ygs(3)
    real :: xcc,ycc
    real :: x2m1,y2m1
    real :: slope
    real :: green_00,green_10,green_01,green_20,green_11,green_02,green_30,green_21,green_12,green_03
    real :: green_40,green_31,green_22,green_13,green_04
    real :: c1_00,c1_10,c1_01,c1_20,c1_11,c1_02
    real :: slope_x21,slope_y21
    real :: xg1_2,xg2_2,yg1_2,yg2_2
    real :: xy1,xy2
    real :: xg1_4,xg2_4,yg1_4,yg2_4
    !*************************************************QC
    real :: xig(3)
    real :: xi_gauss1,xi_gauss2,xi_gauss3, yi_gauss1,yi_gauss2,yi_gauss3
    real :: dxi_gauss1,dxi_gauss2,dxi_gauss3, dyi_gauss1,dyi_gauss2, dyi_gauss3
    real :: dx_xi

    real :: x_base,y_base
    real :: temp2,temp3
    
    real :: xg3_2,yg3_2,xy3,xg3_4,yg3_4
    ! P = 0, Q = x
    x_base = one_third*one_third*( Dij_star(i,j)%node_star(1,1) + Dij_star(i,j)%node_star(2,1)  &
        + Dij_star(i,j)%node_star(3,1)  + Dij_star(i,j)%node_star(4,1) +Dij_star(i,j)%node_star(5,1) + Dij_star(i,j)%node_star(6,1)  &
        + Dij_star(i,j)%node_star(7,1)  + Dij_star(i,j)%node_star(8,1) + Dij_star(i,j)%node_star(9,1) )
    y_base = one_third*one_third*( Dij_star(i,j)%node_star(1,2) + Dij_star(i,j)%node_star(2,2)  &
        + Dij_star(i,j)%node_star(3,2)  + Dij_star(i,j)%node_star(4,2) +Dij_star(i,j)%node_star(5,2) + Dij_star(i,j)%node_star(6,2)  &
        + Dij_star(i,j)%node_star(7,2)  + Dij_star(i,j)%node_star(8,2) + Dij_star(i,j)%node_star(9,2) )

    do nm = 1 ,6
        sum = 0.
        if(nm == 1)then
            aa(1) = 1.
            aa(2:6) = 0.
        else
            if(nm==2)then
                vert_temp(1,1) = ( Dij_star(i,j)%node_star(1,1) - x_base )*overdx
                vert_temp(2,1) = ( Dij_star(i,j)%node_star(2,1) - x_base )*overdx
                vert_temp(3,1) = ( Dij_star(i,j)%node_star(3,1) - x_base )*overdx
                vert_temp(4,1) = ( Dij_star(i,j)%node_star(4,1) - x_base )*overdx
                vert_temp(5,1) = ( Dij_star(i,j)%node_star(5,1) - x_base )*overdx
                vert_temp(6,1) = ( Dij_star(i,j)%node_star(6,1) - x_base )*overdx
                vert_temp(7,1) = ( Dij_star(i,j)%node_star(7,1) - x_base )*overdx
                vert_temp(8,1) = ( Dij_star(i,j)%node_star(8,1) - x_base )*overdx
                vert_temp(9,1) = ( Dij_star(i,j)%node_star(9,1) - x_base )*overdx

                vert_temp(1,2) = ( Dij_star(i,j)%node_star(1,2) - y_base )*overdy
                vert_temp(2,2) = ( Dij_star(i,j)%node_star(2,2) - y_base )*overdy
                vert_temp(3,2) = ( Dij_star(i,j)%node_star(3,2) - y_base )*overdy
                vert_temp(4,2) = ( Dij_star(i,j)%node_star(4,2) - y_base )*overdy
                vert_temp(5,2) = ( Dij_star(i,j)%node_star(5,2) - y_base )*overdy
                vert_temp(6,2) = ( Dij_star(i,j)%node_star(6,2) - y_base )*overdy
                vert_temp(7,2) = ( Dij_star(i,j)%node_star(7,2) - y_base )*overdy
                vert_temp(8,2) = ( Dij_star(i,j)%node_star(8,2) - y_base )*overdy
                vert_temp(9,2) = ( Dij_star(i,j)%node_star(9,2) - y_base )*overdy

                vert_temp(1:9,3) = vert_temp(1:9,1)*vert_temp(1:9,1) - over12
                vert_temp(1:9,4) = vert_temp(1:9,1)*vert_temp(1:9,2)
                vert_temp(1:9,5) = vert_temp(1:9,2)*vert_temp(1:9,2) - over12

                call get_matrix2_a( vert_temp(1:9,1:5),  A_temp(1:6,1:6) )
                call doolittle(A_temp,temp_L,temp_U,6)

                storage_L(i,j,1:6,1:6) = temp_L(1:6,1:6)
                storage_U(i,j,1:6,1:6) = temp_U(1:6,1:6)
                vert(i,j,1:9,1:5) = vert_temp(1:9,1:5)
            endif
            psi(1) = fphi( nm, (Dij(i,j)%node(1,1)-x(i) )*overdx,(Dij(i,j)%node(1,2)-y(j))*overdy ) ;
            psi(2) = fphi( nm, (Dij(i,j)%node(2,1)-x(i) )*overdx,(Dij(i,j)%node(2,2)-y(j))*overdy ) ;
            psi(3) = fphi( nm, (Dij(i,j)%node(3,1)-x(i) )*overdx,(Dij(i,j)%node(3,2)-y(j))*overdy ) ;
            psi(4) = fphi( nm, (Dij(i,j)%node(4,1)-x(i) )*overdx,(Dij(i,j)%node(4,2)-y(j))*overdy ) ;
            psi(5) = fphi( nm, (Dij(i,j)%node(5,1)-x(i) )*overdx,(Dij(i,j)%node(5,2)-y(j))*overdy ) ;
            psi(6) = fphi( nm, (Dij(i,j)%node(6,1)-x(i) )*overdx,(Dij(i,j)%node(6,2)-y(j))*overdy ) ;
            psi(7) = fphi( nm, (Dij(i,j)%node(7,1)-x(i) )*overdx,(Dij(i,j)%node(7,2)-y(j))*overdy ) ;
            psi(8) = fphi( nm, (Dij(i,j)%node(8,1)-x(i) )*overdx,(Dij(i,j)%node(8,2)-y(j))*overdy ) ;
            psi(9) = fphi( nm, (Dij(i,j)%node(9,1)-x(i) )*overdx,(Dij(i,j)%node(9,2)-y(j))*overdy ) ;
            temp_L(1:6,1:6) = storage_L(i,j,1:6,1:6)
            temp_U(1:6,1:6) = storage_U(i,j,1:6,1:6)
            vert_temp(1:9,1:5) =vert(i,j,1:9,1:5)
            call get_vector2_b( vert_temp(1:9,1:5),psi(1:9),b(1:6) )
            call solve17(temp_L,temp_U,b,aa,6)

        endif

        !aaa(i,j,:) = aa(:)
        do kk = 1,Dij_star(i,j)%nsub_outer
            if(nm == 1)then
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
                xshift = ( x(idx) - x_base ) *overdx
                yshift = ( y(idy) - y_base ) *overdy

                c1_20 = Dij(idx,idy)%umodal(4)
                c1_11 = Dij(idx,idy)%umodal(5)
                c1_02 = Dij(idx,idy)%umodal(6)
                c1_00 = Dij(idx,idy)%umodal(1) -over12* ( c1_20  +c1_02 )
                c1_10 =  Dij(idx,idy)%umodal(2)
                c1_01 =  Dij(idx,idy)%umodal(3)
                xcc = ( x1+x2 )*0.5
                ycc = ( y1+y2 )*0.5


                x2m1 = x2-x1
                y2m1 = y2-y1

                !Dij_star(i,j)%segment_outer(kk)%icurve = 0
                if( abs(x2m1)< 1d-12 .and. abs(y2m1)<1d-12 )then
                    Dij_star(i,j)%segment_outer(kk)%c_ab(1:6) = 0.
                else
                    if( Dij_star(i,j)%segment_outer(kk)%icurve == 0 )then
                        !xgs(1) = xcc + x2m1 * gauss2(1,1)
                        !xgs(2) = xcc + x2m1 * gauss2(2,1)
                        !ygs(1) = ycc + y2m1 * gauss2(1,1)
                        !ygs(2) = ycc + y2m1 * gauss2(2,1)

                        xgs(1) = xcc + x2m1 * gauss3(1,1)
                        xgs(2) = xcc + x2m1 * gauss3(2,1)
                        xgs(3) = xcc + x2m1 * gauss3(3,1)
                        ygs(1) = ycc + y2m1 * gauss3(1,1)
                        ygs(2) = ycc + y2m1 * gauss3(2,1)
                        ygs(3) = ycc + y2m1 * gauss3(3,1)

                        !******************************************
                        !xgs(1) = xcc + x2m1 * gauss2(1,1)
                        !xgs(2) = xcc + x2m1 * gauss2(2,1)
                        !ygs(1) = ycc + y2m1 * gauss2(1,1)
                        !ygs(2) = ycc + y2m1 * gauss2(2,1)
                        !xgs(3) = 0.
                        !ygs(3)= 0.
                        !gauss3(1,2) = 0.5
                        !gauss3(2,2) =  0.5
                        !gauss3(3,2) = 0.
                        !**************************
                        xg1_2 = xgs(1)*xgs(1)
                        xg2_2 = xgs(2)*xgs(2)
                        xg3_2 = xgs(3)*xgs(3)
                        yg1_2 = ygs(1)*ygs(1)
                        yg2_2 = ygs(2)*ygs(2)
                        yg3_2 = ygs(3)*ygs(3)

                        xy1 = xgs(1)*ygs(1)
                        xy2 = xgs(2)*ygs(2)
                        xy3 = xgs(3)*ygs(3)
                        xg1_4 = xg1_2*xg1_2
                        xg2_4 = xg2_2 *xg2_2
                        xg3_4 = xg3_2 *xg3_2

                        yg1_4 = yg1_2 * yg1_2
                        yg2_4 = yg2_2 * yg2_2
                        yg3_4 = yg3_2 * yg3_2
                        if( abs(x2m1) > abs(y2m1) )then
                            slope = (y2m1)/(x2m1)
                            slope_x21 = slope*x2m1
                            green_00 = ( xgs(1) *gauss3(1,2)+xgs(2) *gauss3(2,2) + xgs(3)*gauss3(3,2) ) * slope_x21
                            green_10 = ( xg1_2 *gauss3(1,2) + xg2_2 *gauss3(2,2) + xg3_2*gauss3(3,2) )*0.5 *slope_x21
                            green_01 = -( yg1_2*gauss3(1,2) + yg2_2*gauss3(2,2) + yg3_2*gauss3(3,2) )*0.5*x2m1
                            green_20 = ( xg1_2*xgs(1)*gauss3(1,2)  + xgs(2)*xg2_2*gauss3(2,2)  +xgs(3)*xg3_2*gauss3(3,2) )*one_third *slope_x21
                            green_11 = ( xg1_2*ygs(1)*gauss3(1,2)  + xg2_2*ygs(2)*gauss3(2,2) +xg3_2*ygs(3)*gauss3(3,2)  )*0.5 *slope_x21
                            green_02 = -( ygs(1)*yg1_2*gauss3(1,2) + ygs(2)*yg2_2*gauss3(2,2)  + ygs(3)*yg3_2*gauss3(3,2) )*one_third*x2m1
                            !
                            green_30 = ( xg1_4*gauss3(1,2) + xg2_4*gauss3(2,2)   + xg3_4*gauss3(3,2) )*0.25* slope_x21
                            green_21 = ( xg1_2*xy1*gauss3(1,2) + xg2_2*xy2*gauss3(2,2)  + xg3_2*xy3*gauss3(3,2) )* one_third*slope_x21
                            green_12 = -( xy1*yg1_2*gauss3(1,2) + xy2*yg2_2*gauss3(2,2)  + xy3*yg3_2*gauss3(3,2) )*one_third*x2m1
                            !
                            green_03 =  -( yg1_4*gauss3(1,2) + yg2_4*gauss3(2,2)  +yg3_4*gauss3(3,2) )*0.25*x2m1
                            !
                            green_40 = ( xg1_4*xgs(1)*gauss3(1,2) + xg2_4*xgs(2)*gauss3(2,2) +xg3_4*xgs(3)*gauss3(3,2) )*0.2* slope_x21
                            green_31 = ( xg1_4*ygs(1)*gauss3(1,2) + xg2_4*ygs(2)*gauss3(2,2) + xg3_4*ygs(3)*gauss3(3,2) )*0.25* slope_x21
                            green_22 = ( xgs(1)*xg1_2*yg1_2*gauss3(1,2) + xgs(2)*xg2_2*yg2_2*gauss3(2,2)  + xgs(3)*xg3_2*yg3_2*gauss3(3,2) )*one_third*slope_x21
                            !
                            green_13 = -( yg1_4*xgs(1)*gauss3(1,2) + yg2_4*xgs(2)*gauss3(2,2)  + yg3_4*xgs(3)*gauss3(3,2) )*0.25*x2m1
                            !
                            green_04 = -( ygs(1)*yg1_4*gauss3(1,2) + ygs(2)*yg2_4*gauss3(2,2)  + ygs(3)*yg3_4*gauss3(3,2) )*0.2*x2m1
                            !*****************
                            !green_10 = ( xgs(1)**2 + xgs(2)**2 )*0.25 *slope_x21
                            !green_01 = -( ygs(1)**2 + ygs(2)**2 )*0.25*x2m1
                            !green_20 = ( xgs(1)**3 + xgs(2)**3 )*0.5*one_third *slope_x21
                            !green_11 = ( xgs(1)**2*ygs(1) + xgs(2)**2*ygs(2) )*0.25 *slope_x21
                            !green_02 = -( ygs(1)**3 + ygs(2)**3 )*0.5*one_third*x2m1
                            !
                            !green_30 = ( xgs(1)**4 + xgs(2)**4 )*0.125* slope_x21
                            !green_21 = ( xgs(1)**3*ygs(1) + xgs(2)**3*ygs(2) )*0.5*one_third *slope_x21
                            !green_12 = -( xgs(1)*ygs(1)**3 + xgs(2)*ygs(2)**3 )*0.5*one_third*x2m1
                            !
                            !green_03 =  -( ygs(1)**4 + ygs(2)**4 )*0.125*x2m1
                            !
                            !green_40 = ( xgs(1)**5 + xgs(2)**5 )*0.1* slope_x21
                            !green_31 = ( xgs(1)**4*ygs(1) + xgs(2)**4*ygs(2) )*0.125* slope_x21
                            !green_22 = ( xgs(1)**3*ygs(1)**2 + xgs(2)**3*ygs(2)**2 )*0.5*one_third *slope_x21
                            !
                            !green_13 = -( ygs(1)**4*xgs(1) + ygs(2)**4*xgs(2) )*0.125*x2m1
                            !
                            !green_04 = -( ygs(1)**5 + ygs(2)**5 )*0.1*x2m1
                            Dij_star(i,j)%segment_outer(kk)%c_ab(1) =  c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
                                + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02

                            temp2 = c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
                                + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12
                            Dij_star(i,j)%segment_outer(kk)%c_ab(2) =  xshift* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + temp2

                            temp3 = c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
                                + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03
                            Dij_star(i,j)%segment_outer(kk)%c_ab(3) = yshift* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                +temp3

                            Dij_star(i,j)%segment_outer(kk)%c_ab(4) =  (xshift*xshift-over12)* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + temp2*2.*xshift &
                                +  c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
                                + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22

                            Dij_star(i,j)%segment_outer(kk)%c_ab(5) = xshift*yshift*Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + yshift * temp2 + xshift * temp3 &
                                + c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
                                + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13

                            Dij_star(i,j)%segment_outer(kk)%c_ab(6) = (yshift*yshift-over12)* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + 2. * yshift * temp3 &
                                +  c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
                                + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04
                        else
                            slope = (x2m1)/(y2m1)
                            slope_y21 = slope*y2m1
                            green_00 = ( xgs(1) *gauss3(1,2)+xgs(2) *gauss3(2,2) + xgs(3)*gauss3(3,2) ) * y2m1
                            green_10 = ( xg1_2 *gauss3(1,2) + xg2_2 *gauss3(2,2) + xg3_2*gauss3(3,2) )*0.5 *y2m1
                            green_01 = -( yg1_2*gauss3(1,2) + yg2_2*gauss3(2,2) + yg3_2*gauss3(3,2) )*0.5*slope_y21
                            green_20 = ( xg1_2*xgs(1)*gauss3(1,2)  + xgs(2)*xg2_2*gauss3(2,2)  +xgs(3)*xg3_2*gauss3(3,2) )*one_third *y2m1
                            green_11 = ( xg1_2*ygs(1)*gauss3(1,2)  + xg2_2*ygs(2)*gauss3(2,2) +xg3_2*ygs(3)*gauss3(3,2)  )*0.5 *y2m1
                            green_02 = -( ygs(1)*yg1_2*gauss3(1,2) + ygs(2)*yg2_2*gauss3(2,2)  + ygs(3)*yg3_2*gauss3(3,2) )*one_third*slope_y21
                            !
                            green_30 = ( xg1_4*gauss3(1,2) + xg2_4*gauss3(2,2)   + xg3_4*gauss3(3,2) )*0.25* y2m1
                            green_21 = ( xg1_2*xy1*gauss3(1,2) + xg2_2*xy2*gauss3(2,2)  + xg3_2*xy3*gauss3(3,2) )* one_third*y2m1
                            green_12 = -( xy1*yg1_2*gauss3(1,2) + xy2*yg2_2*gauss3(2,2)  + xy3*yg3_2*gauss3(3,2) )*one_third*slope_y21
                            !
                            green_03 =  -( yg1_4*gauss3(1,2) + yg2_4*gauss3(2,2)  +yg3_4*gauss3(3,2) )*0.25*slope_y21
                            !
                            green_40 = ( xg1_4*xgs(1)*gauss3(1,2) + xg2_4*xgs(2)*gauss3(2,2) +xg3_4*xgs(3)*gauss3(3,2) )*0.2* y2m1
                            green_31 = ( xg1_4*ygs(1)*gauss3(1,2) + xg2_4*ygs(2)*gauss3(2,2) + xg3_4*ygs(3)*gauss3(3,2) )*0.25* y2m1
                            green_22 = ( xgs(1)*xg1_2*yg1_2*gauss3(1,2) + xgs(2)*xg2_2*yg2_2*gauss3(2,2)  + xgs(3)*xg3_2*yg3_2*gauss3(3,2) )*one_third*y2m1
                            !
                            green_13 = -( yg1_4*xgs(1)*gauss3(1,2) + yg2_4*xgs(2)*gauss3(2,2)  + yg3_4*xgs(3)*gauss3(3,2) )*0.25*slope_y21
                            !
                            green_04 = -( ygs(1)*yg1_4*gauss3(1,2) + ygs(2)*yg2_4*gauss3(2,2)  + ygs(3)*yg3_4*gauss3(3,2) )*0.2*slope_y21
                            Dij_star(i,j)%segment_outer(kk)%c_ab(1) =  c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
                                + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02

                            temp2 = c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
                                + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12
                            Dij_star(i,j)%segment_outer(kk)%c_ab(2) =  xshift* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + temp2

                            temp3 = c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
                                + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03
                            Dij_star(i,j)%segment_outer(kk)%c_ab(3) = yshift* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                +temp3

                            Dij_star(i,j)%segment_outer(kk)%c_ab(4) =  (xshift*xshift-over12)* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + temp2*2.*xshift &
                                +  c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
                                + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22

                            Dij_star(i,j)%segment_outer(kk)%c_ab(5) = xshift*yshift*Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + yshift * temp2 + xshift * temp3 &
                                + c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
                                + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13

                            Dij_star(i,j)%segment_outer(kk)%c_ab(6) = (yshift*yshift-over12)* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                                + 2. * yshift * temp3 &
                                +  c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
                                + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04
                        endif
                    else
                        xxi1 = ( Dij_star( i,j )%segment_outer(kk) %v11(1) -x(idx) )/dx ;
                        eta1 = ( Dij_star( i,j )%segment_outer(kk) %v11(2) -y(idy) )/dy ;
                        xxi2 = ( Dij_star( i,j )%segment_outer(kk) %v33(1) -x(idx) )/dx ;
                        eta2 = ( Dij_star( i,j )%segment_outer(kk) %v33(2) -y(idy) )/dy ;
                        xxi3 = ( Dij_star( i,j )%segment_outer(kk) %v22(1) -x(idx) )/dx ;
                        eta3 = ( Dij_star( i,j )%segment_outer(kk) %v22(2) -y(idy) )/dy ;
                        call trans( xxi1,eta1,xxi3,eta3,a_trans,b_trans,c_trans,d_trans )
                        e_x2 = trans_to_ex( a_trans ,b_trans , c_trans, xxi2,eta2  )
                        e_y2 = trans_to_ey( a_trans ,b_trans , d_trans, xxi2,eta2  )


                        ac1 = (eta3-eta1)/2.*e_y2/( e_x2**2 -1. )
                        bc1 = (xxi3 - xxi1)/2.
                        cc1 = (xxi1+xxi3)/2. - (eta3-eta1)/2.*e_y2/( e_x2**2 -1. )

                        ac2 = -(xxi3-xxi1)/2.*e_y2/( e_x2**2 -1. )
                        bc2 = ( eta3 - eta1 )/2.
                        cc2 = (eta1+eta3)/2. + (xxi3-xxi1)/2.*e_y2/( e_x2**2 -1. )

                        exx1 = trans_to_ex( a_trans ,b_trans , c_trans, x1,y1 )
                        exx2 = trans_to_ex( a_trans ,b_trans , c_trans, x2,y2 )
                        dx_xi = exx2 -exx1

                        !xig(1) = (exx1+exx2)*0.5 + dx_xi * gauss2(1,1)
                        !xig(2) = (exx1+exx2)*0.5 + dx_xi * gauss2(2,1)
                        xig(1) = (exx1+exx2)*0.5 + dx_xi * gauss3(1,1)
                        xig(2) = (exx1+exx2)*0.5 + dx_xi * gauss3(2,1)
                        xig(3) = (exx1+exx2)*0.5 + dx_xi * gauss3(3,1)
                    !******************************************
                    !    xig(1) = (exx1+exx2)*0.5 + dx_xi * gauss2(1,1)
                    !    xig(2) = (exx1+exx2)*0.5 + dx_xi * gauss2(2,1)
                    !   xig(3) = 0.
                    !gauss3(1,2) = 0.5
                    !gauss3(2,2) =  0.5
                    !gauss3(3,2) = 0.
                    !**************************
                        !****************************************************************************************************
                        xi_gauss1 = ac1*xig(1)*xig(1) + bc1*xig(1) + cc1
                        dxi_gauss1 = 2.*ac1*xig(1) + bc1
                        xi_gauss2 = ac1*xig(2)*xig(2) + bc1*xig(2) + cc1
                        dxi_gauss2 = 2.*ac1*xig(2) + bc1
                        xi_gauss3 = ac1*xig(3)*xig(3) + bc1*xig(3) + cc1
                        dxi_gauss3 = 2.*ac1*xig(3) + bc1

                        yi_gauss1 = ac2*xig(1)*xig(1) + bc2*xig(1) + cc2
                        dyi_gauss1 = 2.*ac2*xig(1) + bc2
                        yi_gauss2 = ac2*xig(2)*xig(2) + bc2*xig(2) + cc2
                        dyi_gauss2 = 2.*ac2*xig(2) + bc2
                        yi_gauss3 = ac2*xig(3)*xig(3) + bc2*xig(3) + cc2
                        dyi_gauss3 = 2.*ac2*xig(3) + bc2
                        
                        green_00 =  xi_gauss1 * dyi_gauss1*gauss3(1,2) + xi_gauss2 * dyi_gauss2*gauss3(2,2)+ xi_gauss3* dyi_gauss3*gauss3(3,2) 
                        green_10 = ( xi_gauss1 * xi_gauss1 * dyi_gauss1*gauss3(1,2) + xi_gauss2 *xi_gauss2 * dyi_gauss2*gauss3(2,2)  &
                        + xi_gauss3* xi_gauss3 * dyi_gauss3*gauss3(3,2)    )*0.5
                        green_01 = -( yi_gauss1*yi_gauss1*dxi_gauss1*gauss3(1,2) + yi_gauss2*yi_gauss2*dxi_gauss2*gauss3(2,2)        &
                        +  yi_gauss3*yi_gauss3*dxi_gauss3*gauss3(3,2)  )*0.5
                        green_20 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * dyi_gauss1*gauss3(1,2)  &
                            + xi_gauss2*xi_gauss2 *xi_gauss2 * dyi_gauss2*gauss3(2,2) &
                            + xi_gauss3*xi_gauss3 *xi_gauss3 * dyi_gauss3*gauss3(3,2)   )*one_third
                        green_11 = (  xi_gauss1 * xi_gauss1*yi_gauss1 * dyi_gauss1 *gauss3(1,2) &
                            + xi_gauss2 *xi_gauss2*yi_gauss2 * dyi_gauss2*gauss3(2,2)   &
                            + xi_gauss3 *xi_gauss3*yi_gauss3 * dyi_gauss3*gauss3(3,2)   )*0.5
                        green_02 = -(  yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gauss3(1,2) &
                            + yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gauss3(2,2)  &
                            + yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gauss3(3,2)    )* one_third
                        !
                        green_30 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1* dyi_gauss1 *gauss3(1,2)  &
                            + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2 * dyi_gauss2*gauss3(2,2) &
                            + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3 * dyi_gauss3*gauss3(3,2)   )*0.25
                        green_21 = ( xi_gauss1*xi_gauss1 * xi_gauss1*yi_gauss1 * dyi_gauss1*gauss3(1,2)  &
                            + xi_gauss2*xi_gauss2 *xi_gauss2*yi_gauss2 * dyi_gauss2*gauss3(2,2)   &
                            + xi_gauss3*xi_gauss3 *xi_gauss3*yi_gauss3* dyi_gauss3*gauss3(3,2)   )*one_third
                        green_12 = -( xi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gauss3(1,2) &
                            + xi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gauss3(2,2)   &
                            + xi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gauss3(3,2)    )*one_third
                        !
                        green_03 =  -(  yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gauss3(1,2) &
                            + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gauss3(2,2)     &
                            + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gauss3(3,2)  )*0.25
                        !
                        green_40 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1* xi_gauss1* dyi_gauss1*gauss3(1,2)  &
                            + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2*xi_gauss2 * dyi_gauss2*gauss3(2,2)  &
                            + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3*xi_gauss3 * dyi_gauss3*gauss3(3,2)   )*0.2
                        green_31 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1*yi_gauss1* dyi_gauss1*gauss3(1,2)  &
                            + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2*yi_gauss2* dyi_gauss2*gauss3(2,2) &
                            + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3*yi_gauss3* dyi_gauss3*gauss3(3,2)   )*0.25
                        green_22 = ( xi_gauss1*xi_gauss1 * xi_gauss1*yi_gauss1*yi_gauss1* dyi_gauss1*gauss3(1,2)  &
                            + xi_gauss2*xi_gauss2 *xi_gauss2*yi_gauss2*yi_gauss2* dyi_gauss2*gauss3(2,2)    &
                            + xi_gauss3*xi_gauss3 *xi_gauss3*yi_gauss3*yi_gauss3* dyi_gauss3*gauss3(3,2)  )*one_third
                        !
                        green_13 = -( yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*xi_gauss1*dxi_gauss1*gauss3(1,2) &
                            + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*xi_gauss2*dxi_gauss2*gauss3(2,2)  &
                            + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*xi_gauss3*dxi_gauss3*gauss3(3,2) )*0.25
                        !
                        green_04 = -( yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gauss3(1,2) &
                            + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gauss3(2,2)        &
                            + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gauss3(3,2)    )*0.2

                        Dij_star(i,j)%segment_outer(kk)%c_ab(1) =  (c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
                            + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02)*dx_xi

                        temp2 = (c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
                            + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12)*dx_xi
                        Dij_star(i,j)%segment_outer(kk)%c_ab(2) =  xshift* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + temp2

                        temp3 = (c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
                            + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03 )*dx_xi
                        Dij_star(i,j)%segment_outer(kk)%c_ab(3) = yshift* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            +temp3

                        Dij_star(i,j)%segment_outer(kk)%c_ab(4) =  (xshift*xshift-over12)* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + temp2*2.*xshift &
                            +  (c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
                            + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22)*dx_xi

                        Dij_star(i,j)%segment_outer(kk)%c_ab(5) = xshift*yshift*Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + yshift * temp2 + xshift * temp3 &
                            + (c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
                            + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13)*dx_xi

                        Dij_star(i,j)%segment_outer(kk)%c_ab(6) = (yshift*yshift-over12)* Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + 2. * yshift * temp3 &
                            +  (c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
                            + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04)*dx_xi

                        !print *,Dij_star(i,j)%segment_outer(kk)%c_ab(1:6)

                    endif

                endif
            endif

            sum = sum + aa(1)*Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                + aa(2)*Dij_star(i,j)%segment_outer(kk)%c_ab(2) + aa(3)*Dij_star(i,j)%segment_outer(kk)%c_ab(3) &
                + aa(4)*Dij_star(i,j)%segment_outer(kk)%c_ab(4) &
                + aa(5)*Dij_star(i,j)%segment_outer(kk)%c_ab(5) + aa(6)*Dij_star(i,j)%segment_outer(kk)%c_ab(6)
        enddo

        do kk = 1,Dij_star(i,j)%nsub_inner
            if(nm == 1)then
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
                xshift = ( x(idx) - x_base ) *overdx
                yshift = ( y(idy) - y_base ) *overdy

                c1_20 = Dij(idx,idy)%umodal(4)
                c1_11 = Dij(idx,idy)%umodal(5)
                c1_02 = Dij(idx,idy)%umodal(6)
                c1_00 = Dij(idx,idy)%umodal(1) -over12* ( c1_20  +c1_02 )
                c1_10 =  Dij(idx,idy)%umodal(2)
                c1_01 =  Dij(idx,idy)%umodal(3)
                xcc = ( x1+x2 )*0.5
                ycc = ( y1+y2 )*0.5


                x2m1 = x2-x1
                y2m1 = y2-y1
                !xgs(1) = xcc + x2m1 * gauss2(1,1)
                !xgs(2) = xcc + x2m1 * gauss2(2,1)
                !ygs(1) = ycc + y2m1 * gauss2(1,1)
                !ygs(2) = ycc + y2m1 * gauss2(2,1)
                !xg1_2 = xgs(1)*xgs(1)
                !xg2_2 = xgs(2)*xgs(2)
                !yg1_2 = ygs(1)*ygs(1)
                !yg2_2 = ygs(2)*ygs(2)
                !xy1 = xgs(1)*ygs(1)
                !xy2 = xgs(2)*ygs(2)
                !xg1_4 = xg1_2*xg1_2
                !xg2_4 = xg2_2 *xg2_2
                !yg1_4 = yg1_2 * yg1_2
                !yg2_4 = yg2_2 * yg2_2
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
                    Dij_star(i,j)%segment_inner(kk)%c_ab(1:6) = 0.
                else
                    !xgs(1) = xcc + x2m1 * gauss2(1,1)
                    !xgs(2) = xcc + x2m1 * gauss2(2,1)
                    !ygs(1) = ycc + y2m1 * gauss2(1,1)
                    !ygs(2) = ycc + y2m1 * gauss2(2,1)

                    xgs(1) = xcc + x2m1 * gauss3(1,1)
                    xgs(2) = xcc + x2m1 * gauss3(2,1)
                    xgs(3) = xcc + x2m1 * gauss3(3,1)
                    ygs(1) = ycc + y2m1 * gauss3(1,1)
                    ygs(2) = ycc + y2m1 * gauss3(2,1)
                    ygs(3) = ycc + y2m1 * gauss3(3,1)

                    !******************************************
                    !xgs(1) = xcc + x2m1 * gauss2(1,1)
                    !xgs(2) = xcc + x2m1 * gauss2(2,1)
                    !ygs(1) = ycc + y2m1 * gauss2(1,1)
                    !ygs(2) = ycc + y2m1 * gauss2(2,1)
                    !xgs(3) = 0.
                    !ygs(3)= 0.
                    !gauss3(1,2) = 0.5
                    !gauss3(2,2) =  0.5
                    !gauss3(3,2) = 0.
                    !**************************
                    xg1_2 = xgs(1)*xgs(1)
                    xg2_2 = xgs(2)*xgs(2)
                    xg3_2 = xgs(3)*xgs(3)
                    yg1_2 = ygs(1)*ygs(1)
                    yg2_2 = ygs(2)*ygs(2)
                    yg3_2 = ygs(3)*ygs(3)

                    xy1 = xgs(1)*ygs(1)
                    xy2 = xgs(2)*ygs(2)
                    xy3 = xgs(3)*ygs(3)
                    xg1_4 = xg1_2*xg1_2
                    xg2_4 = xg2_2 *xg2_2
                    xg3_4 = xg3_2 *xg3_2

                    yg1_4 = yg1_2 * yg1_2
                    yg2_4 = yg2_2 * yg2_2
                    yg3_4 = yg3_2 * yg3_2
                    if( abs(x2m1) > abs(y2m1) )then
                        slope = (y2m1)/(x2m1)
                        slope_x21 = slope*x2m1
                        green_00 = ( xgs(1) *gauss3(1,2)+xgs(2) *gauss3(2,2) + xgs(3)*gauss3(3,2) ) * slope_x21
                        green_10 = ( xg1_2 *gauss3(1,2) + xg2_2 *gauss3(2,2) + xg3_2*gauss3(3,2) )*0.5 *slope_x21
                        green_01 = -( yg1_2*gauss3(1,2) + yg2_2*gauss3(2,2) + yg3_2*gauss3(3,2) )*0.5*x2m1
                        green_20 = ( xg1_2*xgs(1)*gauss3(1,2)  + xgs(2)*xg2_2*gauss3(2,2)  +xgs(3)*xg3_2*gauss3(3,2) )*one_third *slope_x21
                        green_11 = ( xg1_2*ygs(1)*gauss3(1,2)  + xg2_2*ygs(2)*gauss3(2,2) +xg3_2*ygs(3)*gauss3(3,2)  )*0.5 *slope_x21
                        green_02 = -( ygs(1)*yg1_2*gauss3(1,2) + ygs(2)*yg2_2*gauss3(2,2)  + ygs(3)*yg3_2*gauss3(3,2) )*one_third*x2m1
                        !
                        green_30 = ( xg1_4*gauss3(1,2) + xg2_4*gauss3(2,2)   + xg3_4*gauss3(3,2) )*0.25* slope_x21
                        green_21 = ( xg1_2*xy1*gauss3(1,2) + xg2_2*xy2*gauss3(2,2)  + xg3_2*xy3*gauss3(3,2) )* one_third*slope_x21
                        green_12 = -( xy1*yg1_2*gauss3(1,2) + xy2*yg2_2*gauss3(2,2)  + xy3*yg3_2*gauss3(3,2) )*one_third*x2m1
                        !
                        green_03 =  -( yg1_4*gauss3(1,2) + yg2_4*gauss3(2,2)  +yg3_4*gauss3(3,2) )*0.25*x2m1
                        !
                        green_40 = ( xg1_4*xgs(1)*gauss3(1,2) + xg2_4*xgs(2)*gauss3(2,2) +xg3_4*xgs(3)*gauss3(3,2) )*0.2* slope_x21
                        green_31 = ( xg1_4*ygs(1)*gauss3(1,2) + xg2_4*ygs(2)*gauss3(2,2) + xg3_4*ygs(3)*gauss3(3,2) )*0.25* slope_x21
                        green_22 = ( xgs(1)*xg1_2*yg1_2*gauss3(1,2) + xgs(2)*xg2_2*yg2_2*gauss3(2,2)  + xgs(3)*xg3_2*yg3_2*gauss3(3,2) )*one_third*slope_x21
                        !
                        green_13 = -( yg1_4*xgs(1)*gauss3(1,2) + yg2_4*xgs(2)*gauss3(2,2)  + yg3_4*xgs(3)*gauss3(3,2) )*0.25*x2m1
                        !
                        green_04 = -( ygs(1)*yg1_4*gauss3(1,2) + ygs(2)*yg2_4*gauss3(2,2)  + ygs(3)*yg3_4*gauss3(3,2) )*0.2*x2m1
                        !*****************
                        !green_10 = ( xgs(1)**2 + xgs(2)**2 )*0.25 *slope_x21
                        !green_01 = -( ygs(1)**2 + ygs(2)**2 )*0.25*x2m1
                        !green_20 = ( xgs(1)**3 + xgs(2)**3 )*0.5*one_third *slope_x21
                        !green_11 = ( xgs(1)**2*ygs(1) + xgs(2)**2*ygs(2) )*0.25 *slope_x21
                        !green_02 = -( ygs(1)**3 + ygs(2)**3 )*0.5*one_third*x2m1
                        !
                        !green_30 = ( xgs(1)**4 + xgs(2)**4 )*0.125* slope_x21
                        !green_21 = ( xgs(1)**3*ygs(1) + xgs(2)**3*ygs(2) )*0.5*one_third *slope_x21
                        !green_12 = -( xgs(1)*ygs(1)**3 + xgs(2)*ygs(2)**3 )*0.5*one_third*x2m1
                        !
                        !green_03 =  -( ygs(1)**4 + ygs(2)**4 )*0.125*x2m1
                        !
                        !green_40 = ( xgs(1)**5 + xgs(2)**5 )*0.1* slope_x21
                        !green_31 = ( xgs(1)**4*ygs(1) + xgs(2)**4*ygs(2) )*0.125* slope_x21
                        !green_22 = ( xgs(1)**3*ygs(1)**2 + xgs(2)**3*ygs(2)**2 )*0.5*one_third *slope_x21
                        !
                        !green_13 = -( ygs(1)**4*xgs(1) + ygs(2)**4*xgs(2) )*0.125*x2m1
                        !
                        !green_04 = -( ygs(1)**5 + ygs(2)**5 )*0.1*x2m1
                        Dij_star(i,j)%segment_inner(kk)%c_ab(1) =  c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
                            + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02

                        temp2 = c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
                            + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12
                        Dij_star(i,j)%segment_inner(kk)%c_ab(2) =  xshift* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + temp2

                        temp3 = c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
                            + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03
                        Dij_star(i,j)%segment_inner(kk)%c_ab(3) = yshift* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            +temp3

                        Dij_star(i,j)%segment_inner(kk)%c_ab(4) =  (xshift*xshift-over12)* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + temp2*2.*xshift &
                            +  c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
                            + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22

                        Dij_star(i,j)%segment_inner(kk)%c_ab(5) = xshift*yshift*Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + yshift * temp2 + xshift * temp3 &
                            + c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
                            + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13

                        Dij_star(i,j)%segment_inner(kk)%c_ab(6) = (yshift*yshift-over12)* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + 2. * yshift * temp3 &
                            +  c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
                            + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04
                    else
                        slope = (x2m1)/(y2m1)
                        slope_y21 = slope*y2m1
                        green_00 = ( xgs(1) *gauss3(1,2)+xgs(2) *gauss3(2,2) + xgs(3)*gauss3(3,2) ) * y2m1
                        green_10 = ( xg1_2 *gauss3(1,2) + xg2_2 *gauss3(2,2) + xg3_2*gauss3(3,2) )*0.5 *y2m1
                        green_01 = -( yg1_2*gauss3(1,2) + yg2_2*gauss3(2,2) + yg3_2*gauss3(3,2) )*0.5*slope_y21
                        green_20 = ( xg1_2*xgs(1)*gauss3(1,2)  + xgs(2)*xg2_2*gauss3(2,2)  +xgs(3)*xg3_2*gauss3(3,2) )*one_third *y2m1
                        green_11 = ( xg1_2*ygs(1)*gauss3(1,2)  + xg2_2*ygs(2)*gauss3(2,2) +xg3_2*ygs(3)*gauss3(3,2)  )*0.5 *y2m1
                        green_02 = -( ygs(1)*yg1_2*gauss3(1,2) + ygs(2)*yg2_2*gauss3(2,2)  + ygs(3)*yg3_2*gauss3(3,2) )*one_third*slope_y21
                        !
                        green_30 = ( xg1_4*gauss3(1,2) + xg2_4*gauss3(2,2)   + xg3_4*gauss3(3,2) )*0.25* y2m1
                        green_21 = ( xg1_2*xy1*gauss3(1,2) + xg2_2*xy2*gauss3(2,2)  + xg3_2*xy3*gauss3(3,2) )* one_third*y2m1
                        green_12 = -( xy1*yg1_2*gauss3(1,2) + xy2*yg2_2*gauss3(2,2)  + xy3*yg3_2*gauss3(3,2) )*one_third*slope_y21
                        !
                        green_03 =  -( yg1_4*gauss3(1,2) + yg2_4*gauss3(2,2)  +yg3_4*gauss3(3,2) )*0.25*slope_y21
                        !
                        green_40 = ( xg1_4*xgs(1)*gauss3(1,2) + xg2_4*xgs(2)*gauss3(2,2) +xg3_4*xgs(3)*gauss3(3,2) )*0.2* y2m1
                        green_31 = ( xg1_4*ygs(1)*gauss3(1,2) + xg2_4*ygs(2)*gauss3(2,2) + xg3_4*ygs(3)*gauss3(3,2) )*0.25* y2m1
                        green_22 = ( xgs(1)*xg1_2*yg1_2*gauss3(1,2) + xgs(2)*xg2_2*yg2_2*gauss3(2,2)  + xgs(3)*xg3_2*yg3_2*gauss3(3,2) )*one_third*y2m1
                        !
                        green_13 = -( yg1_4*xgs(1)*gauss3(1,2) + yg2_4*xgs(2)*gauss3(2,2)  + yg3_4*xgs(3)*gauss3(3,2) )*0.25*slope_y21
                        !
                        green_04 = -( ygs(1)*yg1_4*gauss3(1,2) + ygs(2)*yg2_4*gauss3(2,2)  + ygs(3)*yg3_4*gauss3(3,2) )*0.2*slope_y21

                        Dij_star(i,j)%segment_inner(kk)%c_ab(1) =  c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
                            + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02

                        temp2 = c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
                            + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12
                        Dij_star(i,j)%segment_inner(kk)%c_ab(2) =  xshift* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + temp2

                        temp3 = c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
                            + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03
                        Dij_star(i,j)%segment_inner(kk)%c_ab(3) = yshift* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            +temp3

                        Dij_star(i,j)%segment_inner(kk)%c_ab(4) =  (xshift*xshift-over12)* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + temp2*2.*xshift &
                            +  c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
                            + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22

                        Dij_star(i,j)%segment_inner(kk)%c_ab(5) = xshift*yshift*Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + yshift * temp2 + xshift * temp3 &
                            + c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
                            + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13

                        Dij_star(i,j)%segment_inner(kk)%c_ab(6) = (yshift*yshift-over12)* Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + 2. * yshift * temp3 &
                            +  c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
                            + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04
                    endif
                endif
            endif

            sum = sum + aa(1)*Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                + aa(2)*Dij_star(i,j)%segment_inner(kk)%c_ab(2) + aa(3)*Dij_star(i,j)%segment_inner(kk)%c_ab(3) &
                + aa(4)*Dij_star(i,j)%segment_inner(kk)%c_ab(4) &
                + aa(5)*Dij_star(i,j)%segment_inner(kk)%c_ab(5) + aa(6)*Dij_star(i,j)%segment_inner(kk)%c_ab(6)

        enddo
        sum16(nm) = sum * ai(nm)

    enddo !nm

    end subroutine green_p2QC_gauss3_eulerian