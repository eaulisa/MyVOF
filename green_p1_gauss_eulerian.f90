    subroutine green_p1_gauss_eulerian(sum16)
    ! The subroutine is the factorization version of
    ! subroutine green_p1_gauss17
    implicit none
    integer :: kk
    real :: x1,y1,x2,y2
    real :: sum16(3)
    real :: sum

    real :: x11,y11,x22,y22
    integer :: idx,idy
    !*********** least square ***********
    real :: A_temp(3,3)
    real :: temp_L(3,3),temp_U(3,3)
    real :: b(3),aa(3),psi(4)
    real :: vert_temp(1:4,1:2)
    !real,allocatable :: storage_L(:,:,:,:)
    !real,allocatable :: storage_U(:,:,:,:)
    !real,allocatable :: vert(:,:,:,:)
    !real,allocatable :: aaa(:,:,:)
    !integer,allocatable :: iswitch(:,:)
    !************************************
    integer :: nm

    real :: xshift,yshift,c1_10,c1_01,c1_00,xcc,ycc,slope
    real :: x2m1,y2m1,xcc2,ycc2

    real :: slope_x21,slope_y21
    real :: green_00,green_10,green_01,green_20,green_11,green_02
    real :: x_base,y_base

    x_base = 0.25*( Dij_star(i,j)%vertex_star(1,1) + Dij_star(i,j)%vertex_star(2,1)  &
        + Dij_star(i,j)%vertex_star(3,1)  + Dij_star(i,j)%vertex_star(4,1)  )
    y_base = 0.25*( Dij_star(i,j)%vertex_star(1,2) + Dij_star(i,j)%vertex_star(2,2)  &
        + Dij_star(i,j)%vertex_star(3,2)  + Dij_star(i,j)%vertex_star(4,2)  )

    do nm = 1 ,3
        sum = 0.
        !
        if(nm == 1)then
            aa(1) = 1.
            aa(2:3) = 0.
        else
            if( nm ==2 )then
                vert_temp(1,1) = ( Dij_star(i,j)%vertex_star(1,1) - x_base )*overdx
                vert_temp(2,1) = ( Dij_star(i,j)%vertex_star(2,1) - x_base )*overdx
                vert_temp(3,1) = ( Dij_star(i,j)%vertex_star(3,1) - x_base )*overdx
                vert_temp(4,1) = ( Dij_star(i,j)%vertex_star(4,1) - x_base )*overdx

                vert_temp(1,2) = ( Dij_star(i,j)%vertex_star(1,2) - y_base )*overdy
                vert_temp(2,2) = ( Dij_star(i,j)%vertex_star(2,2) - y_base )*overdy
                vert_temp(3,2) = ( Dij_star(i,j)%vertex_star(3,2) - y_base )*overdy
                vert_temp(4,2) = ( Dij_star(i,j)%vertex_star(4,2) - y_base )*overdy
                call get_matrix_a( vert_temp(1:4,1:2),  A_temp(1:3,1:3) )

                call doolittle(A_temp,temp_L,temp_U,3)

                storage_L(i,j,1:3,1:3) = temp_L(1:3,1:3)
                storage_U(i,j,1:3,1:3) = temp_U(1:3,1:3)
                vert(i,j,1:4,1:2) = vert_temp(1:4,1:2)
            endif

            psi(1) = fphi( nm, (Dij(i,j)%vertex(1,1)-x(i) )*overdx,(Dij(i,j)%vertex(1,2)-y(j))*overdy ) ;
            psi(2) = fphi( nm, (Dij(i,j)%vertex(2,1)-x(i) )*overdx,(Dij(i,j)%vertex(2,2)-y(j))*overdy ) ;
            psi(3) = fphi( nm, (Dij(i,j)%vertex(3,1)-x(i) )*overdx,(Dij(i,j)%vertex(3,2)-y(j))*overdy ) ;
            psi(4) = fphi( nm, (Dij(i,j)%vertex(4,1)-x(i) )*overdx,(Dij(i,j)%vertex(4,2)-y(j))*overdy ) ;
            temp_L(1:3,1:3) = storage_L(i,j,1:3,1:3)
            temp_U(1:3,1:3) = storage_U(i,j,1:3,1:3)
            vert_temp(1:4,1:2) =vert( i,j ,1:4,1:2)
            call get_vector_b( vert_temp(1:4,1:2),psi(1:4),b(1:3) )
            call solve17(temp_L,temp_U,b,aa,3)
        endif


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
                c1_00 = Dij(idx,idy)%umodal(1)
                c1_10 = Dij(idx,idy)%umodal(2)
                c1_01 = Dij(idx,idy)%umodal(3)

                xcc = ( x1+x2 )*0.5
                ycc = ( y1+y2 )*0.5
                xcc2 = xcc*xcc
                ycc2 = ycc*ycc
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
                    Dij_star(i,j)%segment_outer(kk)%c_ab(1:3) = 0.
                else
                    if( abs(x2m1) > abs(y2m1) )then
                        slope = (y2m1)/(x2m1)
                        
                        slope_x21 = slope*x2m1
                        !green_00 = xcc* slope_x21
                        !green_10 = xcc2 *0.5*slope_x21
                        !green_01 = -ycc2*0.5*x2m1
                        !green_20 = xcc2*xcc*one_third*slope_x21
                        !green_11 = xcc2*ycc*0.5* slope_x21
                        !green_02 = -ycc2*ycc*one_third*x2m1
                        
                        green_00 = xcc* slope_x21
                        green_10 = xcc2 *0.5*slope_x21
                        green_01 = xcc*ycc*slope_x21
                        green_20 = xcc2*xcc*one_third*slope_x21
                        green_11 = xcc2*ycc*0.5* slope_x21
                        green_02 = -ycc2*ycc*one_third*x2m1
                        Dij_star(i,j)%segment_outer(kk)%c_ab(1) = c1_00 * green_00 + c1_10*green_10 + c1_01*green_01
                        !print *, c1_00 * green_00 + c1_10*green_10 + c1_01*green_01
                        !pause
                        Dij_star(i,j)%segment_outer(kk)%c_ab(2) = xshift * Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + c1_00*green_10 + c1_10*green_20 + c1_01*green_11
                        Dij_star(i,j)%segment_outer(kk)%c_ab(3) = yshift * Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + c1_00*green_01 + c1_10*green_11 + c1_01*green_02
                    else
                        slope = (x2m1)/(y2m1)
                        slope_y21 = slope* y2m1
                        !green_00 = xcc* y2m1
                        !green_10 = xcc2 *0.5* y2m1
                        !green_01 = -ycc2*0.5* slope_y21
                        !green_20 = xcc2*xcc*one_third* y2m1
                        !green_11 = xcc2*ycc*0.5* y2m1
                        !green_02 = -ycc2*ycc*one_third* slope_y21

                        green_00 = xcc* y2m1
                        green_10 = xcc2 *0.5* y2m1
                        green_01 = xcc*ycc* y2m1
                        green_20 = xcc2*xcc*one_third* y2m1
                        green_11 = xcc2*ycc*0.5* y2m1
                        green_02 = -ycc2*ycc*one_third* slope_y21
                        
                        Dij_star(i,j)%segment_outer(kk)%c_ab(1) = c1_00 * green_00 + c1_10*green_10 + c1_01*green_01
                        Dij_star(i,j)%segment_outer(kk)%c_ab(2) = xshift * Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + c1_00*green_10 + c1_10*green_20 + c1_01*green_11
                        Dij_star(i,j)%segment_outer(kk)%c_ab(3) = yshift * Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                            + c1_00*green_01 + c1_10*green_11 + c1_01*green_02
                    endif
                endif
            endif

            sum = sum + aa(1)*Dij_star(i,j)%segment_outer(kk)%c_ab(1) &
                + aa(2)*Dij_star(i,j)%segment_outer(kk)%c_ab(2) + aa(3)*Dij_star(i,j)%segment_outer(kk)%c_ab(3)

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
                yshift = ( y(idy) - y_base  ) *overdy
                c1_00 = Dij(idx,idy)%umodal(1)
                c1_10 = Dij(idx,idy)%umodal(2)
                c1_01 = Dij(idx,idy)%umodal(3)

                xcc = ( x1+x2 )*0.5
                ycc = ( y1+y2 )*0.5
                xcc2 = xcc*xcc
                ycc2 = ycc*ycc
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
                    Dij_star(i,j)%segment_inner(kk)%c_ab(1:3) = 0.
                else
                    if( abs(x2m1) > abs(y2m1) )then
                        slope = (y2m1)/(x2m1)
                        slope_x21 = slope*x2m1
                        !green_00 = xcc* slope_x21
                        !green_10 = xcc2 *0.5*slope_x21
                        !green_01 = -ycc2*0.5*x2m1
                        !green_20 = xcc2*xcc*one_third*slope_x21
                        !green_11 = xcc2*ycc*0.5* slope_x21
                        !green_02 = -ycc2*ycc*one_third*x2m1
                        
                        green_00 = xcc* slope_x21
                        green_10 = xcc2 *0.5*slope_x21
                        green_01 = xcc*ycc*slope_x21
                        green_20 = xcc2*xcc*one_third*slope_x21
                        green_11 = xcc2*ycc*0.5* slope_x21
                        green_02 = -ycc2*ycc*one_third*x2m1
                        Dij_star(i,j)%segment_inner(kk)%c_ab(1) = c1_00 * green_00 + c1_10*green_10 + c1_01*green_01
                        Dij_star(i,j)%segment_inner(kk)%c_ab(2) = xshift * Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + c1_00*green_10 + c1_10*green_20 + c1_01*green_11
                        Dij_star(i,j)%segment_inner(kk)%c_ab(3) = yshift * Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + c1_00*green_01 + c1_10*green_11 + c1_01*green_02
                    else
                        slope = (x2m1)/(y2m1)
                        slope_y21 = slope* y2m1
                        !green_00 = xcc* y2m1
                        !green_10 = xcc2 *0.5* y2m1
                        !green_01 = -ycc2*0.5* slope_y21
                        !green_20 = xcc2*xcc*one_third* y2m1
                        !green_11 = xcc2*ycc*0.5* y2m1
                        !green_02 = -ycc2*ycc*one_third* slope_y21

                        green_00 = xcc* y2m1
                        green_10 = xcc2 *0.5* y2m1
                        green_01 = xcc*ycc* y2m1
                        green_20 = xcc2*xcc*one_third* y2m1
                        green_11 = xcc2*ycc*0.5* y2m1
                        green_02 = -ycc2*ycc*one_third* slope_y21
                        Dij_star(i,j)%segment_inner(kk)%c_ab(1) = c1_00 * green_00 + c1_10*green_10 + c1_01*green_01
                        Dij_star(i,j)%segment_inner(kk)%c_ab(2) = xshift * Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + c1_00*green_10 + c1_10*green_20 + c1_01*green_11
                        Dij_star(i,j)%segment_inner(kk)%c_ab(3) = yshift * Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                            + c1_00*green_01 + c1_10*green_11 + c1_01*green_02
                    endif
                endif
            endif

            sum = sum + aa(1)*Dij_star(i,j)%segment_inner(kk)%c_ab(1) &
                + aa(2)*Dij_star(i,j)%segment_inner(kk)%c_ab(2) + aa(3)*Dij_star(i,j)%segment_inner(kk)%c_ab(3)
            

        enddo
        sum16(nm) = sum * ai(nm)
        
        

    enddo !nm

    end subroutine green_p1_gauss_eulerian
