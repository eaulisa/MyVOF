    subroutine search_QC_formula
    implicit none
    real :: xmax,xmin,ymax,ymin
    integer :: im

    integer :: isx(10),isy(10)
    integer :: iam,ibm,icm,idm
    integer :: kside,kside1
    real :: vleft, vright,vleft_y,vright_y
    real :: vleft1, vright1, vleft_y1,vright_y1
    integer :: ia,ib,mx
    integer :: ic,id,my
    integer :: ia1,ib1,ic1,id1

    real :: zz(8,2)
    real :: point_inner_x( 7,2,2 ),point_inner_y( 7,2,2 )
    real :: point_inner_x_temp( 7,4,2 ), point_inner_y_temp( 7,4,2 )
    integer :: ix,iy
    integer :: ii,jj

    integer :: idx,idy,nby,nbx,ix0,iy0
    real ::  umod_temp(200,200,6)

    !
    integer :: inter
    integer :: inner_convex_x, inner_convex_y
    integer :: inner_concave_x(10),inner_concave_y(10)
    !*****
    real :: x11,y11,x22,y22,x33,y33,x44,y44
    real :: a1,b1,c1,a2,b2,c2

    real :: a_trans,b_trans,c_trans,d_trans
    real :: e_x2,e_y2,e_yc
    real :: at1,bt1,ct1,at2,bt2,ct2
    real :: atemp1,btemp1,atemp2,btemp2
    real :: x00,y00
    integer :: kn
    real :: ac,bc,cc
    real :: e_xi(2),e_ta(2)
    integer :: ie
    real :: point_inner_convex_x(20,2,2),point_inner_convex_y(20,2,2)
    real :: point_inner_concave_x(20,4,2),point_inner_concave_y(20,4,2)
    real :: compare1,compare2
    real :: zztemp(2)
    real :: e_temp1,e_temp2,e_tempxc,e_tempyc
    real :: xid,yid
    real :: temp(6,2),temp_store(2)
    integer :: ix_concave,iy_concave
    real :: point_segment(20,2)

    integer :: i1111,i2222

    integer :: isearch_extend

    real :: xleft_id,xright_id
    real :: yleft_id,yright_id

    real :: eps_2d8

    real :: a16,b16,c16,distance ! switch to straight line
    real :: distance1,distance2

    integer :: ii17,jj17

    integer :: icurve_temp


    integer :: n_root

    isearch_extend =  0

    eps_2d8 = 2d-7

    do i = 1-1 , nx +1
        do j = 1 -1 , ny+1
            x11 = Dij_star(i,j)%vertex_star(1,1)
            y11 = Dij_star(i,j)%vertex_star(1,2)
            x22 = Dij_star(i,j)%vertex_star(2,1)
            y22 = Dij_star(i,j)%vertex_star(2,2)
            x33 = Dij_star(i,j)%vertex_star(3,1)
            y33 = Dij_star(i,j)%vertex_star(3,2)
            x44 = Dij_star(i,j)%vertex_star(4,1)
            y44 = Dij_star(i,j)%vertex_star(4,2)
            call get_line( x11,y11,x33,y33,  a1,b1,c1 )
            call get_line( x22,y22,x44,y44,  a2,b2,c2 )
            center(i,j,1) = (b1*c2-b2*c1)/(a1*b2-a2*b1)
            center(i,j,2) = -(a1*c2-a2*c1)/(a1*b2-a2*b1)

            x_id_up(i,j,1) = min(x11,x22,x33,x44)
            x_id_up(i,j,2) = max(x11,x22,x33,x44)

            y_id_up(i,j,1) = min(y11,y22,y33,y44)
            y_id_up(i,j,2) = max(y11,y22,y33,y44)
        enddo
    enddo


    do i = 1, nx
        do j = 1 , ny
            ! find outer line segments and inner line segments.

            ! preparation for inner segments
            xmax = -1000.
            ymax = -1000.
            do im = 1,4
                xmax = max( xmax, Dij_star(i,j)%vertex_star(im,1)  )
                ymax = max( ymax, Dij_star(i,j)%vertex_star(im,2)  )
            enddo

            xmin = 1000.
            ymin = 1000.
            do im = 1,4
                xmin = min( xmin, Dij_star(i,j)%vertex_star(im,1)  )
                ymin = min( ymin, Dij_star(i,j)%vertex_star(im,2)  )
            enddo


            call id_get17( xmin,dx,i,iam,1 )
            call id_get17( xmax,dx,i,ibm,1 )
            !iam = id_get( (xmin-xleft)/dx  ) ; ibm = id_get( (xmax-xleft)/dx ) ;

            isx(:) = 0

            call id_get17( ymin,dy,j,icm,2 );
            call id_get17( ymax,dy,j,idm,2 );

            isy(:) = 0
            !******************************************
            Dij_star(i,j)%nsub_outer = 0
            Dij_star(i,j)%nsub_inner = 0
            inner_convex_x = 0; inner_convex_y = 0;
            inner_concave_x(:) = 0; inner_concave_y(:) = 0;

            do kside = 1,4
                ! connect sides: side1, 1--2, side2, 2--3, side3, 3--4, side4, 4--1
                kside1 = kside + 1
                if(kside ==4 )then
                    kside1 = 1
                endif

                !
                ! Step 1: find intersections between side_k and x=x_i

                ! project side_k and x=x_i to the coordinate e_x - e_y
                !
                x11 = Dij_star(i,j)%vertex_star(kside,1)
                y11 = Dij_star(i,j)%vertex_star(kside,2)
                x22 = Dij_star(i,j)%vertex_star(kside1,1)
                y22 = Dij_star(i,j)%vertex_star(kside1,2)

                x33 = Dij_star(i,j)%node_star(2*kside,1)
                y33 = Dij_star(i,j)%node_star(2*kside,2)
                call trans( x11,y11,x22,y22,a_trans,b_trans,c_trans,d_trans )
                !a_trans = 2.*( x22-x11 )/( (x11-x22)**2 + (y11-y22)**2 )
                !b_trans = 2.*( y22-y11 )/( (x11-x22)**2 + (y11-y22)**2 )
                !c_trans = (x11**2-x22**2 +y11**2-y22**2 )/( (x11-x22)**2 + (y11-y22)**2 )
                !d_trans = 2.*( x22*y11 - x11*y22 )/( (x11-x22)**2 + (y11-y22)**2 )
                !if( abs(x22-x11)< abs(y22-y11) )then
                !    a_trans = 2.*( (x22-x11)/(y22-y11)/(y22-y11) ) / (   ( (x22-x11)/(y22-y11) )**2  + 1.    )
                !    b_trans = 2./( y22-y11 )  /  (   ( (x22-x11)/(y22-y11) )**2  + 1.    )
                !    c_trans = (  (x11-x22)/(y11-y22)*(x11+x22)/(y11-y22)  + (y11+y22)/(y11-y22)  ) /  (   ( (x22-x11)/(y22-y11) )**2  + 1.    )
                !    d_trans = (  (x11+x22)/(y11-y22)  +  (y11+y22)/(y22-y11)*(x22-x11)/(y22-y11)   )  /  (   ( (x22-x11)/(y22-y11) )**2  + 1.    )
                !else
                !    a_trans = 2./(x22-x11) / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
                !    b_trans = 2.*(y22-y11)/(x22-x11)/(x22-x11)   / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
                !    c_trans = (  (x11+x22)/(x11-x22) + (y11-y22)/(x11-x22)*(y11+y22)/(x11-x22)  )  / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
                !    d_trans = (  (x11+x22)/(x11-x22)*(y11-y22)/(x11-x22)  + (y11+y22)/(x22-x11)  )  / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
                !endif


                !********************* for identify trouble side ****************************
                !         if the side is trouble, we go back to choose a smaller time step
                !             or do it by straight segment.
                !****************************************************************************
                !
                ! three points, ( -1,0 ),( 1,0 ),( e_x2,e_y2 )

                e_x2 = trans_to_ex( a_trans ,b_trans , c_trans, x33,y33  )
                e_y2 = trans_to_ey( a_trans ,b_trans , d_trans, x33,y33  )
                !******************open direction of parabola

                e_yc = trans_to_ey( a_trans ,b_trans , d_trans, center(i,j,1),center(i,j,2) )
                if(  e_yc* e_y2/(e_x2**2-1. ) >0.   )then
                    if(kside==1)then
                        x33 = Dij_star(i,j-1)%vertex_star(2,1)
                        y33 = Dij_star(i,j-1)%vertex_star(2,2)
                        x44 = Dij_star(i,j-1)%vertex_star(1,1)
                        y44 = Dij_star(i,j-1)%vertex_star(1,2)
                    elseif(kside==2)then
                        x33 = Dij_star(i+1,j)%vertex_star(3,1)
                        y33 = Dij_star(i+1,j)%vertex_star(3,2)
                        x44 = Dij_star(i+1,j)%vertex_star(2,1)
                        y44 = Dij_star(i+1,j)%vertex_star(2,2)
                    elseif(kside==3)then
                        x33 = Dij_star(i,j+1)%vertex_star(4,1)
                        y33 = Dij_star(i,j+1)%vertex_star(4,2)
                        x44 = Dij_star(i,j+1)%vertex_star(3,1)
                        y44 = Dij_star(i,j+1)%vertex_star(3,2)
                    elseif(kside==4)then
                        x33 = Dij_star(i-1,j)%vertex_star(1,1)
                        y33 = Dij_star(i-1,j)%vertex_star(1,2)
                        x44 = Dij_star(i-1,j)%vertex_star(4,1)
                        y44 = Dij_star(i-1,j)%vertex_star(4,2)
                    endif
                else
                    if(kside==1)then
                        x33 = Dij_star(i,j)%vertex_star(3,1)
                        y33 = Dij_star(i,j)%vertex_star(3,2)
                        x44 = Dij_star(i,j)%vertex_star(4,1)
                        y44 = Dij_star(i,j)%vertex_star(4,2)
                    elseif(kside==2)then
                        x33 = Dij_star(i,j)%vertex_star(4,1)
                        y33 = Dij_star(i,j)%vertex_star(4,2)
                        x44 = Dij_star(i,j)%vertex_star(1,1)
                        y44 = Dij_star(i,j)%vertex_star(1,2)
                    elseif(kside==3)then
                        x33 = Dij_star(i,j)%vertex_star(1,1)
                        y33 = Dij_star(i,j)%vertex_star(1,2)
                        x44 = Dij_star(i,j)%vertex_star(2,1)
                        y44 = Dij_star(i,j)%vertex_star(2,2)
                    elseif(kside==4)then
                        x33 = Dij_star(i,j)%vertex_star(2,1)
                        y33 = Dij_star(i,j)%vertex_star(2,2)
                        x44 = Dij_star(i,j)%vertex_star(3,1)
                        y44 = Dij_star(i,j)%vertex_star(3,2)
                    endif
                endif
                call get_line( x11,y11,x33,y33,  at1,bt1,ct1 )
                call get_line( x22,y22,x44,y44,  at2,bt2,ct2 )

                atemp1 = at1*(x22-x11)/2. + bt1*(y22-y11)/2.
                btemp1 = at1*(y22-y11)/2. - bt1*(x22-x11)/2.
                atemp2 = at2*(x22-x11)/2. + bt2*(y22-y11)/2.
                btemp2 = at2*(y22-y11)/2. - bt2*(x22-x11)/2.

                if(  abs( (2.*e_y2)/(e_x2**2-1. ) )<=abs( atemp1/btemp1 )   &
                    .or.  abs( (2.*e_y2)/(e_x2**2-1. ) )<=abs( atemp2/btemp2 ) ) then
                !!!!!!!!!!!!!!!
                ! pause
                endif
                !
                !****************************************************************************

                if(    e_yc* e_y2/(e_x2**2-1. ) >0.   )then
                    !
                    if(kside==1)then
                        x00 = center(i,j-1,1)
                        y00 = center(i,j-1,2)

                        xleft_id = x_id_up(i,j-1,1)
                        xright_id = x_id_up(i,j-1,2)
                        yleft_id = y_id_up(i,j-1,1)
                        yright_id = y_id_up(i,j-1,2)
                    elseif(kside==2)then
                        x00 = center(i+1,j,1)
                        y00 = center(i+1,j,2)

                        xleft_id = x_id_up(i+1,j,1)
                        xright_id = x_id_up(i+1,j,2)
                        yleft_id = y_id_up(i+1,j,1)
                        yright_id = y_id_up(i+1,j,2)
                    elseif(kside==3)then
                        x00 = center(i,j+1,1)
                        y00 = center(i,j+1,2)

                        xleft_id = x_id_up(i,j+1,1)
                        xright_id = x_id_up(i,j+1,2)
                        yleft_id = y_id_up(i,j+1,1)
                        yright_id = y_id_up(i,j+1,2)
                    elseif(kside==4)then
                        x00 = center(i-1,j,1)
                        y00 = center(i-1,j,2)

                        xleft_id = x_id_up(i-1,j,1)
                        xright_id = x_id_up(i-1,j,2)
                        yleft_id = y_id_up(i-1,j,1)
                        yright_id = y_id_up(i-1,j,2)
                    endif
                    !
                else
                    x00 = center(i,j,1)
                    y00 = center(i,j,2)

                    xleft_id = x_id_up(i,j,1)
                    xright_id = x_id_up(i,j,2)
                    yleft_id = y_id_up(i,j,1)
                    yright_id = y_id_up(i,j,2)
                endif



                vleft = min( x00,Dij_star(i,j)%vertex_star( kside,1),Dij_star(i,j)%vertex_star( kside1,1) )
                vleft = min( vleft,Dij_star(i,j)%node_star( 2* kside,1) )
                vright = max( x00,Dij_star(i,j)%vertex_star( kside,1),Dij_star(i,j)%vertex_star( kside1,1) )
                vright = max( vright,Dij_star(i,j)%node_star( 2* kside,1) )

                vleft_y = min( y00,Dij_star(i,j)%vertex_star( kside,2),Dij_star(i,j)%vertex_star( kside1,2) )
                vleft_y = min( vleft_y , Dij_star(i,j)%node_star( 2*kside,2) )
                vright_y = max( y00,Dij_star(i,j)%vertex_star( kside,2),Dij_star(i,j)%vertex_star( kside1,2) )
                vright_y = max( vright_y,  Dij_star(i,j)%node_star( 2*kside,2) )

                call id_get17( vleft,dx,i,i1111,1 )
                ia = i1111 - isearch_extend;

                call id_get17( vright,dx,i,i1111,1 )
                ib = i1111 + isearch_extend

                call id_get17( vleft_y,dy,j,i1111,2 )
                ic = i1111 - isearch_extend

                call id_get17( vright_y,dy,j,i1111,2 )
                id = i1111 + isearch_extend

                mx = ib-ia
                my = id-ic

                vleft1 = min( Dij_star(i,j)%vertex_star( kside,1),Dij_star(i,j)%vertex_star( kside1,1) )
                vright1 = max( Dij_star(i,j)%vertex_star( kside,1),Dij_star(i,j)%vertex_star( kside1,1) )

                vleft_y1 = min( Dij_star(i,j)%vertex_star( kside,2),Dij_star(i,j)%vertex_star( kside1,2) )
                vright_y1 = max( Dij_star(i,j)%vertex_star( kside,2),Dij_star(i,j)%vertex_star( kside1,2) )

                call id_get17( vleft1,dx,i,ia1,1 );
                call id_get17( vright1,dx,i,ib1,1 );

                call id_get17( vleft_y1,dy,j,ic1,2 )
                call id_get17( vright_y1,dy,j,id1,2 )

                !******************************************************************************
                ! outter segment start from zz(1,:) = Dij_star(i,j)%vertex_star( kside, : )
                ! outter segment end at zz(kn+1,:) = Dij_star(i,j)%vertex_star( kside1, : )
                !******************************************************************************
                zz( 1, :) = Dij_star(i,j)%vertex_star( kside, : )

                kn = 1

                x33 = Dij_star(i,j)%node_star(2*kside,1)! for distance
                y33 = Dij_star(i,j)%node_star(2*kside,2)! for distance
                !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                ! for obtain "icurve"
                call get_line( x11,y11,x22,y22,  a16, b16,c16 )
                distance = abs( ( a16*x33+b16*y33+c16 )/( sqrt(a16**2+b16**2) ) )


                icurve_temp = 1
                !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


                if(mx .ne. 0)then
                    do k = 1 , mx
                        inter = ia + k
                        if(  abs(x22-x11) <= abs(y22-y11)  )then
                            Ac = e_y2/( e_x2**2 - 1.  )
                            Bc = (x22-x11)/(y22-y11)
                            Cc = -e_y2/( e_x2**2 - 1. ) - 2.* ( x0(inter) - (x11+x22)/2. ) / ( y22-y11 )


                            call get_quadratic_root( ac,bc,cc, e_xi(1:2),n_root )
                            if( n_root == 0 )then

                            elseif( n_root == 1 )then
                                if( abs( e_xi(1) )<= 1.+eps_2d8 )then
                                    zz(1+kn,1) = x0(inter)
                                    e_ta(1) = - (x22-x11)/(y22-y11) * e_xi(1) + 2.*( x0(inter)  - (x11+x22)/2. )/(y22-y11)


                                    zz(1+kn,2) = reverse_to_y( x11,y11,x22,y22,  e_xi(1),e_ta(1)   )
                                    if( inter <= ib1 .and. inter > ia1 )then
                                        ix = inter - iam
                                        if( isx(ix)==0  )then
                                            isx(ix) = 1
                                            point_inner_x_temp(ix,1,:) = zz(1+kn,:)
                                        elseif( isx(ix)==1 )then
                                            !
                                            point_inner_x_temp(ix,2,:) = zz(1+kn,:)
                                            isx(ix) = 2
                                        elseif( isx(ix) ==2 )then
                                            !
                                            point_inner_x_temp(ix,3,:) = zz(1+kn,:)
                                            isx(ix) = 3
                                        elseif( isx(ix) ==3 )then
                                            !
                                            point_inner_x_temp(ix,4,:) = zz(1+kn,:)
                                            isx(ix) = 4
                                        endif ! if( isx(ix)==0  )then

                                    endif !if( inter < ib .and. inter > ia )then
                                    kn = kn + 1

                                endif !if( abs(e_xi1)<= 1. )then
                            elseif( n_root == 2 )then
                                do ie = 1,2
                                    if(abs(e_xi(ie) )<=1. +eps_2d8 )then
                                        zz(1+kn,1) = x0(inter)
                                        e_ta(ie) = - (x22-x11)/(y22-y11)*e_xi(ie)  + 2.*( x0(inter) -(x11+x22)/2. ) / (y22-y11)
                                        zz(1+kn,2) = reverse_to_y( x11,y11,x22,y22,e_xi(ie),e_ta(ie) )

                                        kn = kn + 1
                                        !
                                    endif
                                enddo

                                ix = inter - iam


                                if(ix ==0)then
                                    if( abs(e_xi(1) )<=1.+eps_2d8 .and. abs(e_xi(2))<=1.+eps_2d8 )then
                                        !
                                        if(inter <= iam .or. inter > ibm )then
                                            inner_convex_x = inner_convex_x + 1
                                            point_inner_convex_x(inner_convex_x,1,:) =zz(kn-1,:)
                                            point_inner_convex_x(inner_convex_x,2,:) =zz(kn,:)
                                        else
                                            if(inner_concave_x(ix)==0 )then
                                                inner_concave_x(ix) = 1
                                                point_inner_concave_x(ix,1,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,2,:) = zz(kn,:)

                                            elseif(inner_concave_x(ix)==1 )then
                                                inner_concave_x(ix) = 2
                                                point_inner_concave_x(ix,3,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,4,:) = zz(kn,:)
                                            endif

                                        endif
                                    endif !e_xi
                                else
                                    !
                                    if( abs( e_xi(1) )<=1. +eps_2d8 .and. abs(e_xi(2) )>1. +eps_2d8  )then
                                        if(isx(ix)==0 )then
                                            isx(ix) = 1
                                            point_inner_x_temp(ix,1,:) = zz(kn,:)
                                        elseif( isx(ix)==1 )then
                                            !
                                            point_inner_x_temp(ix,2,:) = zz(kn,:)
                                            isx(ix) = 2
                                            !isx(ix)=0
                                        elseif(isx(ix)==2 )then
                                            !
                                            point_inner_x_temp(ix,3,:) = zz(kn,:)
                                            isx(ix) = 3
                                        elseif( isx(ix)==3 )then
                                            !
                                            point_inner_x_temp(ix,4,:) = zz(kn,:)
                                            isx(ix) = 4
                                        endif

                                    endif
                                    if( abs(e_xi(2) )<=1. +eps_2d8 .and. abs(e_xi(1) )>1. +eps_2d8 )then
                                        !
                                        if(isx(ix)==0)then
                                            isx(ix) = 1
                                            point_inner_x_temp(ix,1,:) = zz(kn,:)

                                        elseif( isx(ix)==1 )then
                                            !
                                            point_inner_x_temp(ix,2,:) = zz(kn,:)

                                            isx(ix) = 2
                                        elseif(isx(ix) ==2 )then
                                            !
                                            point_inner_x_temp(ix,3,:) = zz(kn,:)
                                            isx(ix) = 3
                                        elseif(isx(ix)==3 )then
                                            !
                                            point_inner_x_temp(ix,4,:) = zz(kn,:)
                                            isx(ix) = 4
                                            !isx(ix)=0
                                        endif
                                        !print *,2
                                        !pause
                                    endif
                                    if( abs(e_xi(1) )<=1.+eps_2d8 .and. abs(e_xi(2))<=1.+eps_2d8 )then
                                        !

                                        if(inter <= iam .or. inter > ibm )then
                                            inner_convex_x = inner_convex_x + 1
                                            point_inner_convex_x(inner_convex_x,1,:) =zz(kn-1,:)
                                            point_inner_convex_x(inner_convex_x,2,:) =zz(kn,:)
                                        else
                                            if(inner_concave_x(ix)==0 )then
                                                inner_concave_x(ix) = 1
                                                point_inner_concave_x(ix,1,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,2,:) = zz(kn,:)

                                            elseif(inner_concave_x(ix)==1 )then
                                                inner_concave_x(ix) = 2
                                                point_inner_concave_x(ix,3,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,4,:) = zz(kn,:)
                                            endif

                                        endif
                                    endif !e_xi
                                endif !ix
                            endif

                        endif !if(  abs(x22-x11) <= abs(y22-y11)  )then
                        if(  abs(y22-y11) < abs(x22-x11)  )then
                            !

                            Ac =  e_y2/(e_x2**2 - 1.) * ( (y22-y11)/(x22-x11) ) * ( (y22-y11)/(x22-x11) ) ;
                            Bc = -1. - 4.*( x0(inter)- (x11+x22)/2.  )/(x22-x11) * (y22-y11)/(x22-x11) *e_y2/( e_x2**2-1. )
                            Cc = e_y2/( e_x2**2 - 1. )*(  4.* ( (x0(inter) -(x11+x22)/2. )/(x22-x11) )**2  -1. )


                            call get_quadratic_root( ac,bc,cc, e_ta(1:2) ,n_root )
                            if( n_root == 0 )then

                            elseif( n_root ==1 )then
                                e_xi(1) = - (y22-y11)/(x22-x11)*e_ta(1) + 2.*( x0(inter) -(x11+x22)/2.  )/(x22-x11)
                                if( abs(e_xi(1) ) <=1. +eps_2d8 )then
                                    zz(1+kn,1) = x0(inter)
                                    zz(1+kn,2) = reverse_to_y(x11,y11,x22,y22, e_xi(1),e_ta(1) )
                                    if(inter<=ib1 .and. inter>ia )then
                                        ix = inter - iam
                                        if(isx(ix)==0 )then
                                            isx(ix) = 1
                                            point_inner_x_temp(ix,1,:) = zz(1+kn,:)
                                        elseif( isx(ix)==1 )then
                                            !
                                            point_inner_x_temp(ix,2,:) = zz(1+kn,:)
                                            isx(ix) = 2
                                        elseif( isx(ix)==2 )then
                                            !
                                            point_inner_x_temp(ix,3,:) = zz(1+kn,:)
                                            isx(ix) = 3
                                        elseif( isx(ix)==3 )then
                                            !
                                            point_inner_x_temp(ix,4,:) = zz(1+kn,:)
                                            isx(ix) = 4

                                        endif

                                    endif
                                    kn = kn + 1
                                endif
                            elseif( n_root ==2 )then
                                e_xi(1) = - (y22-y11)/(x22-x11)*e_ta(1) + 2.*( x0(inter) -(x11+x22)/2.  )/(x22-x11)
                                e_xi(2) = - (y22-y11)/(x22-x11)*e_ta(2) + 2.*( x0(inter) -(x11+x22)/2.  )/(x22-x11)

                                do ie = 1 ,2
                                    if( abs( e_xi(ie) )<= 1. +eps_2d8 )then
                                        zz(1+kn,1) = x0(inter)
                                        zz(1+kn,2) = reverse_to_y( x11,y11,x22,y22, e_xi(ie), e_ta(ie) )
                                        kn = kn + 1
                                    endif
                                enddo
                                ix = inter - iam
                                if(ix==0)then
                                    if( abs( e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) ) <=1.+eps_2d8 )then
                                        !
                                        if(inter <= iam .or. inter > ibm )then
                                            inner_convex_x = inner_convex_x + 1
                                            point_inner_convex_x(inner_convex_x,1,:) = zz(kn-1,:)
                                            point_inner_convex_x(inner_convex_x,2,:) = zz(kn,:)
                                        else
                                            if( inner_concave_x(ix) == 0 )then
                                                inner_concave_x(ix) = 1
                                                point_inner_concave_x(ix,1,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,2,:) = zz(kn,:)
                                            elseif( inner_concave_x(ix) == 1 )then
                                                inner_concave_x(ix) = 2
                                                point_inner_concave_x(ix,3,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,4,:) = zz(kn,:)
                                            endif ! inner_concave_x(ix)
                                        endif !inter
                                    endif ! e_xi
                                else
                                    !
                                    if( abs( e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) )>1. +eps_2d8 )then
                                        !
                                        if( isx(ix)==0 )then
                                            isx(ix) = 1
                                            point_inner_x_temp(ix,1,:) = zz( kn, : )

                                        elseif( isx(ix) == 1 )then
                                            !
                                            point_inner_x_temp(ix,2,:) = zz( kn,: )
                                            !
                                            isx(ix) = 2
                                        elseif( isx(ix)==2 )then
                                            !
                                            point_inner_x_temp(ix,3,:) = zz( kn,: )
                                            isx(ix) = 3
                                        elseif( isx(ix) ==3 )then
                                            !
                                            point_inner_x_temp(ix,4,:) = zz( kn,: )
                                            isx(ix) = 4
                                        endif
                                    endif
                                    if( abs( e_xi(2) )<=1.+eps_2d8 .and. abs( e_xi(1) )>1. +eps_2d8 )then
                                        if( isx(ix)==0 )then
                                            isx(ix) = 1
                                            point_inner_x_temp(ix,1,:) = zz( kn, : )
                                        elseif( isx(ix)==1 )then
                                            !
                                            point_inner_x_temp(ix,2,:) = zz( kn,: )
                                            !
                                            isx(ix) = 2
                                        elseif( isx(ix) == 2 )then
                                            !
                                            point_inner_x_temp(ix,3,:) = zz( kn,: )
                                            isx(ix) = 3
                                        elseif(isx(ix)==3 )then
                                            !
                                            point_inner_x_temp(ix,4,:) = zz( kn,: )
                                            isx(ix) = 4
                                            !isx(ix)=0
                                        endif
                                    endif
                                    if( abs( e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) ) <=1.+eps_2d8 )then
                                        !
                                        if(inter <= iam .or. inter > ibm )then
                                            inner_convex_x = inner_convex_x + 1
                                            point_inner_convex_x(inner_convex_x,1,:) = zz(kn-1,:)
                                            point_inner_convex_x(inner_convex_x,2,:) = zz(kn,:)
                                        else
                                            if( inner_concave_x(ix) == 0 )then
                                                inner_concave_x(ix) = 1
                                                point_inner_concave_x(ix,1,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,2,:) = zz(kn,:)
                                            elseif( inner_concave_x(ix) == 1 )then
                                                inner_concave_x(ix) = 2
                                                point_inner_concave_x(ix,3,:) = zz(kn-1,:)
                                                point_inner_concave_x(ix,4,:) = zz(kn,:)
                                            endif ! inner_concave_x(ix)
                                        endif !inter
                                    endif ! e_xi
                                endif!ix
                            endif

                        endif  ! if(  abs(y22-y11) < abs(x22-x11)  )then
                    enddo!do k = 1 , mx
                endif!               if(mx .ne. 0)then



                !********************************************************************************************************
                if(my .ne. 0 )then
                    do k = 1 , my
                        inter = ic + k


                        if( abs(x22-x11) <= abs(y22-y11) )then
                            !

                            Ac = e_y2/(e_x2**2 - 1. ) *  (  (x22-x11)/(y22-y11)   )**2
                            Bc = -1. + 4.*( y0(inter) - (y11+y22)/2. )/(y22-y11) *(x22-x11)/(y22-y11)*e_y2/(e_x2**2-1.)
                            Cc = (  4.*(  ( y0(inter) - (y11+y22)/2. )/(y22-y11)  )**2 - 1.  )* e_y2/(e_x2**2-1.)


                            call get_quadratic_root( ac,bc,cc, e_ta(1:2),n_root )

                            if( n_root == 0 )then

                            elseif( n_root == 1 )then
                                e_xi(1) = (x22-x11)/(y22-y11)*e_ta(1) + 2.*( y0(inter)  -(y11+y22)/2. )/( y22-y11 )
                                if(  abs(e_xi(1) ) <= 1. +eps_2d8  )then
                                    zz(1+kn,2) = y0(inter)
                                    zz(1+kn,1) = reverse_to_x( x11,y11,x22,y22,e_xi(1),e_ta(1) )
                                    if( inter<= id1 .and. inter>ic1 )then
                                        iy = inter - icm
                                        if(isy(iy)==0 )then
                                            isy(iy) = 1
                                            point_inner_y_temp(iy,1,:) = zz(1+kn,:)
                                        elseif( isy(iy)==1 )then
                                            point_inner_y_temp(iy,2,:) = zz(1+kn,:)
                                            !
                                            isy(iy) = 2
                                        elseif( isy(iy) ==2 )then
                                            !
                                            point_inner_y_temp(iy,3,:) = zz(1+kn,:)
                                            isy(iy) = 3
                                        elseif( isy(iy)==3 )then
                                            !
                                            point_inner_y_temp(iy,4,:) = zz(1+kn,:)
                                            isy(iy) = 4
                                        endif !isy
                                    endif!inter
                                    kn = kn + 1
                                endif!exi
                            elseif( n_root == 2 )then
                                e_xi(1) = (x22-x11)/(y22-y11) * e_ta(1) + 2.*( y0(inter) -(y11+y22)/2. )/(y22-y11)
                                e_xi(2) = (x22-x11)/(y22-y11) * e_ta(2) + 2.*( y0(inter) -(y11+y22)/2. )/(y22-y11)

                                do ie = 1,2
                                    if( abs(e_xi(ie) ) <=1. +eps_2d8 )then
                                        zz(1+kn,2) = y0(inter)
                                        zz(1+kn,1) = reverse_to_x(x11,y11,x22,y22,e_xi(ie),e_ta(ie) )

                                        kn = kn+1
                                    endif
                                enddo

                                iy = inter -icm
                                if(iy==0)then
                                    if( abs(e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) )<= 1.+eps_2d8  )then
                                        !
                                        if(inter <= icm .or. inter >idm)then
                                            !
                                            inner_convex_y = inner_convex_y + 1
                                            point_inner_convex_y( inner_convex_y,1,: ) = zz(kn-1,:)
                                            point_inner_convex_y( inner_convex_y,2,: ) = zz(kn,:)
                                        else
                                            iy = inter - icm
                                            if(inner_concave_y(iy)==0 )then
                                                inner_concave_y(iy) = 1
                                                point_inner_concave_y( iy,1,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,2,: ) = zz(kn,:)
                                            elseif( inner_concave_y(iy) == 1 )then
                                                inner_concave_y(iy) = 2
                                                point_inner_concave_y( iy,3,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,4,: ) = zz(kn,:)
                                            endif
                                        endif! inter

                                    endif ! abs(e_xi)
                                else
                                    !
                                    if( abs(e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) ) > 1. +eps_2d8  )then
                                        if( isy(iy)==0 )then
                                            isy(iy) = 1
                                            point_inner_y_temp(iy,1,:) = zz(kn,:)
                                        elseif( isy(iy) == 1 )then
                                            point_inner_y_temp(iy,2,:) = zz(kn,:)
                                            !
                                            isy(iy) = 2
                                        elseif( isy(iy) == 2 )then
                                            point_inner_y_temp(iy,3,:) = zz(kn,:)
                                            isy(iy) = 3
                                        elseif( isy(iy) == 3 )then
                                            point_inner_y_temp(iy,4,:) = zz(kn,:)
                                            isy(iy) = 4
                                        endif
                                    endif
                                    if( abs( e_xi(2) ) <=1.+eps_2d8 .and. abs( e_xi(1) ) > 1. +eps_2d8 )then
                                        if(isy(iy)==0 )then
                                            isy(iy) = 1
                                            point_inner_y_temp(iy,1,:) = zz(kn,:)
                                        elseif(isy(iy)==1)then
                                            !
                                            point_inner_y_temp(iy,2,:) = zz(kn,:)
                                            isy(iy) = 2
                                        elseif(isy(iy)==2)then
                                            !
                                            point_inner_y_temp(iy,3,:) = zz(kn,:)
                                            isy(iy) = 3
                                        elseif( isy(iy)==3 )then
                                            !
                                            point_inner_y_temp(iy,4,:) = zz(kn,:)
                                            isy(iy) = 4
                                        endif

                                    endif !exi2
                                    !
                                    if( abs(e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) )<= 1.+eps_2d8  )then
                                        !
                                        if(inter <= icm .or. inter >idm)then
                                            !
                                            inner_convex_y = inner_convex_y + 1
                                            point_inner_convex_y( inner_convex_y,1,: ) = zz(kn-1,:)
                                            point_inner_convex_y( inner_convex_y,2,: ) = zz(kn,:)
                                        else
                                            iy = inter - icm
                                            if(inner_concave_y(iy)==0 )then
                                                inner_concave_y(iy) = 1
                                                point_inner_concave_y( iy,1,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,2,: ) = zz(kn,:)
                                            elseif( inner_concave_y(iy) == 1 )then
                                                inner_concave_y(iy) = 2
                                                point_inner_concave_y( iy,3,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,4,: ) = zz(kn,:)
                                            endif
                                        endif! inter
                                    endif ! abs(e_xi)
                                endif ! iy==0

                            endif
                        endif ! |x2-x1|  <  |y2-y1|
                        !*************************************************************************************************************
                        if( abs( y22-y11 ) <abs(x22-x11)  )then
                            Ac = e_y2/( e_x2**2 - 1. )
                            Bc = - (y22-y11)/(x22-x11)
                            Cc = - e_y2/( e_x2**2 - 1. ) - 2.*(  (y11+y22)/2. - y0(inter)   )/( x22 - x11 )

                            call get_quadratic_root(  ac,bc,cc, e_xi(1:2) ,n_root )

                            if( n_root == 0 )then

                            elseif( n_root == 1 )then
                                e_ta(1) = (y22- y11)/(x22-x11)*e_xi(1) + 2.* ( (y11+y22)/2. - y0(inter)  )/( x22 -x11 )
                                if(  abs( e_xi(1) ) <= 1. +eps_2d8)then
                                    zz(1+kn, 2) = y0(inter)
                                    zz(1+kn, 1) = reverse_to_x( x11,y11,x22,y22, e_xi(1) ,e_ta(1) )
                                    if( inter<= id1 .and. inter >ic1 )then
                                        iy = inter -icm
                                        if(isy(iy) ==0 )then
                                            isy(iy) = 1
                                            point_inner_y_temp( iy,1,: ) = zz(1+kn, :)
                                        elseif( isy(iy) == 1 )then
                                            !
                                            point_inner_y_temp( iy,2,: ) = zz(1+kn, :)
                                            !
                                            isy(iy) = 2
                                        elseif( isy(iy) == 2 )then
                                            point_inner_y_temp( iy,3,: ) = zz(1+kn,: )
                                            isy(iy) = 3
                                        elseif( isy(iy) ==3 )then
                                            point_inner_y_temp( iy,4,: ) = zz(1+kn,:)
                                            isy(iy) = 4
                                        endif

                                    endif!inter
                                    kn = kn + 1
                                endif!exi1
                            elseif( n_root == 2 )then

                                e_ta(1) = (y22- y11)/(x22-x11)*e_xi(1) + 2.* ( (y11+y22)/2. - y0(inter) )/( x22 -x11 )
                                e_ta(2) = (y22- y11)/(x22-x11)*e_xi(2) + 2.* ( (y11+y22)/2. - y0(inter) )/( x22 -x11 )

                                do ie = 1,2
                                    if( abs( e_xi(ie) ) <= 1. +eps_2d8 )then
                                        zz(1+kn,2) = y0(inter)
                                        zz(1+kn,1) = reverse_to_x( x11,y11,x22,y22,e_xi(ie) , e_ta(ie) )
                                        kn = kn + 1
                                    endif

                                enddo
                                iy = inter -icm

                                if(iy == 0)then

                                    if( abs(e_xi(1) )<=1. +eps_2d8 .and. abs( e_xi(2) )<= 1.+eps_2d8  )then
                                        !
                                        if(inter <= icm .or. inter >idm)then
                                            inner_convex_y = inner_convex_y + 1
                                            point_inner_convex_y( inner_convex_y,1,: ) = zz(kn-1,:)
                                            point_inner_convex_y( inner_convex_y,2,: ) = zz(kn,:)
                                        else
                                            iy = inter - icm
                                            if(inner_concave_y(iy)==0 )then
                                                inner_concave_y(iy) = 1
                                                point_inner_concave_y( iy,1,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,2,: ) = zz(kn,:)
                                            elseif( inner_concave_y(iy) == 1 )then
                                                inner_concave_y(iy) = 2
                                                point_inner_concave_y( iy,3,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,4,: ) = zz(kn,:)
                                            endif

                                        endif!inter
                                    endif ! abs(e_xi)

                                else
                                    !
                                    if( abs(e_xi(1) )<=1.+eps_2d8 .and. abs( e_xi(2) ) > 1. +eps_2d8  )then
                                        if( isy(iy)==0 )then
                                            isy(iy) = 1
                                            point_inner_y_temp(iy,1,:) = zz(kn,:)
                                        elseif( isy(iy) == 1 )then
                                            !
                                            point_inner_y_temp(iy,2,:) = zz(kn,:)
                                            !
                                            isy(iy) = 2
                                        elseif( isy(iy) ==2 )then
                                            point_inner_y_temp(iy,3,:) = zz(kn,:)
                                            isy(iy) = 3
                                        elseif(isy(iy)==3 )then
                                            point_inner_y_temp(iy,4,:) = zz(kn,:)
                                            isy(iy) = 4
                                        endif

                                    endif!exi1
                                    if( abs( e_xi(2) ) <=1.+eps_2d8 .and. abs( e_xi(1) ) > 1. +eps_2d8 )then
                                        if(isy(iy)==0 )then                    !!!!!!!!!!!!!!!!!
                                            isy(iy) = 1
                                            point_inner_y_temp(iy,1,:) = zz(kn,:)
                                        elseif(isy(iy)==1)then
                                            point_inner_y_temp(iy,2,:) = zz(kn,:)
                                            isy(iy)=2
                                        elseif(isy(iy)==2 )then
                                            point_inner_y_temp(iy,3,:) = zz(kn,:)
                                            isy(iy) = 3
                                        elseif( isy(iy)==3 )then
                                            point_inner_y_temp(iy,4,:) = zz(kn,:)
                                            isy(iy) = 4
                                        endif

                                    endif
                                    !
                                    if( abs(e_xi(1) )<=1. +eps_2d8 .and. abs( e_xi(2) )<= 1.+eps_2d8  )then
                                        !
                                        if(inter <= icm .or. inter >idm)then
                                            inner_convex_y = inner_convex_y + 1
                                            point_inner_convex_y( inner_convex_y,1,: ) = zz(kn-1,:)
                                            point_inner_convex_y( inner_convex_y,2,: ) = zz(kn,:)
                                        else
                                            iy = inter - icm
                                            if(inner_concave_y(iy)==0 )then
                                                inner_concave_y(iy) = 1
                                                point_inner_concave_y( iy,1,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,2,: ) = zz(kn,:)
                                            elseif( inner_concave_y(iy) == 1 )then
                                                inner_concave_y(iy) = 2
                                                point_inner_concave_y( iy,3,: ) = zz(kn-1,:)
                                                point_inner_concave_y( iy,4,: ) = zz(kn,:)
                                            endif

                                        endif!inter
                                    endif ! abs(e_xi)
                                endif !iy

                            endif


                        endif !abs( y22-y11 ) <abs(x22-x11)
                    enddo  ! my
                endif ! my .ne. 0
                !*****************************
                zz( kn+1,: ) = Dij_star( i,j ) % vertex_star( kside1,: )

                if( kn > 2 )then
                    do ii = 1 , kn - 2
                        do jj = 2 , kn -ii
                            compare1 = trans_to_ex( a_trans,b_trans,c_trans, zz(jj,1),zz(jj,2) )
                            compare2 = trans_to_ex( a_trans,b_trans,c_trans, zz(jj+1,1),zz(jj+1,2) )
                            if( compare1 > compare2 )then
                                zztemp( : ) = zz( jj+1 , : )
                                zz( jj+1,: ) = zz( jj,: )
                                zz( jj,: ) = zztemp( : )
                            endif
                        enddo
                    enddo

                endif! kn>2
                !  get final outer segments!!!
                do ii = Dij_star( i,j )% nsub_outer + 1, Dij_star( i,j )% nsub_outer + kn
                    jj = ii - Dij_star(i,j)%nsub_outer
                    Dij_star(i,j)%segment_outer(ii)%vl(:) = zz(jj,:)
                    Dij_star(i,j)%segment_outer(ii)%vr(:) = zz(jj+1,:)


                    if(icurve_temp == 0)then
                        xid = ( zz(jj,1)+zz(jj+1,1) )/2.
                        yid = ( zz(jj,2)+zz(jj+1,2) )/2.
                        Dij_star(i,j)%segment_outer(ii)%icurve = 0
                    else
                        xid = ( zz(jj,1)+zz(jj+1,1) )/2.
                        yid = ( zz(jj,2)+zz(jj+1,2) )/2.

                        e_temp1 = trans_to_ex( a_trans,b_trans,c_trans, zz(jj,1),zz(jj,2)  )
                        e_temp2 = trans_to_ex( a_trans,b_trans,c_trans, zz(jj+1,1),zz(jj+1,2) )
                        e_tempxc = ( e_temp1+e_temp2 )/2.
                        e_tempyc = e_y2/( e_x2**2 -1. ) *( e_tempxc**2 - 1. )
                        xid = reverse_to_x( x11,y11,x22,y22, e_tempxc, e_tempyc )
                        yid = reverse_to_y( x11,y11,x22,y22, e_tempxc, e_tempyc )
                        Dij_star(i,j)%segment_outer(ii)%icurve = 1
                    endif

                    call id_get17( xid,dx,i,idx,1 )
                    call id_get17( yid,dy,j,idy,2 )


                    Dij_star( i,j )%segment_outer(ii) %id(1) = idx
                    Dij_star( i,j )%segment_outer(ii) %id(2) = idy
                    Dij_star( i,j )%segment_outer(ii) %v11(:) = Dij_star(i,j)%vertex_star(kside,:)
                    Dij_star( i,j )%segment_outer(ii) %v22(:) = Dij_star(i,j)%vertex_star(kside1,:)
                    Dij_star( i,j )%segment_outer(ii) %v33(:) = Dij_star(i,j)%node_star(kside*2,:)
                enddo

                Dij_star(i,j)%nsub_outer = Dij_star( i,j ) % nsub_outer + kn

            enddo !do kside = 1,4


            do ix = 1 , (ibm-iam)
                !
                if( isx(ix) == 1 )then

                    point_inner_x(ix,1,:) = point_inner_x_temp(ix,1,:)
                    point_inner_x(ix,2,:) = point_inner_x_temp(ix,1,:)
                elseif(isx(ix) == 2 )then
                    point_inner_x(ix,1,:) = point_inner_x_temp(ix,1,:)
                    point_inner_x(ix,2,:) = point_inner_x_temp(ix,2,:)
                elseif(isx(ix)==3 )then
                    distance1 = distance_f(  point_inner_x_temp(ix,1,:), point_inner_x_temp(ix,2,:) )
                    distance2 = distance_f(  point_inner_x_temp(ix,1,:), point_inner_x_temp(ix,3,:) )
                    if( distance1 > distance2 )then
                        point_inner_x(ix,1,:) = point_inner_x_temp(ix,1,:)
                        point_inner_x(ix,2,:) = point_inner_x_temp(ix,2,:)
                    else
                        point_inner_x(ix,1,:) = point_inner_x_temp(ix,1,:)
                        point_inner_x(ix,2,:) = point_inner_x_temp(ix,3,:)
                    endif
                elseif(isx(ix)==4)then
                    distance1 = distance_f(  point_inner_x_temp(ix,1,:), point_inner_x_temp(ix,2,:) )
                    distance2 = distance_f(  point_inner_x_temp(ix,1,:), point_inner_x_temp(ix,3,:) )
                    if( distance1 > distance2 )then
                        point_inner_x(ix,1,:) = point_inner_x_temp(ix,1,:)
                        point_inner_x(ix,2,:) = point_inner_x_temp(ix,2,:)
                    else
                        point_inner_x(ix,1,:) = point_inner_x_temp(ix,1,:)
                        point_inner_x(ix,2,:) = point_inner_x_temp(ix,3,:)
                    endif
                endif

                if( inner_concave_x(ix) == 0 )then
                    temp(1,:) = point_inner_x( ix,1,: )
                    temp(2,:) = point_inner_x( ix,2,: )
                elseif( inner_concave_x(ix) == 1 )then
                    temp(1,:) = point_inner_x( ix,1,: )
                    temp(2,:) = point_inner_concave_x( ix,1,: )
                    temp(3,:) = point_inner_concave_x( ix,2,: )
                    temp(4,:) = point_inner_x( ix,2,: )
                    do ii17 = 1 , 3
                        do jj17 = 1, 4- ii17
                            if( temp(jj17,2) > temp(jj17+1,2) )then
                                zztemp(:) = temp(jj17+1,:)
                                temp(jj17+1,:) = temp(jj17,:)
                                temp(jj17,:) = zztemp(:)
                            endif
                        enddo
                    enddo
                elseif( inner_concave_x(ix) == 2 )then
                    temp( 1,: ) = point_inner_x( ix,1,: )
                    temp( 2,: ) = point_inner_concave_x(ix,1,:)
                    temp( 3,: ) = point_inner_concave_x(ix,2,:)
                    temp( 4,: ) = point_inner_concave_x(ix,3,:)
                    temp( 5,: ) = point_inner_concave_x(ix,4,:)
                    temp( 6,: ) = point_inner_x( ix,2,: )
                    do ii17 = 1 , 5
                        do jj17 = 1, 6- ii17
                            if( temp(jj17,2) > temp(jj17+1,2) )then
                                zztemp(:) = temp(jj17+1,:)
                                temp(jj17+1,:) = temp(jj17,:)
                                temp(jj17,:) = zztemp(:)
                            endif

                        enddo
                    enddo
                endif

                do ix_concave = 1, inner_concave_x(ix) + 1
                    call id_get17( temp(2*ix_concave-1,2), dy, j,i1111, 2 )
                    call id_get17( temp(2*ix_concave,2),   dy, j,i2222, 2 )

                    nby = i1111  -  i2222


                    temp_store( : ) = temp( 2*ix_concave, : )
                    point_segment(1,:) = temp( 2*ix_concave-1 , : )

                    call id_get17( temp(2*ix_concave-1,2), dy , j,iy0,2 )

                    if( nby>0 )then
                        do ii = 1 , nby
                            point_segment(ii+1,1) = temp( 2*ix_concave-1,1 )
                            point_segment(ii+1,2) = y0(iy0+1-ii)

                        enddo

                    elseif( nby<0)then
                        do ii = 1 , abs(nby)
                            point_segment(ii+1,1 ) = temp( 2*ix_concave-1,1 )
                            point_segment(ii+1,2) = y0(iy0+ ii)
                        enddo

                    else
                    endif
                    point_segment(2+abs(nby),: ) = temp_store(:)

                    do ii = 1,1+abs(nby)
                        jj = Dij_star(i,j)% nsub_inner + ii

                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii+1,:)


                        call id_get17( ( point_segment(ii,2)+point_segment(ii+1,2) )/2., dy,j,idy,2 )

                        Dij_star(i,j)% segment_inner(jj)%id(2) = idy
                        if( point_segment(ii,2)>point_segment(ii+1,2)  )then

                            call id_get_up_down( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        else

                            call id_get_down_up( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx

                        endif

                    enddo !ii

                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nby) + 1

                    do ii = 1 , 1 + abs(nby)
                        jj = Dij_star(i,j)% nsub_inner + ii
                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii+1,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii,:)

                        call id_get17( ( point_segment(ii,2)+point_segment(ii+1,2) )/2., dy,j,idy,2 )

                        Dij_star(i,j)% segment_inner(jj)%id(2) = idy
                        if( point_segment(ii+1,2)>point_segment(ii,2)  )then

                            call id_get_up_down( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        else

                            call id_get_down_up( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        endif
                        !

                    enddo!ii

                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nby) + 1
                enddo!ix_concave = 1, inner_concave_x(ix) + 1

            enddo ! ix
            !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            do iy = 1 ,(idm-icm)
                !
                if(isy(iy)==1)then
                    !print *,i,j,'iy=',iy
                    !pause
                    point_inner_y(iy,1,:) = point_inner_y_temp(iy,1,:)
                    point_inner_y(iy,2,:) = point_inner_y_temp(iy,1,:)
                elseif(isy(iy) == 2 )then
                    point_inner_y(iy,1,:) = point_inner_y_temp(iy,1,:)
                    point_inner_y(iy,2,:) = point_inner_y_temp(iy,2,:)
                elseif(isy(iy)==3 )then
                    distance1 = distance_f(  point_inner_y_temp(iy,1,:), point_inner_y_temp(iy,2,:) )
                    distance2 = distance_f(  point_inner_y_temp(iy,1,:), point_inner_y_temp(iy,3,:) )
                    if( distance1 > distance2 )then
                        point_inner_y(iy,1,:) = point_inner_y_temp(iy,1,:)
                        point_inner_y(iy,2,:) = point_inner_y_temp(iy,2,:)
                    else
                        point_inner_y(iy,1,:) = point_inner_y_temp(iy,1,:)
                        point_inner_y(iy,2,:) = point_inner_y_temp(iy,3,:)
                    endif
                elseif(isy(iy)==4)then
                    distance1 = distance_f(  point_inner_y_temp(iy,1,:), point_inner_y_temp(iy,2,:) )
                    distance2 = distance_f(  point_inner_y_temp(iy,1,:), point_inner_y_temp(iy,3,:) )
                    if( distance1 > distance2 )then
                        point_inner_y(iy,1,:) = point_inner_y_temp(iy,1,:)
                        point_inner_y(iy,2,:) = point_inner_y_temp(iy,2,:)
                    else
                        point_inner_y(iy,1,:) = point_inner_y_temp(iy,1,:)
                        point_inner_y(iy,2,:) = point_inner_y_temp(iy,3,:)
                    endif
                endif


                if( inner_concave_y(iy) == 0 )then
                    temp(1,:) = point_inner_y( iy,1,: )
                    temp(2,:) = point_inner_y( iy,2,: )
                elseif( inner_concave_y(iy) == 1 )then
                    temp(1,:) = point_inner_y( iy,1,: )
                    temp(2,:) = point_inner_concave_y( iy,1,: )
                    temp(3,:) = point_inner_concave_y( iy,2,: )
                    temp(4,:) = point_inner_y( iy,2,: )

                    do ii17 = 1 , 3
                        do jj17 = 1, 4- ii17
                            if( temp(jj17,1) > temp(jj17+1,1) )then
                                zztemp(:) = temp(jj17+1,:)
                                temp(jj17+1,:) = temp(jj17,:)
                                temp(jj17,:) = zztemp(:)
                            endif

                        enddo
                    enddo
                elseif( inner_concave_y(iy) == 2 )then
                    temp( 1,: ) = point_inner_y( iy,1,: )
                    temp( 2,: ) = point_inner_concave_y(iy,1,:)
                    temp( 3,: ) = point_inner_concave_y(iy,2,:)
                    temp( 4,: ) = point_inner_concave_y(iy,3,:)
                    temp( 5,: ) = point_inner_concave_y(iy,4,:)
                    temp( 6,: ) = point_inner_y( iy,2,: )

                    do ii17 = 1 , 5
                        do jj17 = 1, 6- ii17
                            if( temp(jj17,1) > temp(jj17+1,1) )then
                                zztemp(:) = temp(jj17+1,:)
                                temp(jj17+1,:) = temp(jj17,:)
                                temp(jj17,:) = zztemp(:)
                            endif

                        enddo
                    enddo
                endif

                do iy_concave = 1, inner_concave_y(iy) + 1
                    call id_get17( temp(2*iy_concave-1,1),dx,i,i1111, 1 )
                    call id_get17( temp(2*iy_concave,1),dx,i,i2222, 1 )
                    nbx = i1111 - i2222

                    temp_store( : ) = temp( 2*iy_concave, : )
                    point_segment(1,:) = temp( 2*iy_concave-1 , : )

                    call id_get17( temp(2*iy_concave-1,1),dx,i,ix0,1 )

                    if( nbx>0 )then
                        do ii = 1 , nbx
                            point_segment(ii+1,2) = temp( 2*iy_concave-1,2 )
                            point_segment(ii+1,1) = x0(ix0+1-ii)
                        enddo
                    elseif( nbx<0)then
                        do ii = 1 , abs(nbx)
                            point_segment(ii+1,2 ) = temp( 2*iy_concave-1,2 )
                            point_segment(ii+1,1) = x0(ix0+ ii)
                        enddo
                    else
                    endif
                    point_segment(2+abs(nbx),: ) = temp_store(:)


                    do ii = 1,1+abs(nbx)
                        jj = Dij_star(i,j)% nsub_inner + ii

                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii+1,:)

                        call id_get17( ( point_segment(ii,1)+point_segment(ii+1,1) )/2.,dx,i,idx,1 )

                        Dij_star(i,j)% segment_inner(jj)%id(1) = idx
                        if( point_segment(ii,1)< point_segment(ii+1,1)  )then

                            call id_get_left_right( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        else

                            call id_get_right_left( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        endif
                    enddo
                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nbx) + 1

                    do ii = 1 , 1 + abs(nbx)
                        jj = Dij_star(i,j)% nsub_inner + ii
                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii+1,:)


                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii,:)

                        call id_get17( ( point_segment(ii,1)+point_segment(ii+1,1) )/2.,dx,i,idx,1 )

                        Dij_star(i,j)% segment_inner(jj)%id(1) = idx
                        if( point_segment(ii+1,1)< point_segment(ii,1)  )then

                            call id_get_left_right( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        else


                            call id_get_right_left( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        endif
                    enddo
                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nbx) + 1
                enddo!
                !endif!if(isy(iy) == 1 )then
            enddo! iy

            If( inner_convex_x >0 )then
                do ic = 1 , inner_convex_x
                    call id_get17( point_inner_convex_x(ic,1,2),dy,j,i1111,2 )
                    call id_get17( point_inner_convex_x(ic,2,2),dy,j,i2222,2 )

                    nby = i1111 - i2222

                    temp_store(:) = point_inner_convex_x(ic,2,:)
                    point_segment( 1,: ) = point_inner_convex_x(ic,1,:)

                    call id_get17( point_segment(1,2),dy,j,iy0,2 )

                    if( nby>0 )then
                        do ii = 1 , nby
                            point_segment(ii+1,1) = point_segment( 1,1 )
                            point_segment(ii+1,2) = y0(iy0+1-ii)
                        enddo
                    elseif( nby<0)then
                        do ii = 1 , abs(nby)
                            point_segment(ii+1,1 ) = point_segment( 1,1 )
                            point_segment(ii+1,2) = y0(iy0+ ii)
                        enddo
                    else
                    endif
                    point_segment(2+abs(nby),: ) = temp_store(:)

                    do ii = 1,1+abs(nby)
                        jj = Dij_star(i,j)% nsub_inner + ii
                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii+1,:)
                        call id_get17( ( point_segment(ii,2)+point_segment(ii+1,2) )/2.,dy,j,idy,2  )


                        Dij_star(i,j)% segment_inner(jj)%id(2) = idy
                        if( point_segment(ii,2)>point_segment(ii+1,2)  )then


                            call id_get_up_down( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        else


                            call id_get_down_up( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        endif
                    enddo
                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nby) + 1
                    do ii = 1 , 1 + abs(nby)
                        jj = Dij_star(i,j)% nsub_inner + ii
                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii+1,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii,:)
                        call id_get17(  ( point_segment(ii,2)+point_segment(ii+1,2) )/2., dy ,j,idy,2 )

                        Dij_star(i,j)% segment_inner(jj)%id(2) = idy
                        if( point_segment(ii+1,2)>point_segment(ii,2)  )then


                            call id_get_up_down( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        else


                            call id_get_down_up( point_segment(ii,1),dx,i,idx )
                            Dij_star(i,j)% segment_inner(jj)% id(1) = idx
                        endif
                    enddo
                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nby) + 1
                enddo
            endif


            if( inner_convex_y >0 )then
                do ic = 1 , inner_convex_y
                    call id_get17( point_inner_convex_y(ic,1,1),dx,i,i1111, 1 )
                    call id_get17( point_inner_convex_y(ic,2,1),dx,i,i2222, 1 )
                    nbx = i1111 - i2222

                    temp_store(:) = point_inner_convex_y( ic,2,: )
                    point_segment( 1,: ) = point_inner_convex_y(ic,1,:)

                    call id_get17( point_segment(1,1),dx,i,ix0,1 )

                    if( nbx>0 )then
                        do ii = 1 , nbx
                            point_segment(ii+1,2) = point_segment( 1,2 )
                            point_segment(ii+1,1) = x0(ix0+1-ii)
                        enddo

                    elseif( nbx<0)then
                        do ii = 1 , abs(nbx)
                            point_segment(ii+1,2 ) = point_segment( 1,2 )
                            point_segment(ii+1,1) = x0(ix0+ ii)
                        enddo
                    else
                    endif
                    point_segment(2+abs(nbx),: ) = temp_store(:)

                    do ii = 1,1+abs(nbx)
                        jj = Dij_star(i,j)% nsub_inner + ii
                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii+1,:)

                        call id_get17( ( point_segment(ii,1)+point_segment(ii+1,1) )/2.,dx,i,idx,1 )

                        Dij_star(i,j)% segment_inner(jj)%id(1) = idx
                        if( point_segment(ii,1)< point_segment(ii+1,1)  )then


                            call id_get_left_right( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        else


                            call id_get_right_left( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        endif

                    enddo
                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nbx) + 1
                    do ii = 1 , 1 + abs(nbx)
                        jj = Dij_star(i,j)% nsub_inner + ii
                        Dij_star(i,j) % segment_inner(jj)% vl(:) = point_segment(ii+1,:)
                        Dij_star(i,j) % segment_inner(jj)% vr(:) = point_segment(ii,:)

                        call id_get17( ( point_segment(ii,1)+point_segment(ii+1,1) )/2., dx,i,idx,1 )


                        Dij_star(i,j)% segment_inner(jj)%id(1) = idx
                        if( point_segment(ii+1,1)< point_segment(ii,1)  )then


                            call id_get_left_right( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        else

                            call id_get_right_left( point_segment(ii,2),dy,j,idy )
                            Dij_star(i,j)% segment_inner(jj)% id(2) = idy
                        endif
                    enddo
                    Dij_star(i,j)%nsub_inner = Dij_star(i,j)% nsub_inner + abs(nbx) + 1
                enddo
            endif


            call green_p2QC_gauss3_eulerian(umod_temp(i,j,1:6) )

        enddo !do j = 1, ny
    enddo  !do i = 1, nx



    do i = 1 ,nx
        do j = 1,ny
            Dij(i,j)%umodal(1:n_moment) =  umod_temp(i,j,1:n_moment)
        enddo
    enddo


    end subroutine search_QC_formula


    real function distance_f( a17,b17 )
    implicit none
    real,intent(in) :: a17(2),b17(2)

    distance_f = sqrt( ( a17(2)-b17(2) )**2 + ( a17(1)-b17(1) )**2 )


    end function distance_f

    !
    real function yside( v11, v22, x)
    implicit none
    real,intent(in) :: v11(2),v22(2)
    real,intent(in) :: x

    yside = ( v22(2)-v11(2) )/(v22(1)-v11(1) ) * ( x- v11(1) ) + v11(2)


    end function yside
    !
    real function xside( v11, v22, y)
    implicit none
    real,intent(in) :: v11(2),v22(2)
    real,intent(in) :: y

    xside = ( v22(1)-v11(1) )/(v22(2)-v11(2) ) * ( y- v11(2) ) + v11(1)


    end function xside

