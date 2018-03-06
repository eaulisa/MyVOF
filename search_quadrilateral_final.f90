    subroutine search_quadrilateral_final
    implicit none
    real :: xmax,xmin,ymax,ymin
    integer :: inter
    integer :: im
    integer :: i1111,i2222
    integer :: isx(18),isy(18)
    integer :: iam,ibm,icm,idm
    integer :: kside,kside1
    real :: vleft, vright,vleft_y,vright_y
    integer :: ia,ib,mx
    integer :: ic,id,my

    real :: zz(18,2)
    real :: point_inner_x( 18,18,2 ),point_inner_y( 18,18,2 )
    integer :: ix,iy
    integer :: ii,jj
    real :: temp(2)

    integer :: idx,idy,nby,nbx,ix0,iy0

    do i = 1 , nx
        do j = 1,ny
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

            !iam = id_get( (xmin-xleft)/dx  ); ibm = id_get( (xmax-xleft)/dx );


            isx(:) = 0

            call id_get17( ymin,dy,j,icm,2 )
            call id_get17( ymax,dy,j,idm,2 )

            !icm = id_get( (ymin-yleft)/dy );  idm = id_get( (ymax-yleft)/dy );

            isy(:) = 0


            Dij_star(i,j)%nsub_outer = 0
            Dij_star(i,j)%nsub_inner = 0

            do kside = 1,4
                ! connect sides: side1, 1--2, side2, 2--3, side3, 3--4, side4, 4--1
                kside1 = kside + 1
                if(kside ==4 ) kside1 = 1
                !
                ! Step 1: find intersections between side_k and x=x_i
                vleft = Dij_star(i,j)%vertex_star( kside,1)
                vright = Dij_star(i,j)%vertex_star( kside1,1)

                call id_get17( vleft,dx,i,ia,1 )
                call id_get17( vright,dx,i,ib,1 )


                vleft_y = Dij_star(i,j)%vertex_star( kside,2)
                vright_y = Dij_star(i,j)%vertex_star( kside1,2)

                call id_get17( vleft_y,dy,j,ic ,2)
                call id_get17( vright_y,dy,j,id ,2)

                mx = abs(ib-ia)
                my = abs(id-ic)

                zz( 1, :) = Dij_star(i,j)%vertex_star( kside, : )
                zz( mx + my + 2 , : ) = Dij_star(i,j)%vertex_star( kside1, : )
                ! store side boundary
                !                        vl ----> zz( 1, :)
                !                        vr ----> zz( mx + my + 2 , : )
                if(mx .ne. 0)then
                    do k = 1 , mx
                        if(ia<ib) then
                            inter = ia  + k
                        else
                            inter = ia + 1 -k
                        endif

                        zz( 1+k, 1 ) = x0( inter )

                        zz( 1+k, 2 ) = yside( zz( 1, :), zz( mx + my+ 2 , : ),  x0(inter) )

                        if( abs(zz( 1, 1)-zz( mx + my+ 2 , 1 ) )< 1d-13)then  !
                            zz( 1+k, 2 ) = zz(1,2)

                        endif


                        ix = inter - iam

                        if( isx(ix) == 0 )then
                            isx(ix) =1
                            point_inner_x( ix,1,: ) = zz(1+k,:)

                        else
                            point_inner_x( ix,2,: ) = zz(1+k,:)

                        endif

                    enddo
                endif

                ! Step2:

                if(my .ne. 0)then
                    !
                    do k = 1,my
                        if(ic<id)then
                            inter = ic + k
                        else
                            inter = ic  + 1 -k
                        endif
                        zz(1+mx+k,2) = y0( inter )
                        zz( 1+mx+k,1) = xside(  zz( 1, :), zz( mx + my+ 2 , : ),  y0(inter)  )

                        if( abs(zz( 1, 2)-zz( mx + my+ 2 , 2 ) )< 1d-13)then
                            zz( 1+mx+k,1) = zz(1,1)
                            !pause
                        endif

                        iy = inter - icm

                        if(isy(iy)==0 )then
                            point_inner_y(iy,1,:) = zz(1+mx+k,:)
                            isy(iy) = 1

                        else
                            point_inner_y(iy,2,:) = zz(1+mx+k,:)

                        endif

                    enddo
                endif


                ! mx+my /= 0
                ! we do a sorting
                if(mx+my>1)then
                    if(mx .ne. 0)then
                        if( vleft < vright )then
                            ! increase
                            do ii = 1, mx +my -1
                                do jj = 2, 1+ mx+my -ii
                                    if(zz(jj,1 )>zz(jj+1,1) )then
                                        temp(:) = zz(jj+1,:)
                                        zz(jj+1,:) = zz(jj,:)
                                        zz(jj,:) = temp(:)
                                    endif
                                enddo
                            enddo
                        else
                            ! decrease

                            do ii = 1, mx +my -1
                                do jj = 2, 1+ mx+my -ii
                                    if(zz(jj,1 )<zz(jj+1,1) )then
                                        temp(:) = zz(jj+1,:)
                                        zz(jj+1,:) = zz(jj,:)
                                        zz(jj,:) = temp(:)
                                    endif

                                enddo

                            enddo
                        endif
                    else ! mx == 0
                        if(my .ne. 0)then
                            if( vleft_y < vright_y )then
                                ! increase
                                do ii = 1, mx +my -1
                                    do jj = 2, 1+ mx+my -ii
                                        if(zz(jj,2 )>zz(jj+1,2) )then
                                            temp(:) = zz(jj+1,:)
                                            zz(jj+1,:) = zz(jj,:)
                                            zz(jj,:) = temp(:)
                                        endif
                                    enddo
                                enddo
                            else
                                ! decrease
                                do ii = 1, mx +my -1
                                    do jj = 2, 1+ mx+my -ii
                                        if(zz(jj,2 )<zz(jj+1,2) )then
                                            temp(:) = zz(jj+1,:)
                                            zz(jj+1,:) = zz(jj,:)
                                            zz(jj,:) = temp(:)
                                        endif
                                    enddo
                                enddo
                            endif
                        endif
                    endif!mx
                endif!mx+my>1

                ! final outer segments
                do ii = Dij_star(i,j)%nsub_outer +1 , Dij_star(i,j)%nsub_outer +  1+ mx +my
                    jj = ii - Dij_star(i,j)%nsub_outer
                    Dij_star(i,j)%segment_outer( ii )%vl(:) = zz( jj,: )
                    Dij_star(i,j)%segment_outer( ii )%vr(:) = zz( jj+1,: )

                    call id_get17( (zz(jj,1) + zz(jj+1,1) )/2.,dx,i,idx,1 )

                    call id_get17( (zz(jj,2) + zz(jj+1,2) )/2. ,dy,j,idy,2 )


                    Dij_star(i,j)%segment_outer( ii )%id(1) = idx
                    Dij_star(i,j)%segment_outer( ii )%id(2) = idy

                enddo

                Dij_star(i,j)%nsub_outer = Dij_star(i,j)%nsub_outer +  1+ mx +my

            enddo !end ksides

            ! find inner segments

            do ix = 1 , (ibm -iam)
                !
                call id_get17( point_inner_x(ix,1,2),dy,j,i1111,2 )
                call id_get17( point_inner_x(ix,2,2),dy,j,i2222 ,2)

                nby = i1111 - i2222



                temp(:) = point_inner_x(ix,2,:)

                call id_get17( point_inner_x(ix,1,2),dy,j,iy0,2 )


                if( nby >0 )then
                    do ii = 1,nby
                        point_inner_x(ix,ii+1,1) = point_inner_x(ix,1,1)
                        point_inner_x(ix,ii+1,2) = y0( iy0 + 1 - ii)
                    enddo
                elseif(nby < 0)then
                    do ii = 1, abs(nby)
                        point_inner_x(ix,ii+1,1) = point_inner_x(ix,1,1)
                        point_inner_x(ix,ii+1,2) = y0( iy0 + ii)
                    enddo
                else

                endif

                point_inner_x(ix,2+abs(nby), :) = temp(:)

                ! final inner segments
                do ii = 1, 1+ abs(nby)
                    jj = Dij_star(i,j)%nsub_inner + ii
                    Dij_star(i,j)%segment_inner(jj)%vl(:) =  point_inner_x(ix,ii, :)
                    Dij_star(i,j)%segment_inner(jj)%vr(:) =  point_inner_x(ix,ii+1, :)

                    call id_get17( (point_inner_x(ix,ii,2)+point_inner_x(ix,ii+1,2) )/2., dy, j, idy ,2 )                 !
                    Dij_star(i,j)%segment_inner(jj)%id(2) = idy
                    if( point_inner_x(ix,ii,2)>point_inner_x(ix,ii+1,2) )then

                        call id_get_up_down( point_inner_x(ix,ii,1),dx,i,idx )


                        Dij_star(i,j)%segment_inner(jj)%id(1) =  idx

                    else

                        call id_get_down_up( point_inner_x(ix,ii,1),dx,i,idx )
                        Dij_star(i,j)%segment_inner(jj)%id(1) =  idx
                        !
                    endif

                enddo
                !
                Dij_star(i,j)%nsub_inner = Dij_star(i,j)%nsub_inner + abs(nby) + 1
                do ii = 1, 1+ abs(nby)
                    jj = Dij_star(i,j)%nsub_inner + ii
                    Dij_star(i,j)%segment_inner(jj)%vl(:) =  point_inner_x(ix,ii+1, :)
                    Dij_star(i,j)%segment_inner(jj)%vr(:) =  point_inner_x(ix,ii, :)

                    call id_get17( (point_inner_x(ix,ii,2)+point_inner_x(ix,ii+1,2) )/2.,dy,j,idy,2 )

                    Dij_star(i,j)%segment_inner(jj)%id(2) = idy
                    if( point_inner_x(ix,ii+1,2)>point_inner_x(ix,ii,2) )then

                        call id_get_up_down( point_inner_x(ix,ii,1),dx,i,idx )
                        Dij_star(i,j)%segment_inner(jj)%id(1) =  idx
                    else

                        call id_get_down_up( point_inner_x(ix,ii,1),dx,i,idx )
                        Dij_star(i,j)%segment_inner(jj)%id(1) =  idx
                    endif

                enddo
                Dij_star(i,j)%nsub_inner = Dij_star(i,j)%nsub_inner + abs(nby) + 1

            enddo  !do ix = 1 , (ibm -iam)
            !******************************************************************************

            do iy = 1 , (idm -icm)
                call id_get17( point_inner_y(iy,1,1),dx,i,i1111 ,1)
                call id_get17( point_inner_y(iy,2,1),dx,i,i2222 ,1)
                nbx = i1111 - i2222

                temp(:) = point_inner_y(iy,2,:)

                call id_get17( point_inner_y(iy,1,1),dx,i,ix0 ,1 )

                if( nbx >0 )then
                    do ii = 1,nbx
                        point_inner_y(iy,ii+1,2) = point_inner_y(iy,1,2)
                        point_inner_y(iy,ii+1,1) = x0( ix0 + 1 - ii)
                    enddo
                elseif(nbx < 0)then
                    do ii = 1, abs(nbx)
                        point_inner_y(iy,ii+1,2) = point_inner_y(iy,1,2)
                        point_inner_y(iy,ii+1,1) = x0( ix0 + ii)
                    enddo
                else

                endif

                point_inner_y(iy,2+abs(nbx), :) = temp(:)

                ! final inner segments
                do ii = 1, 1+ abs(nbx)
                    jj = Dij_star(i,j)%nsub_inner + ii
                    Dij_star(i,j)%segment_inner(jj)%vl(:) =  point_inner_y(iy,ii, :)
                    Dij_star(i,j)%segment_inner(jj)%vr(:) =  point_inner_y(iy,ii+1, :)

                    call id_get17(  (point_inner_y(iy,ii,1)+point_inner_y(iy,ii+1,1) )/2., dx, i,idx,1 )


                    Dij_star(i,j)%segment_inner(jj)%id(1) = idx
                    if( point_inner_y(iy,ii,1) < point_inner_y(iy,ii+1,1) )then


                        call id_get_left_right( point_inner_y(iy,ii,2),dy,j,idy )
                        Dij_star(i,j)%segment_inner(jj)%id(2) =  idy

                    else

                        call id_get_right_left( point_inner_y(iy,ii,2),dy,j,idy )
                        Dij_star(i,j)%segment_inner(jj)%id(2) =  idy

                    endif

                enddo
                !

                Dij_star(i,j)%nsub_inner = Dij_star(i,j)%nsub_inner + abs(nbx) + 1

                do ii = 1, 1+ abs(nbx)
                    jj = Dij_star(i,j)%nsub_inner + ii
                    Dij_star(i,j)%segment_inner(jj)%vl(:) =  point_inner_y(iy,ii+1, :)
                    Dij_star(i,j)%segment_inner(jj)%vr(:) =  point_inner_y(iy,ii, :)

                    call id_get17( (point_inner_y(iy,ii,1)+point_inner_y(iy,ii+1,1) )/2.,dx,i,idx ,1 )

                    Dij_star(i,j)%segment_inner(jj)%id(1) = idx
                    if( point_inner_y(iy,ii+1,1) < point_inner_y(iy,ii,1) )then

                        call id_get_left_right( point_inner_y(iy,ii,2),dy,j,idy )
                        Dij_star(i,j)%segment_inner(jj)%id(2) =  idy

                    else

                        call id_get_right_left( point_inner_y(iy,ii,2),dy,j,idy )
                        Dij_star(i,j)%segment_inner(jj)%id(2) =  idy

                    endif
                enddo
                !
                Dij_star(i,j)%nsub_inner = Dij_star(i,j)%nsub_inner + abs(nbx) + 1

            enddo  !do iy = 1 , (idm -icm)

            if(n_moment == 1)then
                call green_p0_gauss_eulerian( umod_temp(i,j,1) )
            elseif(n_moment == 3)then
                call green_p1_gauss_eulerian(umod_temp(i,j,1:3) )
            elseif(n_moment ==6)then
                call green_p2_gauss3_eulerian(umod_temp(i,j,1:6) )
            endif
        enddo
    enddo



    do i = 1 ,nx
        do j = 1,ny
            Dij(i,j)%umodal(1:n_moment) =  umod_temp(i,j,1:n_moment)
        enddo
    enddo

    end subroutine search_quadrilateral_final
