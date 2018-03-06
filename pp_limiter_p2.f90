    subroutine pp_limiter_p2
    implicit none
    real :: pmin
    real :: theta
    integer :: im

    real :: uhat1,uhat2,uhat3,uhat4,uhat5,uhat6
    real :: xmix,ymix
    real :: x111,y111,x222,y222,x333,y333,x444,y444

    do i = 1 , nx
        do j = 1 ,ny
            pmin= 1000.
            !
            do im = 1,4 
                pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), Dij(i,j)%vertex(im,1),x(i),dx,Dij(i,j)%vertex(im,2),y(j),dy,n_moment ) )
            enddo


            !x = -(2*dx*uhat2*uhat6-dx*uhat3*uhat5-4*uhat4*uhat6*x(i)+uhat5**2*x(i))/(4*uhat4*uhat6-uhat5**2),
            !y = (dy*uhat2*uhat5-2*dy*uhat3*uhat4+4*uhat4*uhat6*y(j)-uhat5**2*y(j))/(4*uhat4*uhat6-uhat5**2)
            uhat1 = Dij(i,j)%umodal(1)
            uhat2 = Dij(i,j)%umodal(2)
            uhat3 = Dij(i,j)%umodal(3)
            uhat4 = Dij(i,j)%umodal(4)
            uhat5 = Dij(i,j)%umodal(5)
            uhat6 = Dij(i,j)%umodal(6)

            if( 4.*uhat4*uhat6-uhat5**2 .ne. 0.  )then
                xmix = -(2.*dx*uhat2*uhat6-dx*uhat3*uhat5-4.*uhat4*uhat6*x(i)+uhat5**2*x(i))/(4.*uhat4*uhat6-uhat5**2)
                ymix = (dy*uhat2*uhat5-2.*dy*uhat3*uhat4+4.*uhat4*uhat6*y(j)-uhat5**2*y(j))/(4.*uhat4*uhat6-uhat5**2)
                if( (xmix-x0(i)+1d-14)*(xmix-x0(i+1)+1d-14)<=0. .and.  (ymix-y0(j)+1d-14)*(ymix-y0(j+1)+1d-14)<=0.   )then
                    !
                    pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), xmix,x(i),dx,ymix,y(j),dy,n_moment ) )
                endif

            endif

            ! x = x_{-1/2}
            ! y = (1/4)*(-2*dy*uhat3+dy*uhat5+4*uhat6*y(j))/uhat6

            ! x = x_{1/2}
            ! y = -(1/4)*(2*dy*uhat3+dy*uhat5-4*uhat6*y(j))/uhat6

            if( uhat6 .ne. 0. )then
                x111 = x0(i) - 1d-14
                y111 = (1./4.)*(-2.*dy*uhat3+dy*uhat5+4.*uhat6*y(j))/uhat6
                if( (y111-y0(j)+1d-14)*(y111-y0(j+1)+1d-14)<=0. )then
                    pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), x111,x(i),dx,y111,y(j),dy,n_moment ) )
                endif

                x222 =  x0(i+1) - 1d-14
                y222 = -(1./4.)*(2.*dy*uhat3+dy*uhat5-4.*uhat6*y(j))/uhat6
                if( (y222-y0(j)+1d-14)*(y222-y0(j+1)+1d-14)<=0. )then
                    pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), x222,x(i),dx,y222,y(j),dy,n_moment ) )
                endif

            endif


            ! x = (1/4)*(-2*dx*uhat2+dx*uhat5+4*uhat4*x(i))/uhat4
            ! y = y_{-1/2}

            ! x = -(1/4)*(2*dx*uhat2+dx*uhat5-4*uhat4*x(i))/uhat4
            ! y = y_{1/2}

            if( uhat4 .ne. 0. )then
                x333 = (1./4.)*(-2*dx*uhat2+dx*uhat5+4*uhat4*x(i))/uhat4
                y333 = y0(j) - 1d-14
                if( (x333-x0(i)+1d-14)*(x333-x0(i+1)+1d-14)<=0. )then
                    pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), x333,x(i),dx,y333,y(j),dy,n_moment ) )       
                endif

                x444 = -(1./4.)*(2*dx*uhat2+dx*uhat5-4*uhat4*x(i))/uhat4
                y444 = y0(j+1) - 1d-14
                if( (x444-x0(i)+1d-14)*(x444-x0(i+1)+1d-14)<=0.  )then
                    pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), x444,x(i),dx,y444,y(j),dy,n_moment ) ) 
                endif 

            endif




            if( pmin ==Dij(i,j)%umodal(1) )then
                theta = 1.
            else
                theta = min(1., abs(  ( Dij(i,j)%umodal(1) ) /(pmin-Dij(i,j)%umodal(1)  )  )       )
            endif



            Dij(i,j)%umodal(2:6)  = theta* Dij(i,j)%umodal(2:6)


        enddo
    enddo



    end subroutine pp_limiter_p2