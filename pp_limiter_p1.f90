    subroutine pp_limiter_p1
    implicit none
    real :: pmin
    real :: theta
    integer :: im

    do i = 1 , nx
        do j = 1 ,ny
            pmin= 1000.
            !
            do im = 1,4 
                pmin = min(pmin,polynomial(Dij(i,j)%umodal(1:n_moment), Dij(i,j)%vertex(im,1),x(i),dx,Dij(i,j)%vertex(im,2),y(j),dy,n_moment ) )
            enddo

            if( pmin ==Dij(i,j)%umodal(1) )then
                theta = 1.
            else
                theta = min(1., abs(  Dij(i,j)%umodal(1) /(pmin-Dij(i,j)%umodal(1)  )  )       )
            endif



            Dij(i,j)%umodal(2:3)  = theta* Dij(i,j)%umodal(2:3)


        enddo
    enddo



    end subroutine pp_limiter_p1