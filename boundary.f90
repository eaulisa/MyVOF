    subroutine boundary
    !-------------------------------------------subroutine  comment
    !  Purpose   :  periodic boundary condition
    !---------------------------------------------------------------
    implicit none

    do i = 1 , nx
        do j = 1 ,ighost
            Dij(i,1-j)%umodal(:) = Dij(i,ny+1 -j)%umodal(:)
            Dij(i,ny +j)%umodal(:) = Dij(i,0+j)%umodal(:)
        enddo
    enddo

    do j = 1 -ighost , ny +ighost
        do i = 1 , ighost
            Dij(1-i,j)%umodal(:) = Dij(nx+1 -i,j)%umodal(:)
            Dij(nx +i,j)%umodal(:) = Dij(0+i,j)%umodal(:)
        enddo
    enddo


    end subroutine boundary