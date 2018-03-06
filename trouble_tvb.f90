    subroutine trouble_tvb
    implicit none

    real :: temp_tvbx,temp_tvby
    real :: temp_x,temp_y

    itrouble(:,:) = 0


    do i = 1 , nx
        do j = 1 , ny
            temp_tvbx = 0.5*Dij(i,j)%umodal(2)
            if( abs(temp_tvbx) < tvb_m * dx**2 )then
                temp_x = temp_tvbx
            else
                temp_x = fminmod3( temp_tvbx , Dij(i+1,j)%umodal(1)-Dij(i,j)%umodal(1), Dij(i,j)%umodal(1)-Dij(i-1,j)%umodal(1)  )
            endif

            temp_tvby = 0.5*Dij(i,j)%umodal(3)
            if( abs(temp_tvby) < tvb_m * dy**2 )then
                temp_y = temp_tvby
            else
                temp_y = fminmod3( temp_tvby , Dij(i,j+1)%umodal(1)-Dij(i,j)%umodal(1), Dij(i,j)%umodal(1)-Dij(i,j-1)%umodal(1)  )
            endif


            if( abs(temp_tvbx - temp_x)>1d-13  .or.  abs(temp_tvby - temp_y)>1d-13 )then
                itrouble(i,j) = 1
            endif

        enddo
    enddo



    end subroutine trouble_tvb

    real function fminmod3(a1,a2,a3)
    implicit none
    real,intent(in) :: a1,a2,a3

    if(sign(1.,a1) .eq. sign(1.,a2) .and. sign(1.,a2) .eq. sign(1.,a3) )then
        fminmod3 = min(abs(a1),abs(a2),abs(a3) ) * sign(1.,a1)
    else
        fminmod3 = 0.
    endif

    end function fminmod3