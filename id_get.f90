

    subroutine id_get_left_right( y16,dy16,id16,id17 )
    implicit none
    real,intent(in) :: y16,dy16
    integer,intent(in) :: id16
    integer,intent(out) :: id17
    integer :: i00,id1703

    id1703 = id_get( (y16-yleft)/dy16 )
    do i00 = id1703 +1-2, id1703 + 2+1
        if(  y16 < y(i00) .and. y16> ygrid(i00) - 1d-14 )then
            id17 = i00
        endif
    enddo

    end subroutine id_get_left_right

    subroutine id_get_right_left( y16,dy16,id16,id17 )
    implicit none
    real,intent(in) :: y16,dy16
    integer,intent(in) :: id16
    integer,intent(out) :: id17
    integer :: i00,id1703

    id1703 = id_get( (y16-yleft)/dy16 )
    do i00 = id1703 +1-2, id1703 + 2+1
        if(  y16 > y(i00) .and. y16 < ygrid(i00+1) + 1d-14 )then
            id17 = i00
        endif
    enddo

    end subroutine id_get_right_left
    !*****************
    subroutine id_get_down_up( x16,dx16,id16,id17 )
    implicit none
    real,intent(in) :: x16,dx16
    integer,intent(in) :: id16
    integer,intent(out) :: id17
    integer :: i00,id1703

    id1703 = id_get( (x16-xleft)/dx16 )
    do i00 = id1703 +1-2, id1703 + 2+1
        if(  x16 > x(i00) .and. x16 < xgrid(i00+1) + 1d-14 )then
            id17 = i00
        endif
    enddo

    end subroutine id_get_down_up

    subroutine id_get_up_down( x16,dx16,id16,id17 )
    implicit none
    real,intent(in) :: x16,dx16
    integer,intent(in) :: id16
    integer,intent(out) :: id17
    integer :: i00,id1703

    id1703 = id_get( (x16-xleft)/dx16 )
    do i00 = id1703 +1-2, id1703 + 2+1
        if(  x16 < x(i00) .and. x16 > xgrid(i00) -1d-14 )then
            id17 = i00
        endif
    enddo

    end subroutine id_get_up_down


    subroutine  id_get17(xy,dxy,id16, id17,ixy)
    implicit none
    real,intent(in) :: xy,dxy
    integer,intent(in) :: id16,ixy
    integer,intent(out) :: id17
    integer :: i00

    integer :: id1703
    if(ixy == 1)then
        id1703 = id_get( (xy-xleft)/dxy )
        do i00 = id1703 +1-2, id1703 + 2+2
            if( abs( xy - xgrid(i00) ) .le. dxy + 3d-16 .and. xy-xgrid(i00) .le. 0. )then
                id17 = i00-1
            endif
        enddo
    else
        id1703 = id_get( (xy-yleft)/dxy )
        do i00 = id1703 +1-2, id1703 +2 +2
            if( abs( xy - ygrid(i00) ) .le. dxy + 3d-16 .and. xy-ygrid(i00) .le. 0. )then
                id17 = i00-1
                !print *,123
            endif
        enddo

    endif



    end subroutine id_get17

    integer function id_get(x16)
    !*******************************************************************************
    !
    !   Purpose : getting the location index of point in a Eulerian cell
    !
    !   Notice  : Background cell and Eulerian cell is different in the implementation
    !             Background grid : x(i) + 1d-14
    !             Eulerian cell : x(i)
    !
    !*******************************************************************************
    implicit none
    real,intent(in) :: x16

    if(x16<0.)then
        id_get = int(x16)
    else
        id_get = int(x16) + 1
    endif

    end function id_get