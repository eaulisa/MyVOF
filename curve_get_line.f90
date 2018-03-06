    subroutine get_line( x1,y1,x2,y2, a16,b16,c16 )
    implicit none
    real,intent(in) :: x1,y1,x2,y2
    real,intent(inout) :: a16,b16,c16


    if( abs(x2-x1) > abs(y2-y1) )then
        ! y = kx+ b
        ! kx - y + b = 0
        a16 = (y2-y1)/(x2-x1)
        b16 = - 1.
        c16 = y1-a16*x1
    else
        ! x = ky + b
        ! ky - x + b =0.
        a16 = - 1.
        b16 = (x2-x1)/(y2-y1)
        c16 = x1 - b16*y1
    endif

    end subroutine get_line