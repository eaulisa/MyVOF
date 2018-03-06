    !**********************************************************
    real function trans_to_ex( a16,b16,c16,x16,y16 )
    implicit none
    real,intent(in) :: a16,b16,c16,x16,y16

    trans_to_ex = a16*x16 + b16*y16 + c16

    end function trans_to_ex
    !**********************************************************
    real function trans_to_ey( a16,b16,d16,x16,y16 )
    implicit none
    real,intent(in) :: a16,b16,d16,x16,y16

    trans_to_ey = b16*x16 - a16*y16 + d16

    end function trans_to_ey
    !**********************************************************
    real function reverse_to_x(x16,y16,x26,y26, exi20,eta20)
    implicit none
    real,intent(in) ::  x16,y16,x26,y26, exi20,eta20

    reverse_to_x = (x26-x16)/2.*exi20 + (y26-y16)/2.*eta20 + (x16+x26)/2.


    end function reverse_to_x
    !**********************************************************
    real function reverse_to_y( x16,y16,x26,y26, exi20,eta20 )
    implicit none  
    real,intent(in) :: x16,y16,x26,y26, exi20,eta20

    reverse_to_y = (y26-y16)/2.*exi20 - (x26-x16)/2.*eta20 + (y16+y26)/2.

    end function reverse_to_y

    subroutine trans( x11,y11,x22,y22,a_trans,b_trans,c_trans,d_trans )
    implicit none
    real,intent(in) :: x11,y11,x22,y22
    real,intent(inout) :: a_trans,b_trans,c_trans,d_trans

    if( abs(x22-x11)< abs(y22-y11) )then
        a_trans = 2.*( (x22-x11)/(y22-y11)/(y22-y11) ) / (   ( (x22-x11)/(y22-y11) )**2  + 1.    ) 
        b_trans = 2./( y22-y11 )  /  (   ( (x22-x11)/(y22-y11) )**2  + 1.    ) 
        c_trans = (  (x11-x22)/(y11-y22)*(x11+x22)/(y11-y22)  + (y11+y22)/(y11-y22)  ) /  (   ( (x22-x11)/(y22-y11) )**2  + 1.    ) 
        d_trans = (  (x11+x22)/(y11-y22)  +  (y11+y22)/(y22-y11)*(x22-x11)/(y22-y11)   )  /  (   ( (x22-x11)/(y22-y11) )**2  + 1.    )
    else
        a_trans = 2./(x22-x11) / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
        b_trans = 2.*(y22-y11)/(x22-x11)/(x22-x11)   / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
        c_trans = (  (x11+x22)/(x11-x22) + (y11-y22)/(x11-x22)*(y11+y22)/(x11-x22)  )  / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
        d_trans = (  (x11+x22)/(x11-x22)*(y11-y22)/(x11-x22)  + (y11+y22)/(x22-x11)  )  / (   1.+(   (y22-y11)/(x22-x11)  )**2    )
    endif    
    end subroutine trans