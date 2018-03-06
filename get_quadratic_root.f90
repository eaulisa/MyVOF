    subroutine get_quadratic_root( a17,b17,c17, xi17 , n_xi )
    implicit none
    real,intent(in) :: a17,b17,c17
    real,intent(out) :: xi17(2)
    integer,intent(out) :: n_xi
    real :: delta,gamma

    if( abs(a17)<3d-14 .and. abs(b17)<1d-15 .and. abs(c17)<3d-14 )then
        !if( a17==0. .and. b17==0. .and. c17==0. )then
        n_xi =0
        !elseif( a17 ==0. .and. b17 == 0. .and. c17 .ne. 0. )then
    elseif( abs(a17)<3d-14  .and. abs(b17)<1d-15 .and. abs(c17)>= 3d-14 )then
        n_xi =0
        !elseif( a17==0. .and. b17 .ne. 0. )then
    elseif( abs(a17)<3d-14 .and. abs(b17) >=1d-15 )then
        n_xi = 1
        xi17(1) = -c17/b17
    else
        delta = b17*b17 - 4.*a17*c17
        if(delta<0.)then
            n_xi =0
        elseif( delta ==0. )then
            n_xi = 1
            xi17(1) = - b17/a17*0.5
        elseif( delta >0. )then
            n_xi = 2
            if(b17>=0.)then
                gamma = 1.
            else
                gamma = -1.
            endif

            xi17(1) = 2.*c17/( -b17 - gamma *sqrt(delta) )
            xi17(2) = ( -b17-gamma*sqrt(delta)  )/a17*0.5
        endif
    endif


    end subroutine get_quadratic_root