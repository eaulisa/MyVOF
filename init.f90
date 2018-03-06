    subroutine init
    implicit none

    real :: utemp
    integer :: lx,ly
    integer :: kk
    real :: cpu1,cpu2
    integer :: icpu


    ! x(i) ----> center
    ! y(j) ----> center
    DO I = -ighost + 1, nx + ighost
        X(I) = XLEFT + (I-0.5) * DX
    ENDDO

    DO J = - ighost + 1, NY+ ighost
        Y(J) = YLEFT + (J-0.5) * DY
    ENDDO



    DO I = -ighost + 1, nx + 1+ ighost
        X0(I) = XLEFT + (I-1.) * DX
    ENDDO

    DO J = - ighost + 1, NY+ 1+ighost
        Y0(J) = YLEFT + (J-1.) * DY
    ENDDO

    DO I = -ighost + 1, nx + 1+ ighost
        Xgrid(I) = XLEFT + (I-1.) * DX
    ENDDO

    DO J = - ighost + 1, NY+ 1+ighost
        ygrid(J) = YLEFT + (J-1.) * DY
    ENDDO



    do i = 1 - 1 , nx + 1
        do j = 1  -1, ny + 1
            do kk = 1 , n_moment
                utemp =0.0
                do lx=1,6
                    do ly=1,6
                        utemp = utemp +exact(x(i)+xg(lx)*dx ,y(j)+xg(ly)*dy,0.) *fphi(kk,xg(lx),xg(ly) )*wg(lx)*wg(ly)
                    enddo
                enddo
                Dij(i,j)%umodal(kk)= utemp * ai(kk)
            enddo
            !
            !anti-clockwise
            Dij(i,j)%vertex(1,1) = xleft + (i-1)*dx
            Dij(i,j)%vertex(1,2) = yleft + (j-1)*dy
            Dij(i,j)%vertex(2,1) = xleft + (i)*dx
            Dij(i,j)%vertex(2,2) = yleft + (j-1)*dy
            Dij(i,j)%vertex(3,1) = xleft + (i)*dx
            Dij(i,j)%vertex(3,2) = yleft + (j)*dy
            Dij(i,j)%vertex(4,1) = xleft + (i-1)*dx
            Dij(i,j)%vertex(4,2) = yleft + (j)*dy

            !*****************************************************
            ! anti-clockwise
            Dij(i,j)%node(1,1) = xleft + (i-1)*dx
            Dij(i,j)%node(1,2) = yleft + (j-1)*dy
            Dij(i,j)%node(2,1) = xleft + (i-0.5)*dx
            Dij(i,j)%node(2,2) = yleft + (j-1)*dy
            Dij(i,j)%node(3,1) = xleft + i*dx
            Dij(i,j)%node(3,2) = yleft + (j-1)*dy

            Dij(i,j)%node(4,1) = xleft + i*dx
            Dij(i,j)%node(4,2) = yleft + (j-0.5)*dy
            Dij(i,j)%node(5,1) = xleft + i*dx
            Dij(i,j)%node(5,2) = yleft + j*dy

            Dij(i,j)%node(6,1) = xleft + (i-0.5)*dx
            Dij(i,j)%node(6,2) = yleft + j*dy
            Dij(i,j)%node(7,1) = xleft + (i-1)*dx
            Dij(i,j)%node(7,2) = yleft + j*dy
            Dij(i,j)%node(8,1) = xleft + (i-1)*dx
            Dij(i,j)%node(8,2) = yleft + (j-0.5)*dy


            Dij(i,j)%node(9,1) = xleft + (i-0.5)*dx
            Dij(i,j)%node(9,2) = yleft + (j-0.5)*dy
        enddo
    enddo


    end subroutine init


    real function exact(x,y,t)
    implicit none
    real, intent(in) :: x,y,t
    real :: rb,rb0

    rb = sqrt( (x-0.15/0.5*pi)**2 + (y)**2  )
    rb0 = 0.3*pi
    if(rb < rb0)then
        exact =  rb0 * cos( pi*rb/(2.*rb0) )**6
    else
        exact = 0.
    endif

    return
    end function exact

    !

    real function fphi(k,x,y)
    implicit none
    integer,intent(in) :: k
    real,intent(in) :: x,y

    if(k .eq. 1)then
        fphi = 1
    elseif(k .eq. 2) then
        fphi = x
    elseif(k .eq. 3)then
        fphi = y
    elseif(k .eq. 4)then
        fphi = x**2 - 1./12.
    elseif(k .eq. 5)then
        fphi = x*y
    elseif(k .eq. 6)then
        fphi = y**2 - 1./12.
    endif


    end function fphi