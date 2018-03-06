    ! Runge-Kutta
    subroutine RK3( vx16,vx_star16,num16,dt16)
    implicit none
    integer,intent(in) :: num16
    real,intent(in) :: dt16
    real,intent(in) :: vx16(num16,2)
    real :: vx1(num16,2),vx2(num16,2),vx3(num16,2)
    real,intent(out) :: vx_star16(num16,2)
    integer :: ii,ik

    real :: rk_k1(num16,2),rk_k2(num16,2),rk_k3(num16,2),rk_k4(num16,2),rk_k5(num16,2),rk_k6(num16,2)


    ! TVD RK3
    if(irk ==3)then
        do ik = 1 , num16
            vx1(ik,1) = vx16(ik,1) - ax( vx16(ik,1), vx16(ik,2) ,tnum ) * dt16
            vx1(ik,2) = vx16(ik,2) - bx( vx16(ik,1), vx16(ik,2) ,tnum ) * dt16
            vx2(ik,1) = 0.75*vx16(ik,1) + 0.25*(  vx1(ik,1) - ax( vx1(ik,1),vx1(ik,2),tnum-dt16 ) * dt16 )
            vx2(ik,2) = 0.75*vx16(ik,2) + 0.25*(  vx1(ik,2) - bx( vx1(ik,1),vx1(ik,2),tnum-dt16 ) * dt16 )
            vx_star16(ik,1) = one_third*vx16(ik,1) + two_third*( vx2(ik,1) - ax( vx2(ik,1),vx2(ik,2),tnum-dt16*0.5 ) * dt16  )
            vx_star16(ik,2) = one_third*vx16(ik,2) + two_third*( vx2(ik,2) - bx( vx2(ik,1),vx2(ik,2),tnum-dt16*0.5 ) * dt16  )
        enddo
    elseif(irk ==4)then
        do ik = 1 , num16
            vx1(ik,1) = vx16(ik,1) - 0.5* ax( vx16(ik,1),vx16(ik,2),tnum ) * dt16
            vx1(ik,2) = vx16(ik,2) - 0.5* bx( vx16(ik,1),vx16(ik,2),tnum ) * dt16
            vx2(ik,1) = vx16(ik,1) - 0.5* ax( vx1(ik,1),vx1(ik,2),tnum-dt16*0.5 ) * dt16
            vx2(ik,2) = vx16(ik,2) - 0.5* bx( vx1(ik,1),vx1(ik,2),tnum-dt16*0.5 ) * dt16
            vx3(ik,1) = vx16(ik,1) - ax( vx2(ik,1),vx2(ik,2),tnum-dt16*0.5 ) *dt16
            vx3(ik,2) = vx16(ik,2) - bx( vx2(ik,1),vx2(ik,2),tnum-dt16*0.5 ) *dt16
            vx_star16(ik,1) = 1./3.*( -vx16(ik,1) + vx1(ik,1) +2.*vx2(ik,1)+vx3(ik,1) ) - 1./6.*dt16* ax( vx3(ik,1),vx3(ik,2),tnum-dt16  )
            vx_star16(ik,2) = 1./3.*( -vx16(ik,2) + vx1(ik,2) +2.*vx2(ik,2)+vx3(ik,2) ) - 1./6.*dt16* bx( vx3(ik,1),vx3(ik,2),tnum-dt16  )
        enddo
    elseif(irk==5)then
        do ik = 1,num16
            rk_k1(ik,1) = ax( vx16(ik,1),vx16(ik,2),tnum ) 
            rk_k1(ik,2) = bx( vx16(ik,1),vx16(ik,2),tnum )

            rk_k2(ik,1) = ax( vx16(ik,1)-0.25*rk_k1(ik,1)*dt16 , vx16(ik,2)-0.25*rk_k1(ik,2)*dt16   ,tnum-0.25*dt16)
            rk_k2(ik,2) = bx( vx16(ik,1)-0.25*rk_k1(ik,1)*dt16 , vx16(ik,2)-0.25*rk_k1(ik,2)*dt16   ,tnum-0.25*dt16)

            rk_k3(ik,1) = ax( vx16(ik,1)-rk_k1(ik,1)*dt16/8.-rk_k2(ik,1)*dt16/8. , &
            vx16(ik,2)-rk_k1(ik,2)*dt16/8.-rk_k2(ik,2)*dt16/8. ,tnum-0.25*dt16)
            rk_k3(ik,2) = bx( vx16(ik,1)-rk_k1(ik,1)*dt16/8.-rk_k2(ik,1)*dt16/8. , &
            vx16(ik,2)-rk_k1(ik,2)*dt16/8.-rk_k2(ik,2)*dt16/8. ,tnum-0.25*dt16)

            rk_k4(ik,1) = ax(  vx16(ik,1)+rk_k2(ik,1)*dt16/2.-rk_k3(ik,1)*dt16 , &
            vx16(ik,2)+rk_k2(ik,2)*dt16/2.-rk_k3(ik,2)*dt16 ,tnum-0.5*dt16)
            rk_k4(ik,2) = bx( vx16(ik,1)+rk_k2(ik,1)*dt16/2.-rk_k3(ik,1)*dt16 , &
            vx16(ik,2)+rk_k2(ik,2)*dt16/2.-rk_k3(ik,2)*dt16 ,tnum-0.5*dt16)

            rk_k5(ik,1) = ax(  vx16(ik,1)-rk_k1(ik,1)*dt16*3./16.-rk_k4(ik,1)*dt16*9./16. , &
            vx16(ik,2)-rk_k1(ik,2)*dt16*3./16.-rk_k4(ik,2)*dt16*9./16. ,tnum-0.75*dt16)
            rk_k5(ik,2) = bx( vx16(ik,1)-rk_k1(ik,1)*dt16*3./16.-rk_k4(ik,1)*dt16*9./16. , &
            vx16(ik,2)-rk_k1(ik,2)*dt16*3./16.-rk_k4(ik,2)*dt16*9./16. ,tnum-0.75*dt16)

            rk_k6(ik,1) = ax( vx16(ik,1)+rk_k1(ik,1)*dt16*3./7.-rk_k2(ik,1)*dt16*2./7.-rk_k3(ik,1)*dt16*12./7.+rk_k4(ik,1)*dt16*12./7. &
            - rk_k5(ik,1)*dt16*8./7. , &
            vx16(ik,2)+rk_k1(ik,2)*dt16*3./7.-rk_k2(ik,2)*dt16*2./7.-rk_k3(ik,2)*dt16*12./7.+rk_k4(ik,2)*dt16*12./7. &
            - rk_k5(ik,2)*dt16*8./7. ,tnum-dt16)
            rk_k6(ik,2) = bx(  vx16(ik,1)+rk_k1(ik,1)*dt16*3./7.-rk_k2(ik,1)*dt16*2./7.-rk_k3(ik,1)*dt16*12./7.+rk_k4(ik,1)*dt16*12./7. &
            - rk_k5(ik,1)*dt16*8./7. , &
            vx16(ik,2)+rk_k1(ik,2)*dt16*3./7.-rk_k2(ik,2)*dt16*2./7.-rk_k3(ik,2)*dt16*12./7.+rk_k4(ik,2)*dt16*12./7. &
            - rk_k5(ik,2)*dt16*8./7. ,tnum-dt16)

            vx_star16(ik,1) = vx16(ik,1) - (7.*rk_k1(ik,1)+32.*rk_k3(ik,1)+12.*rk_k4(ik,1)+32.*rk_k5(ik,1)+7.*rk_k6(ik,1) )*dt16/90.
            vx_star16(ik,2) = vx16(ik,2) - (7.*rk_k1(ik,2)+32.*rk_k3(ik,2)+12.*rk_k4(ik,2)+32.*rk_k5(ik,2)+7.*rk_k6(ik,2) )*dt16/90.
        enddo

    endif



    end subroutine RK3


    ! problem 
    real function  ax( x16 ,y16 ,t16)
    implicit none

    real,intent(in) :: x16,y16,t16

    ax = -cos(x16/2)**2 * sin(y16) * cos(pi*t16/tprint )*pi


    end function ax

    real function  bx( x16 ,y16 ,t16)
    implicit none

    real,intent(in) :: x16,y16,t16

    bx = sin(x16)*cos(y16/2.)**2 * cos(pi*t16/tprint )*pi

    end function bx
