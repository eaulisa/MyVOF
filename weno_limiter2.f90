    subroutine wenolimiter2
    implicit none
    real :: beta0,beta1,beta2,beta3,beta4
    real :: ome0,ome1,ome2,ome3,ome4
    real :: osum
    real :: omega0,omega1,omega2,omega3,omega4


    do i = 1 , nx
        do j = 1 ,ny
            if(itrouble(i,j)==1)then

                beta0 = Dij(i,j)%umodal(2)**2 + Dij(i,j)%umodal(3)**2   &
                    +4./3.* Dij(i,j)%umodal(4)**2  +  5./12.*Dij(i,j)%umodal(5)**2 + 4./3.*Dij(i,j)%umodal(6)**2

                beta1 = 17./12.*Dij(i-1,j)%umodal(5)**2 +4./3.*Dij(i-1,j)%umodal(6)**2 + 16./3.*Dij(i-1,j)%umodal(4)**2 &
                    +Dij(i-1,j)%umodal(2)**2 +4.*Dij(i-1,j)%umodal(2)*Dij(i-1,j)%umodal(4)         &
                    +Dij(i-1,j)%umodal(3)**2 +2.*Dij(i-1,j)%umodal(3)*Dij(i-1,j)%umodal(5)


                beta2 = 17./12.*Dij(i+1,j)%umodal(5)**2 +4./3.*Dij(i+1,j)%umodal(6)**2 + 16./3.*Dij(i+1,j)%umodal(4)**2 &
                    +Dij(i+1,j)%umodal(2)**2 - 4.*Dij(i+1,j)%umodal(2)*Dij(i+1,j)%umodal(4)         &
                    +Dij(i+1,j)%umodal(3)**2 - 2.*Dij(i+1,j)%umodal(3)*Dij(i+1,j)%umodal(5)


                beta3 = 17./12.*Dij(i,j-1)%umodal(5)**2 + 16./3.*Dij(i,j-1)%umodal(6)**2 +4./3.*Dij(i,j-1)%umodal(4)**2 &
                    + Dij(i,j-1)%umodal(2)**2 + 2.*Dij(i,j-1)%umodal(2)*Dij(i,j-1)%umodal(5)  &
                    +Dij(i,j-1)%umodal(3)**2  +4.* Dij(i,j-1)%umodal(3)*Dij(i,j-1)%umodal(6)



                beta4 = 17./12.*Dij(i,j+1)%umodal(5)**2 + 16./3.*Dij(i,j+1)%umodal(6)**2 +4./3.*Dij(i,j+1)%umodal(4)**2 &
                    + Dij(i,j+1)%umodal(2)**2 - 2.*Dij(i,j+1)%umodal(2)*Dij(i,j+1)%umodal(5)  &
                    +Dij(i,j+1)%umodal(3)**2  - 4.* Dij(i,j+1)%umodal(3)*Dij(i,j+1)%umodal(6)



                ome0= r0_dgweno/(eps+ beta0 )**2
                ome1= r1_dgweno/(eps+ beta1 )**2
                ome2= r2_dgweno/(eps+ beta2 )**2
                ome3= r3_dgweno/(eps+ beta3 )**2
                ome4= r4_dgweno/(eps+ beta4 )**2

                osum = ome0+ome1+ome2+ome3+ome4
                omega0 = ome0/osum
                omega1 = ome1/osum
                omega2 = ome2/osum
                omega3 = ome3/osum
                omega4 = ome4/osum


                Dij(i,j)%temp(1) = Omega0 * Dij(i,j)%umodal(2) &
                    +Omega1 * ( Dij(i-1,j)%umodal(2)+2.*Dij(i-1,j)%umodal(4)  )   &
                    +Omega2 * ( Dij(i+1,j)%umodal(2)-2.*Dij(i+1,j)%umodal(4)   )  &
                    +Omega3 * (   Dij(i,j-1)%umodal(2)+Dij(i,j-1)%umodal(5)    )        &
                    +Omega4 * (  Dij(i,j+1)%umodal(2) - Dij(i,j+1)%umodal(5)   )


                Dij(i,j)%temp(2) = Omega0 *   Dij(i,j)%umodal(3)   &
                    +Omega1 * (  Dij(i-1,j)%umodal(3)+Dij(i-1,j)%umodal(5)   ) &
                    +Omega2 * (  Dij(i+1,j)%umodal(3)-Dij(i+1,j)%umodal(5)   ) &
                    +Omega3 * (  Dij(i,j-1)%umodal(3)+2.*Dij(i,j-1)%umodal(6) )  &
                    +Omega4 *  (  Dij(i,j+1)%umodal(3)-2.*Dij(i,j+1)%umodal(6)   )

                Dij(i,j)%temp(3) = Omega0 * Dij(i,j)%umodal(4) +Omega1 * Dij(i-1,j)%umodal(4)+Omega2 * Dij(i+1,j)%umodal(4) &
                    +Omega3 * Dij(i,j-1)%umodal(4)+Omega4 * Dij(i,j+1)%umodal(4)
                Dij(i,j)%temp(4) = Omega0 * Dij(i,j)%umodal(5) +Omega1 * Dij(i-1,j)%umodal(5)+Omega2 * Dij(i+1,j)%umodal(5) &
                    +Omega3 * Dij(i,j-1)%umodal(5)+Omega4 * Dij(i,j+1)%umodal(5)
                Dij(i,j)%temp(5) = Omega0 * Dij(i,j)%umodal(6) +Omega1 * Dij(i-1,j)%umodal(6)+Omega2 * Dij(i+1,j)%umodal(6) &
                    +Omega3 * Dij(i,j-1)%umodal(6)+Omega4 * Dij(i,j+1)%umodal(6)
            endif
        enddo
    enddo

    do i = 1 , nx
        do j = 1 ,ny
            if(itrouble(i,j)==1)then
                Dij(i,j)%umodal(2:6)  = Dij(i,j)%temp(1:5)
                itrouble(i,j) = 0
            endif
        enddo
    enddo

    end subroutine wenolimiter2