    subroutine wenolimiter1
    implicit none
    real :: beta0,beta1,beta2,beta3,beta4
    real :: ome0,ome1,ome2,ome3,ome4
    real :: osum
    real :: omega0,omega1,omega2,omega3,omega4



    do i = 1 , nx
        do j = 1 ,ny
            if(itrouble(i,j)==1)then
                beta0 = Dij(i,j)%umodal(2)**2 + Dij(i,j)%umodal(3)**2
                beta1 = Dij(i-1,j)%umodal(2)**2 + Dij(i-1,j)%umodal(3)**2
                beta2 = Dij(i+1,j)%umodal(2)**2 + Dij(i+1,j)%umodal(3)**2
                beta3 = Dij(i,j-1)%umodal(2)**2 + Dij(i,j-1)%umodal(3)**2
                beta4 = Dij(i,j+1)%umodal(2)**2 + Dij(i,j+1)%umodal(3)**2

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

                Dij(i,j)%temp(1) = Omega0 * Dij(i,j)%umodal(2) +Omega1 * Dij(i-1,j)%umodal(2)+Omega2 * Dij(i+1,j)%umodal(2) &
                    +Omega3 * Dij(i,j-1)%umodal(2)+Omega4 * Dij(i,j+1)%umodal(2)
                Dij(i,j)%temp(2) = Omega0 * Dij(i,j)%umodal(3) +Omega1 * Dij(i-1,j)%umodal(3)+Omega2 * Dij(i+1,j)%umodal(3) &
                    +Omega3 * Dij(i,j-1)%umodal(3)+Omega4 * Dij(i,j+1)%umodal(3)
            endif

        enddo
    enddo

    do i = 1 , nx
        do j = 1 ,ny
            if(itrouble(i,j)==1)then
                Dij(i,j)%umodal(2:3)  = Dij(i,j)%temp(1:2)
                itrouble(i,j) = 0
            endif
        enddo
    enddo


    end subroutine wenolimiter1