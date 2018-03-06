    subroutine order_DG
    implicit none
    integer :: kk0,kk1
    !  real,allocatable :: uct(:,:)
    real :: den
    integer :: lx,ly
    integer :: inf_x,inf_y
    real :: xrg,yrg



    error1=0.0
    error2=0.0
    error3=0.0

    do kk0=1,nx
        do kk1=1,ny
            do lx=1,6
                do ly=1,6
                    den=abs( exact(x(kk0)+xg(lx)*dx ,y(kk1)+xg(ly)*dy,tprint)  &
                        - polynomial(Dij(kk0,kk1)%umodal(1:n_moment), x(kk0)+xg(lx)*dx,x(kk0),dx,y(kk1)+xg(ly)*dy,y(kk1),dy,n_moment ) )
                    error1=error1+den*wg(lx)*wg(ly)*dx*dy
                    error2=error2+den*den*wg(lx)*wg(ly)*dx*dy
                enddo
            enddo



            do lx=1,6
                do ly=1,6
                    xrg =  x(kk0)+xg(lx)*dx
                    yrg =  y(kk1)+xg(ly)*dy
                    error3=max(error3,abs( exact(xrg ,yrg,tprint) - polynomial( Dij(kk0,kk1)%umodal(1:n_moment),xrg,x(kk0),dx,yrg,y(kk1),dy,n_moment )  ))
                enddo
            enddo

        enddo
    enddo
    error1=error1/( (xright-xleft)*(yright-yleft) )
    error2=sqrt(error2/( (xright-xleft)*(yright-yleft) ))
    if(kkkk.eq.1) write(123,103) nx,ny,error1,error2,error3
    write(*,*) error1,error2,error3
    if(kkkk.gt.1) then
        rr1=log(er11/error1)/log(2.0)
        rr2=log(er22/error2)/log(2.0)
        rr3=log(er33/error3)/log(2.0)
        write(123,102) nx,ny,error1,rr1,error2, rr2,error3, rr3
        write(*,*) nx,ny,rr1,rr2,rr3
    endif
    er11=error1
    er22=error2
    er33=error3


    close(111)


111 format(4(1x,f12.4))
102 format(i6,'*',i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'&',i8,'&',i8,'\\',1x,'\hline')
103 format(i6,'*',i6,1x,3('&',1x,es12.2E2,1x,'&',1x) ,'&',i8,'&',i8,'\\',1x,'\hline')

    end subroutine order_DG

    real function polynomial(a00,x16,xc16,dx16,y16,yc16,dy16,k)
    implicit none
    real,intent(in) :: x16,xc16,dx16
    real,intent(in) :: y16,yc16,dy16
    integer,intent(in) :: k
    real,intent(in) :: a00(k)


    if(k==1)then
        polynomial = a00(1)
    elseif(k==3)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16
    elseif(k==6)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16 &
            + a00(4)*(  ((x16-xc16)/dx16)**2 -1./12. ) + a00(5)*(x16 - xc16 )/dx16*(y16 - yc16 )/dy16 &
            + a00(6)*(  ((y16-yc16)/dy16)**2 -1./12. )
    endif


    end function polynomial