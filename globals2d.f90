    module globals2d
    implicit none

    real :: pi,one_third,two_third,over12,over6

    real :: xleft,xright,yleft,yright,xlength,ylength
    real :: dx,dy,overdx,overdy
    integer :: i,j,nx,ny,ighost,kkkk

    real,allocatable :: x(:),y(:),x0(:),y0(:),xgrid(:),ygrid(:)
    integer,allocatable :: itrouble(:,:)
    ! allocate temporary variables
    real,allocatable :: umod_temp(:,:,:)
    ! allocate temporary variables for least square
    real,allocatable :: storage_L(:,:,:,:),storage_U(:,:,:,:)
    real,allocatable :: vert(:,:,:,:)
    real,allocatable :: aaa(:,:,:)
    integer,allocatable :: iswitch(:,:)


    real :: xg(6),wg(6)

    real :: tnum,tprint,cfl,dt
    integer :: nt,irk
    integer :: k

    real :: error1,error2,error3,rr1,er1,rr2,er2,rr3,er3
    real :: er11,er22,er33


    real :: ai(1:10)

    integer :: n_moment

    real :: para_dx3 ,para_limit                     !switch to straight line

    integer :: isearch

    integer :: limiter

    real :: eps
    real :: tvb_m,r0_dgweno,r1_dgweno,r2_dgweno,r3_dgweno,r4_dgweno

    real :: gauss2(2,2),gauss3(3,2)
    !**********************************************************************
    real,allocatable :: center(:,:,:)
    real,allocatable :: x_id_up(:,:,:)
    real,allocatable :: y_id_up(:,:,:)
    end module globals2d