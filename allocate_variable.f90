    subroutine allocate_variable
    implicit none

    !*****************************************************************************
    ! global variables
    ! x(i) ----> center
    ! y(j) ----> center
    allocate( x(-ighost + 1 : nx + ighost) )
    allocate( y(-ighost + 1 : ny + ighost) )

    allocate( x0(-ighost + 1: nx + 1+ighost) )
    allocate( y0(-ighost + 1: ny + 1+ighost) )

    allocate( xgrid(-ighost + 1: nx + 1+ighost) )
    allocate( ygrid(-ighost + 1: ny + 1+ighost) )
    if( limiter == 1 )then
        allocate( itrouble( 1:nx , 1:ny) )
        itrouble(1:nx,1:ny) = 0
    endif

    ! Eulerian cells
    allocate( Dij( 1-ighost:nx+ighost , 1-ighost:ny+ighost) )
    ! Upstream cells
    allocate( Dij_star( 1-ighost:nx+ighost , 1-ighost:ny+ighost) )
    !*****************************************************************************
    ! temporary variables
    allocate( umod_temp( 1-ighost:nx+ighost , 1-ighost:ny+ighost , 1:n_moment )  )
    !*****************************************************************************
    ! local variable in Green for least square
    if(n_moment==3)then
        allocate( aaa(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,3) )
        allocate( iswitch(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost) )
        iswitch(:,:) = 0
        allocate( storage_L(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,1:3,1:3) )
        allocate( storage_U(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,1:3,1:3) )
        allocate( vert(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,1:4,1:2) )
    elseif(n_moment ==6)then
        allocate( aaa(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,6) )
        allocate( iswitch(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost) )
        iswitch(:,:) = 0
        allocate( storage_L(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,1:6,1:6) )
        allocate( storage_U(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,1:6,1:6) )
        allocate( vert(-ighost + 1: nx + ighost,-ighost + 1: ny + ighost,1:9,1:5) )
    endif
    !*****************************************************************************
    
    if( isearch== 2 )then
        allocate( center(-1 + 1 : nx +1, -1+1:ny+1,2 ) )
        allocate( x_id_up(-1 + 1 : nx +1, -1+1:ny+1,2 ) )
        allocate( y_id_up(-1 + 1 : nx +1, -1+1:ny+1,2 ) )
    endif
    end subroutine allocate_variable

    subroutine deallocate_variable
    implicit none

    !*****************************************************************************
    ! global variables
    ! x(i) ----> center
    ! y(j) ----> center
    deallocate( x  )
    deallocate( y )

    deallocate( x0 )
    deallocate( y0 )

    deallocate( xgrid )
    deallocate( ygrid )
    if( limiter == 1 )then
        deallocate( itrouble  )
    endif

    ! Eulerian cells
    deallocate( Dij )
    ! Upstream cells
    deallocate( Dij_star )
    !*****************************************************************************
    ! temporary variables
    deallocate( umod_temp  )
    !*****************************************************************************
    ! local variable in Green for least square
    if(n_moment>1)then
        deallocate(aaa)
        deallocate(iswitch)
        deallocate(storage_L)
        deallocate(storage_U)
        deallocate(vert)
    endif
    !*****************************************************************************
    
    if( isearch== 2 )then
        deallocate( center  )
        deallocate( x_id_up  )
        deallocate( y_id_up  )
    endif
    end subroutine deallocate_variable