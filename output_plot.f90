    subroutine output_tecplot
    implicit none
    integer :: lx,ly,ic,jc
    real :: xrg,yrg

    open(121,file='tecplot_2d.plt')
    write(121,*)'zone    ','i=',nx*3,',    j=',ny*3

    DO  J=1,NY
        do jc = 1,3
            DO  I=1 ,NX
                do ic = 1,3
                    xrg =  x(i)+ (ic-2.) *dx*0.25
                    yrg =  y(j)+ (jc-2.) *dy*0.25
                    WRITE(121,*) xrg,yrg,  polynomial( Dij(i,j)%umodal(1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )
                enddo
            enddo
        enddo
    enddo

    CLOSE(121)


    end subroutine output_tecplot
    !**************************************************
    subroutine output_paraview (kkkk, iteration)
    integer kkkk
    integer iteration
    integer :: lx,ly,ic,jc
    real :: xrg,yrg
    
    character(len=1024) :: filename
    character(len=100) :: format_string
    if(iteration < 10) then
      format_string = "(a21, i1, a1, i1, a4)"
    else if (iteration < 100) then
      format_string = "(a21, i1, a1, i2, a4)"
    else if (iteration < 1000) then
      format_string = "(a21, i1, a1, i3, a4)"
    else if (iteration < 1000) then
      format_string = "(a21, i1, a1, i4, a4)"
    else if (iteration < 10000) then
      format_string = "(a21, i1, a1, i5, a4)"
    else if (iteration < 100000) then
      format_string = "(a21, i1, a1, i6, a4)"
    endif
    write (filename,format_string) "./output/paraview_2d.", kkkk,".", iteration, ".vtk"
    open(121,file=trim(filename))
    !open(121,file='./output/paraview_2d.vtk')
    !OPEN(unit=121, status="old", file=fileplace//"paraview_2d.vtk")
    write(121,'(a)')'# vtk DataFile Version 2.0'
    write(121,'(a)')'Cube example'
    write(121,*)
    write(121,'(a)')'ASCII'
    write(121,'(a)')'DATASET RECTILINEAR_GRID'
    write(121,'(a)')'DIMENSIONS'
    write(121,'(i6,i6,i6)') nx*3,ny*3,1
    write(121,*) 'X_COORDINATES',nx*3,'float'
    do i = 1 , nx
        do ic = 1 ,3
            xrg = x(i) + (ic-2.)*dx*0.25
            WRITE(121,'(3e18.6)') xrg
        enddo
    enddo
    write(121,*)'Y_COORDINATES',ny*3,'float'
    do j= 1 ,ny
        do jc = 1 ,3 
            yrg = y(j) + (jc-2.)*dy*0.25
            write(121,*) yrg
        enddo
    enddo
     write(121,*)'Z_COORDINATES',1,'float'
     write(121,*) 0.
     
     write(121,*) 
     write(121,*)
   write(121,*)  'CELL_DATA', (3*nx-1)*(3*ny-1)
write(121,*) 'POINT_DATA',3*nx*3*ny
write(121,*) 'SCALARS nodal float'
write(121,*) 'LOOKUP_TABLE default'
    
    DO  J=1,NY
        do jc = 1,3
            DO  I=1 ,NX
                do ic = 1,3
                    xrg =  x(i)+ (ic-2.) *dx*0.25
                    yrg =  y(j)+ (jc-2.) *dy*0.25
                    WRITE(121,*)  polynomial( Dij(i,j)%umodal(1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )
                enddo
            enddo
        enddo
    enddo

    CLOSE(121)


    end subroutine output_paraview
    !**************************************************
    subroutine output_matlab_plot
    implicit none
    integer :: lx,ly,ic,jc
    real :: xrg,yrg

    open(121,file='matlab_plot2d.plt')

    DO  J=1,NY
        do jc = 1,3
            DO  I=1 ,NX
                do ic = 1,3
                    xrg =  x(i)+ (ic-2.) *dx*0.25
                    yrg =  y(j)+ (jc-2.) *dy*0.25
                    WRITE(121,*) xrg,yrg,  polynomial( Dij(i,j)%umodal(1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )
                enddo
            enddo
        enddo
    enddo

    CLOSE(121)


    end subroutine output_matlab_plot