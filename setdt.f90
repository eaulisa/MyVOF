    subroutine setdt
    implicit none

    dt = cfl/( pi*overdx + pi*overdy )

    if(tnum+dt>tprint) dt = tprint- tnum

    tnum = tnum +dt

    end subroutine setdt