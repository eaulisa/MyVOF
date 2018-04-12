    subroutine setdt
    implicit none

    real denominator
    denominator = abs( cos(pi*tnum/tprint ) )
    if (denominator < 1.e-5) denominator = 1.e-05
    !denominator = 1.;
    
    dt = cfl/( pi*overdx + pi*overdy ) / denominator

    if(tnum+dt>tprint) dt = tprint- tnum

    tnum = tnum +dt

    end subroutine setdt