    subroutine VOFlimiter
    implicit none
    real :: c1,c2,c3,c4,c5,c6
    real :: minValue,maxValue
    real :: temp, x, y, den, d 
   
    do i = 1 , nx
      do j = 1 ,ny
        c1 = Dij(i,j)%umodal(1)
        c2 = Dij(i,j)%umodal(2)
        c3 = Dij(i,j)%umodal(3)
        c4 = Dij(i,j)%umodal(4)
        c5 = Dij(i,j)%umodal(5)
        c6 = Dij(i,j)%umodal(6)
    
       
    
	! left - bottom Corner
        temp = 1./12. * (12. * c1 - 6.* c2 - 6.* c3 + 2. * c4 + 3. * c5 + 2. * c6)
	minValue = temp
	maxValue = temp
	
	! right - bottom Corner
	temp = 1./12. * (12. * c1 + 6. * c2 - 6. * c3 + 2.* c4 - 3. * c5 + 2. * c6)
	minValue = min (minValue, temp)
	maxValue = max (maxValue, temp)
	
	! right - top Corner
	temp = 1./12. * (12. * c1 + 6. * c2 + 6. * c3 + 2.* c4 + 3. * c5 + 2. * c6)
	minValue = min (minValue, temp)
	maxValue = max (maxValue, temp)
	
	! left - top Corner
	temp = 1./12. * (12. * c1 - 6. * c2 + 6. * c3 + 2.* c4 - 3. * c5 + 2. * c6)
	minValue = min (minValue, temp)
	maxValue = max (maxValue, temp)
	
	
	if ( abs(c4) > 1.0e-10 ) then
	  ! bottom edge
	  x = (-2.* c2 + c5)/(4. * c4)
	  if( x > -0.5 .and. x < 0.5) then
	    temp = -(3. * (-2. * c2 + c5) * (-2. * c2 + c5) + 4. * c4 * (-12.* c1 + 6. * c3 + c4 - 2. * c6)) / (48. * c4)
	    minValue = min (minValue, temp)
	    maxValue = max (maxValue, temp)
	  end if
	  ! top edge
	  x = (-2.* c2 - c5)/(4. * c4)
	  if( x > -0.5 .and. x < 0.5) then
	    temp = (-3. * (2. * c2 + c5) * (2. * c2 + c5) + 4. * c4 * (12.* c1 + 6. * c3 - c4 + 2. * c6)) / (48. * c4)
	    minValue = min (minValue, temp)
	    maxValue = max (maxValue, temp)
	  end if
	end if
	
	if ( abs(c6) > 1.0e-10 ) then
	  ! left edge
	  y = (-2.* c3 + c5)/(4. * c6)
	  if( y > -0.5 .and. y < 0.5) then
	    temp = -(3. * (-2. * c3 + c5) * (-2. * c3 + c5) + 4. * c6 * (-12.* c1 + 6. * c2 + c6 - 2. * c4)) / (48. * c6)
	    minValue = min (minValue, temp)
	    maxValue = max (maxValue, temp)
	  end if
	  ! right edge
	  y = (-2.* c3 - c5)/(4. * c6)
	  if( y > -0.5 .and. y < 0.5) then
	    temp = (-3. * (2. * c3 + c5) * (2. * c3 + c5) + 4. * c6 * (12.* c1 + 6. * c2 - c6 + 2. * c4)) / (48. * c6)
	    minValue = min (minValue, temp)
	    maxValue = max (maxValue, temp)
	  end if
	end if
	
	
	den = c5 * c5 - 4.* c4 * c6;
	if (abs(den) > 1.0e-10 ) then    
	  x = - (c3 * c5 - 2. * c2 * c6) / den
	  y = - (c2 * c5 - 2. * c3 * c4) / den
	  if( x > -0.5 .and. x < 0.5 .and. y > -0.5 .and. y < 0.5) then
	    temp = (12. * c3 * c3 * c4 - 12. * c2 * c3 * c5 +  12. * c2 * c2 * c6 + (12. * c1 - c4 - c6) * (c5 * c5 - 4. * c4 * c6))/( 12. * den )
	    minValue = min (minValue, temp)
	    maxValue = max (maxValue, temp)
	  end if
	end if
	!if(maxValue > 0.5 ) then
	  !print*, c1," ",c2," ",c3," ",c4," ",c5," ",c6
	  !print*, minValue, " ", maxValue, "AAAAAAAAAAAAA"
	!end if
	
	minValue = minValue - 0.5
	maxValue = maxValue - 0.5;
            
        d = 1.    
	if( minValue < 0. .and. maxValue > 0.) then
	  if (minValue > -0.5 .and. maxValue < 0.5) then
	    d = min(1.5 , 0.5 / max( abs(minValue) , maxValue))
	  endif
	else if( minValue >= 0. .and. maxValue < 0.5) then
	  d = min(1.5 , 0.5 / maxValue)
	else if( maxValue <= 0. .and. minValue >-0.5) then
	  d = min(1.5 , 0.5 / abs(minValue) )
	else if( maxValue > 0.5 .and. minValue < -0.5) then
	  d = max(1./1.5 , 0.5 / max( abs(minValue) , maxValue) )
	  print*, minValue, " ",maxValue," ", d
	else if( maxValue > 0.5 .and. minValue > -0.5) then
	  d = max(1./1.5 , 0.5 /  maxValue )
	  print*, minValue, " ",maxValue," ", d
	else if( maxValue < 0.5 .and. minValue < -0.5) then
	  d = max(1./1.5 , 0.5 / abs( minValue )  )
	  print*, minValue, " ",maxValue," ", d
	endif
	
	
	if(d .ne. 1.) then
	  Dij(i,j)%umodal(1:6)  = Dij(i,j)%umodal(1:6) * d
	  Dij(i,j)%umodal(1) = Dij(i,j)%umodal(1) + 0.5 * (1. - d)
	end if
	
	if(minValue > 0.4) then
	  Dij(i,j)%umodal(2:6)  = 0
	  Dij(i,j)%umodal(1) = 1
	else if(maxValue < -0.4) then
	  Dij(i,j)%umodal(1:6)  = 0.
	end if
	
	enddo
      enddo
      

! 
!     do i = 1 , nx
!         do j = 1 ,ny
!             if(itrouble(i,j)==1)then
!                 Dij(i,j)%umodal(2:3)  = Dij(i,j)%temp(1:2)
!                 itrouble(i,j) = 0
!             endif
!         enddo
!     enddo


    end subroutine VOFlimiter