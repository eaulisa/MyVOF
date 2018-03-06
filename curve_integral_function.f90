    real function curve_dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16


    curve_dx =  (e16**2*a1+e16*b1+c1)*(2.*e16 *a2+b2)

    end function curve_dx
    !
    real function curve_xdx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

    curve_xdx = (1./2.)*(e16**2*a1+e16 *b1+c1)**2*(2.*e16 *a2+b2)

    end function curve_xdx
    !
    real function curve_ydx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16


    curve_ydx = (1./2.)*(e16**2*a2+e16 *b2+c2)**2*(2.*e16*a1+b1)

    curve_ydx = -curve_ydx

    end function curve_ydx
    !
    real function curve_x2dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

  
    curve_x2dx = (1./3.)*(e16**2*a1+e16*b1+c1)**3*(2*e16*a2+b2)

    end function curve_x2dx
    !
    real function curve_xydx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

   
    curve_xydx = (1./2.)*(e16**2*a1+e16*b1+c1)**2*(e16**2*a2+e16*b2+c2)*(2.*e16*a2+b2)
    
    end function curve_xydx
    !
    real function curve_y2dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16


    curve_y2dx = (1./3.)*(e16**2*a2+e16*b2+c2)**3*(2.*e16*a1+b1)

    curve_y2dx = - curve_y2dx



    end function curve_y2dx
    !
    real function curve_x3dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16



    curve_x3dx = (1./4.)*(e16**2*a1+e16*b1+c1)**4*(2.*e16*a2+b2)

    end function curve_x3dx
    !
    real function curve_x2ydx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

    

    curve_x2ydx = (1./3.)*(e16**2*a1+e16*b1+c1)**3*(e16**2*a2+e16*b2+c2)*(2.*e16*a2+b2)


    end function curve_x2ydx
    !
    real function curve_xy2dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16


    curve_xy2dx = (1./3.)*(e16**2*a2+e16*b2+c2)**3*(e16**2*a1+e16*b1+c1)*(2.*e16*a1+b1)        

    curve_xy2dx =  -curve_xy2dx


    end function curve_xy2dx
    !
    real function curve_y3dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

  

    curve_y3dx =    (1./4.)*(e16**2*a2+e16*b2+c2)**4*(2.*e16*a1+b1)

    curve_y3dx = - curve_y3dx 

    end function curve_y3dx
    !
    real function curve_x4dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16


    curve_x4dx =   (1./5.)*(e16**2*a1+e16*b1+c1)**5*(2.*e16*a2+b2)
    
    end function curve_x4dx
    !
    real function curve_x3ydx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

  

    curve_x3ydx =     (1./4.)*(e16**2*a1+e16*b1+c1)**4*(e16**2*a2+e16*b2+c2)*(2.*e16*a2+b2)


    end function curve_x3ydx
    !
    real function curve_x2y2dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16
    
    
    curve_x2y2dx = (1./3.)*(e16**2*a1+e16*b1+c1)**3*(e16**2*a2+e16*b2+c2)**2*(2.*e16*a2+b2)


    end function curve_x2y2dx
    !
    real function curve_xy3dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

  
    curve_xy3dx = (1./4.)*(e16**2*a2+e16*b2+c2)**4*(e16**2*a1+e16*b1+c1)*(2.*e16*a1+b1)
    curve_xy3dx = -curve_xy3dx 

    end function curve_xy3dx
    !
    real function curve_y4dx( a1,b1,c1,a2,b2,c2,e16 )
    implicit none
    real,intent(in) :: a1,b1,c1,a2,b2,c2,e16

    curve_y4dx =    (1./5.)*(e16 **2*a2+e16*b2+c2)**5*(2.*e16*a1+b1)
    curve_y4dx = - curve_y4dx



    end function curve_y4dx
