    subroutine get_matrix( ver,ps,at,bt )
    !-------------------------------------------subroutine  comment
    !  Purpose   :  get matrix for least square
    !               Ax=b
    !--------------------------------------------------------------
    !  Input  parameters  :
    !       1.    A        -----------> ver
    !       2.    b        -----------> b
    !  Output parameters  :
    !       1.   A^T A      ----------> at
    !       2.   A^T b      ----------> bt
    !---------------------------------------------------------------
    implicit none
    real,intent(in) :: ver(1:4,1:2),ps(1:4)
    real,intent(out) :: at(1:3,1:3),bt(1:3)

    at(1,1) = 4.; at(1,2) = ver(1,1)+ver(2,1)+ver(3,1)+ver(4,1); at(1,3)=ver(1,2)+ver(2,2)+ver(3,2)+ver(4,2); 
    at(2,1) = at(1,2);
    at(2,2) = ver(1,1)**2+ver(2,1)**2+ver(3,1)**2+ver(4,1)**2;
    at(2,3) = ver(1,1)*ver(1,2)+ver(2,1)*ver(2,2)+ver(3,1)*ver(3,2)+ver(4,1)*ver(4,2);
    at(3,1) = at(1,3);
    at(3,2) = at(2,3);
    at(3,3) = ver(1,2)**2+ver(2,2)**2+ver(3,2)**2+ver(4,2)**2;

    !---------------------------------------------------------------

    bt(1) = ps(1) + ps(2)+ps(3)+ps(4);
    bt(2) = ps(1)*ver(1,1) + ps(2)*ver(2,1)+ps(3)*ver(3,1)+ps(4)*ver(4,1);
    bt(3) = ps(1)*ver(1,2) + ps(2)*ver(2,2)+ps(3)*ver(3,2)+ps(4)*ver(4,2);

    end subroutine  get_matrix
