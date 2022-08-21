subroutine left_index(nodes,n,minn)
    implicit none
    integer, intent(in)                                 :: n
    double precision,dimension(n,3), intent(in)         :: nodes
    integer                                             :: i
    integer,intent(out)                                 :: minn
    minn = 1
    do i=1,n
        if (nodes(i,2)<nodes(minn,2)) then
            minn = i
        else if (nodes(i,2)==nodes(minn,2)) then
            if (nodes(i,3)>nodes(minn,3)) then
                minn = i
            end if
        end if
    end do
end subroutine left_index
    