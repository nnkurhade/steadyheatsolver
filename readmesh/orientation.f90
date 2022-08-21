subroutine orientation(p,q,r,nodes,n,orient)
    implicit none
    integer, intent(in)::p,q,r,n
    double precision, dimension(n,3), intent(in)::nodes
    integer, intent(out):: orient
    double precision :: val
    val=((nodes(q,3)-nodes(p,3))*(nodes(r,2)-nodes(q,2)))-((nodes(q,2)-nodes(p,2))*(nodes(r,3)-nodes(q,3)))
    !print*,p,q,r,val
    if (val==0.0) then
        orient = 0
    else if(val>0.0) then
        orient = 1
    else
        orient = 2
    end if
end subroutine orientation
        