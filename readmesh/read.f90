program reading
    implicit none
    !declarations
    integer :: i,j,k,n_entity,n_count,min_count,max_count,node_d,sur_id,n_para,n_psudo,count,count_in,count_grad
    integer:: ele_d,ele_sur_id,ele_type_id,ele_count
    double precision,dimension (:,:),allocatable :: nodes,centroids,face_centroids
    integer,dimension(:,:),allocatable:: elements,faces_fake,faces
    integer, dimension(:),allocatable :: psudoid,mask_face
    double precision,dimension(:,:),allocatable :: phi,phi_f,gma,delta_eta,delta_zhi,af_dot_af,af,e_zhi,af_dot_e_zhi,delv,af_vector,grad_phi,del_r,grad_phi_dash,grad_phi_f,e_zhi_cell
    double precision, dimension(:,:),allocatable:: anb_f,sgrad,phi_dash,sc,sp,ap,phi_dash_gs,b,r_zhi_cell,dummy,m,mt,d,dummy2,dummy3,d_dash,q_flux,line_fake,line,phi_star
    double precision :: x1,y1,x2,y2,x0,y0,x,y,error,error_in,x_q,y_q,y_avg,tolerance_avg,alpha
    integer :: minn,p,q,orient,n,fx
    integer, dimension(:,:), allocatable::hull_dummy,hull
    logical :: neuman,outer_dirichlet,inner_dirichlet
    !###########################################################################################
    !node file read operation
    !counter for node counting
    print*, 'please keep patience, configuring geometry will take some time depending upon the mesh complexity'
    k=1
    open(2,file = 'nodes_c.txt',status = 'old')
    read(2,*) n_entity,n_count,min_count,max_count
    allocate(nodes(n_count,4))
    do i=1,n_entity
        read(2,*) node_d,sur_id,n_para,n_psudo
        allocate(psudoid(n_psudo))
        do j=1,n_psudo
            read(2,*) psudoid(j)
            nodes(k,1) = psudoid(j)
            k = k+1
        end do
        k = k-n_psudo
        do j=1,n_psudo
            read(2,*) nodes(k,2),nodes(k,3),nodes(k,4)
            k = k+1
        end do
        deallocate(psudoid)
    end do
    close(2)
    print*, n_count,'nodes detected'
    !node file read operation end
    !##################################################################################################
    !optional code to view nodes
    !do i = 1,n_count
    !    print*,nodes(i,:)
    !end do
    !###################################################################################################
    !code to determine the outer boundary forming nodes
    !__________________________________________________________________________________________________
    !find the top left point (that lies certainly on the hull
    call left_index(nodes,size(nodes(:,1)),minn)
    !__________________________________________________________________________________________________
    !optional code to view the top left point
    !print*,minn
    !__________________________________________________________________________________________________
    !allocate a large array of outermost points, p.s. why fortran whyyyy
    allocate(hull_dummy(size(nodes(:,1)),1))
    hull_dummy = 0
    !___________________________________________________________________________________________________
    !modified jarvis march algorithm
    p = minn
    q = 1
    i=1
    n = size(nodes(:,1))
    do while (.true.)
        hull_dummy(i,1)=p
        q=mod(p+1,n)
        !print*,p,q
        do j =1,n
            !print*,j
            call orientation(p,j,q,nodes,n,orient)
            if(orient==2) then
                q = j
            end if
            !modicication to the jarvis march algorithm found in literature, this fixes a bug it has with handling multiple points on straight edge
            if((orient==0.0).and.(j/=q).and.(p/=q).and.(p/=j)) then
                !print*,'yes'
                if(norm2(nodes(p,2:3)-nodes(j,2:3))+norm2(nodes(j,2:3)-nodes(q,2:3))==norm2(nodes(p,2:3)-nodes(q,2:3)))then
                    if(norm2(nodes(p,2:3)-nodes(j,2:3))<norm2(nodes(p,2:3)-nodes(q,2:3))) then
                        !print*,'yes'
                        q=j
                    end if
                end if
            end if
        end do
        p = q
        i=i+1
        if(p==minn) then
            exit
        end if
    end do
    !____________________________________________________________________________________________________________________
    k=1
    !find length of the hull dummy. Either I'm dumb or the fortran is!!
    do while (.true.)
       if(hull_dummy(k,1) /=0)then
           k=k+1
       else
           exit
       end if
    end do
    !_____________________________________________________________________________________________________________________
    !final hull points are given as
    allocate(hull(k-1,1))
    do i = 1,k-1
        hull(i,1) = hull_dummy(i,1)
    end do
    !______________________________________________________________________________________________________________________
    !optional code to view hull
    !do i=1,size(hull(:,1))
    !    print*,hull(i,:)
    !end do
    !######################################################################################################################    
    !###################################################################################################
    !element file read operation
    open(3,file = 'elements_c.txt',status = 'old')
    read(3,*) ele_d,ele_sur_id,ele_type_id,ele_count
    allocate(elements(ele_count,11))
    do i=1,ele_count
        read(3,*) elements(i,1),elements(i,2),elements(i,3),elements(i,4)
    end do
    do i=1,ele_count
       elements(i,1)=i
    end do
    close(3)
    !element file read operation end
    !##################################################################################################
    !optional code to view elements
    !do i = 1,ele_count
    !    print*,elements(i,:)
    !end do
    !###################################################################################################
    !identify faces
    allocate(faces_fake(ele_count*3,6))
    k=1
    do i=1,ele_count
        faces_fake(k,1) = k
        faces_fake(k,2) = elements(i,2)
        faces_fake(k,3) = elements(i,3)
        faces_fake(k,4:5) = elements(i,1)
        faces_fake(k+1,1) = k+1
        faces_fake(k+1,2) = elements(i,3)
        faces_fake(k+1,3) = elements(i,4)
        faces_fake(k+1,4:5) = elements(i,1)
        faces_fake(k+2,1) = k+2
        faces_fake(k+2,2) = elements(i,2)
        faces_fake(k+2,3) = elements(i,4)
        faces_fake(k+2,4:5) = elements(i,1)
        k = k+3
    end do
    faces_fake(:,6) = 0
    !determining the faces with outer boundaries
    do i=1,size(faces_fake(:,1))
        do j = 1,size(hull(:,1))
            if((hull(j,1)==faces_fake(i,2))) then
                faces_fake(i,6) = faces_fake(i,6)+1
            end if
        end do
        do j = 1,size(hull(:,1))
            if((hull(j,1)==faces_fake(i,3))) then
                faces_fake(i,6) = faces_fake(i,6)+1
            end if
        end do
    end do
    !end determining faces with outer boundaries
    allocate(mask_face(size(faces_fake(:,1))))
    mask_face = 0
    do i = 1,size(faces_fake(:,1))
        do j = i+1,size(faces_fake(:,1))
            if(((faces_fake(i,2)==faces_fake(j,2)).and.(faces_fake(i,3)==faces_fake(j,3))).or.((faces_fake(i,2)==faces_fake(j,3)).and.(faces_fake(i,3)==faces_fake(j,2)))) then
                mask_face(j) = 1
                faces_fake(i,5)=faces_fake(j,4)
            end if            
        end do
    end do
    k=0
    do i = 1,size(mask_face)
        if(mask_face(i) == 0) then
            k = k+1
        end if
    end do
    allocate(faces(k,6))
    j=1
    do i = 1,size(faces_fake(:,1))
        if(mask_face(i)==0) then
            faces(j,:) = faces_fake(i,:)
            j = j+1
        end if
    end do
    do i = 1,size(faces(:,1))
        faces(i,1) = i
    end do
    do i = 1,size(faces(:,1))
        if(faces(i,4)==faces(i,5)) then
            faces(i,5) = 0
        end if
    end do
    print*,size(hull(:,1)), 'nodes form the outer boundary'
    print*,size(faces(:,1)), 'faces detected'
    !face identification end
    !##################################################################################################
    !optional code to view faces
    !do i = 1,size(faces(:,1))
    !    print*,faces(i,:)
    !end do
    !###################################################################################################
    !update element boundary status
    elements(:,5) = 0
    do i=1,ele_count
        do j=1,size(faces(:,1))
            if((faces(j,4)==elements(i,1)).and.(faces(j,5)==0)) then
                elements(i,5) = elements(i,5)+1
            end if
        end do
    end do
    !update element neighbors
    elements(:,6:8) = 0
    k = 6
    do i=1,ele_count
        do j = 1,size(faces(:,1))
            if(faces(j,4)==elements(i,1)) then
                elements(i,k) = faces(j,5)
                k = k+1
                !print*,i,j,k
            end if
            if(faces(j,5)==elements(i,1)) then
                elements(i,k) = faces(j,4)
                k = k+1
                !print*,i,j,k
            end if
        end do
        k = k-3
    end do
    !update element faces
    k=9
    do i=1,ele_count
        do j=1,size(faces(:,1))
            if((faces(j,4)==elements(i,1)).or.(faces(j,5)==elements(i,1))) then
                elements(i,k) = faces(j,1)
                k=k+1
            end if
        end do
        k=k-3
    end do
    print*,size(elements(:,1)),'elements detected'                
    
    !##################################################################################################
    !optional code to view elements
    !do i = 1,size(elements(:,1))
    !    print*,elements(i,:)
    !end do
            
    !##################################################################################################
    !get coordinates for cell centroids
    allocate(centroids(ele_count,3))
    do i=1,ele_count
        centroids(i,1) = elements(i,1)
        centroids(i,2) = (1.0/3.0)*(nodes(elements(i,2),2)+nodes(elements(i,3),2)+nodes(elements(i,4),2))
        centroids(i,3) = (1.0/3.0)*(nodes(elements(i,2),3)+nodes(elements(i,3),3)+nodes(elements(i,4),3))
    end do
    !centroid estimation end
    !##################################################################################################
    !optional code to view cell centroids
    !do i = 1,size(centroids(:,1))
    !    print*,centroids(i,:)
    !end do
    !##################################################################################################
    !get coordinates for face centroids
    allocate(face_centroids(size(faces(:,1)),3))
    do i=1,size(faces(:,1))
        face_centroids(i,1) = faces(i,1)
        face_centroids(i,2) = (1.0/2.0)*(nodes(faces(i,2),2)+nodes(faces(i,3),2))
        face_centroids(i,3) = (1.0/2.0)*(nodes(faces(i,2),3)+nodes(faces(i,3),3))
    end do
    !face centroid estimation end
    !##################################################################################################
    !optional code to view face centroids
    !do i = 1,size(face_centroids(:,1))
    !    print*,face_centroids(i,:)
    !end do
    !##################################################################################################
    !face based data structures
    !delta_eta computation
    allocate(delta_eta(size(faces(:,1)),2))
    do i = 1,size(faces(:,1))
        delta_eta(i,1)=i
        delta_eta(i,2) = sqrt((nodes(faces(i,2),3)-nodes(faces(i,3),3))**2+(nodes(faces(i,2),2)-nodes(faces(i,3),2))**2)
    end do
    !end delta_eta computation
    !#####################################################################################################
    !optional code to view delta_eta
    !do i=1,size(delta_eta(:,1))
    !    print*,delta_eta(i,:)
    !end do
    !#####################################################################################################
    !delta_zhi computation
    allocate(delta_zhi(size(faces(:,1)),2))
    do i = 1,size(faces(:,1))
        delta_zhi(i,1)=i
        if((faces(i,4)/=0).and.faces(i,5)/=0) then
            delta_zhi(i,2) = sqrt((centroids(faces(i,4),3)-centroids(faces(i,5),3))**2+(centroids(faces(i,4),2)-centroids(faces(i,5),2))**2)
        else if(faces(i,4)==0) then
            delta_zhi(i,2) = sqrt((face_centroids(faces(i,1),3)-centroids(faces(i,5),3))**2+(face_centroids(faces(i,1),2)-centroids(faces(i,5),2))**2)
        else if (faces(i,5)==0) then
            delta_zhi(i,2) = sqrt((face_centroids(faces(i,1),3)-centroids(faces(i,4),3))**2+(face_centroids(faces(i,1),2)-centroids(faces(i,4),2))**2)
        end if
    end do
    !end delta_zhi computation
    !#####################################################################################################
    !optional code to view delta_zhi
    !do i=1,size(delta_zhi(:,1))
    !    print*,delta_zhi(i,:)
    !end do
    !#####################################################################################################
    !compute af_dot_af
    allocate(af_dot_af(size(faces(:,1)),2))
    do i=1,size(faces(:,1))
        af_dot_af(i,1) = i
    end do
    af_dot_af(:,2) = delta_eta(:,2)**2
    !end compute af_dot_af
    !########################################################################################################
    !optional code to view af_dot_af
    !do i=1,size(af_dot_af(:,1))
    !    print*,af_dot_af(i,:)
    !end do
    !#####################################################################################################
    !compute af.....which only represents magnitudes of components of face normal
    allocate(af(size(faces(:,1)),3))
    do i=1,size(faces(:,1))
        af(i,1) = i
        af(i,2) = abs(nodes(faces(i,2),3)-nodes(faces(i,3),3))
        af(i,3) = abs(nodes(faces(i,2),2)-nodes(faces(i,3),2))
    end do
    !end compute af
    !########################################################################################################
    !optional code to view af
    !do i=1,size(af(:,1))
    !    print*,af(i,:)
    !end do
    !#########################################################################################################
    !compute e_zhi
    allocate(e_zhi(size(faces(:,1)),3))
    do i=1,size(faces(:,1))
        e_zhi(i,1) = i
        if((faces(i,4)/=0).and.faces(i,5)/=0) then
            e_zhi(i,2) = (centroids(faces(i,4),2)-centroids(faces(i,5),2))/delta_zhi(i,2)
            e_zhi(i,3) = (centroids(faces(i,4),3)-centroids(faces(i,5),3))/delta_zhi(i,2)
        else if(faces(i,4)==0) then
            e_zhi(i,2) = (face_centroids(faces(i,1),2)-centroids(faces(i,5),2))/delta_zhi(i,2)
            e_zhi(i,3) = (face_centroids(faces(i,1),3)-centroids(faces(i,5),3))/delta_zhi(i,2)
        else if (faces(i,5)==0) then
            e_zhi(i,2) = (centroids(faces(i,4),2)-face_centroids(faces(i,1),2))/delta_zhi(i,2)
            e_zhi(i,3) = (centroids(faces(i,4),3)-face_centroids(faces(i,1),3))/delta_zhi(i,2)
        end if
    end do
    !end compute e_zhi
    !########################################################################################################
    !optional code to view e_zhi
    !do i=1,size(e_zhi(:,1))
    !    print*,e_zhi(i,:)
    !end do
    !#########################################################################################################
    !compute af_dot_e_zhi
    allocate(af_dot_e_zhi(size(faces(:,1)),2))
    do i = 1,size(faces(:,1))
        af_dot_e_zhi(i,1) = i
        af_dot_e_zhi(i,2) = abs(af(i,2)*e_zhi(i,2))+abs(af(i,3)*e_zhi(i,3))
    end do
    !end compute af_dot_e_zhi
    !########################################################################################################
    !optional code to view af_dot_e_zhi
    !do i=1,size(af_dot_e_zhi(:,1))
    !    print*,af_dot_e_zhi(i,:)
    !end do
    !#########################################################################################################
    !code to compute the volume of the cell
    allocate(delv(ele_count,2))
    do i=1,ele_count
        delv(i,1)=i
        delv(i,2)=0.5*abs(nodes(elements(i,2),2)*(nodes(elements(i,3),3)-nodes(elements(i,4),3))+nodes(elements(i,3),2)*(nodes(elements(i,4),3)-nodes(elements(i,2),3))+nodes(elements(i,4),2)*(nodes(elements(i,2),3)-nodes(elements(i,3),3)))
    end do
    !########################################################################################################
    !optional code to view delv
    !do i=1,size(delv(:,1))
    !    print*,delv(i,:)
    !end do
    !#########################################################################################################
    !calculate af_vector as seen from every cell centroid
    allocate(af_vector(ele_count,7))
    do i = 1,ele_count
        af_vector(i,1)=i
        
        x2 = nodes(faces(elements(i,9),2),2)
        y2 = nodes(faces(elements(i,9),2),3)
        x1 = nodes(faces(elements(i,9),3),2)
        y1 = nodes(faces(elements(i,9),3),3)
        x0 = centroids(i,2)
        y0 = centroids(i,3)
        if((y2-y1/=0.0).and.(x2-x1/=0.0)) then
            x = ((x0*(x2-x1)/(y2-y1))+(x1*(y2-y1)/(x2-x1))+(y0-y1)) / (((y2-y1)/(x2-x1))+((x2-x1)/(y2-y1)))
            y = y1+((y2-y1)/(x2-x1))*(x-x1)
            af_vector(i,2) =sqrt(((y2-y1)**2+(x2-x1)**2)/((x-x0)**2+(y-y0)**2))*(x-x0)
            af_vector(i,3) = sqrt(((y2-y1)**2+(x2-x1)**2)/((x-x0)**2+(y-y0)**2))*(y-y0)
        else if ((y2-y1==0.0)) then
            if((y2>=y0).and.(y1>=y0)) then
                af_vector(i,2) = 0.0
                af_vector(i,3) = abs(x2-x1)
            else if((y2<=y0).and.(y1<=y0))then
                af_vector(i,2) = 0.0
                af_vector(i,3) = -1.0*abs(x2-x1)
            end if
        else if(x2-x1==0.0) then
            if((x2>=x0).and.(x1>=x0)) then
                af_vector(i,2) = abs(y2-y1)
                af_vector(i,3) = 0.0
            else if((x2<=x0).and.(x1<=x0))then
                af_vector(i,2) = -1.0*abs(y2-y1)
                af_vector(i,3) = 0.0
            end if
        end if
        
            
        
        x2 = nodes(faces(elements(i,10),2),2)
        y2 = nodes(faces(elements(i,10),2),3)
        x1 = nodes(faces(elements(i,10),3),2)
        y1 = nodes(faces(elements(i,10),3),3)
        x0 = centroids(i,2)
        y0 = centroids(i,3)
        if((y2-y1/=0.0).and.(x2-x1/=0.0)) then
            x = ((x0*(x2-x1)/(y2-y1))+(x1*(y2-y1)/(x2-x1))+(y0-y1)) / (((y2-y1)/(x2-x1))+((x2-x1)/(y2-y1)))
            y = y1+((y2-y1)/(x2-x1))*(x-x1)
            af_vector(i,4) =sqrt(((y2-y1)**2+(x2-x1)**2)/((x-x0)**2+(y-y0)**2))*(x-x0)
            af_vector(i,5) = sqrt(((y2-y1)**2+(x2-x1)**2)/((x-x0)**2+(y-y0)**2))*(y-y0)
        else if ((y2-y1==0.0)) then
            if((y2>=y0).and.(y1>=y0)) then
                af_vector(i,4) = 0.0
                af_vector(i,5) = abs(x2-x1)
            else if((y2<=y0).and.(y1<=y0))then
                af_vector(i,4) = 0.0
                af_vector(i,5) = -1.0*abs(x2-x1)
            end if
        else if(x2-x1==0.0) then
            if((x2>=x0).and.(x1>=x0)) then
                af_vector(i,4) = abs(y2-y1)
                af_vector(i,5) = 0.0
            else if((x2<=x0).and.(x1<=x0))then
                af_vector(i,4) = -1.0*abs(y2-y1)
                af_vector(i,5) = 0.0
            end if
        end if
        
        x2 = nodes(faces(elements(i,11),2),2)
        y2 = nodes(faces(elements(i,11),2),3)
        x1 = nodes(faces(elements(i,11),3),2)
        y1 = nodes(faces(elements(i,11),3),3)
        x0 = centroids(i,2)
        y0 = centroids(i,3)
        if((y2-y1/=0.0).and.(x2-x1/=0.0)) then
            x = ((x0*(x2-x1)/(y2-y1))+(x1*(y2-y1)/(x2-x1))+(y0-y1)) / (((y2-y1)/(x2-x1))+((x2-x1)/(y2-y1)))
            y = y1+((y2-y1)/(x2-x1))*(x-x1)
            af_vector(i,6) =sqrt(((y2-y1)**2+(x2-x1)**2)/((x-x0)**2+(y-y0)**2))*(x-x0)
            af_vector(i,7) = sqrt(((y2-y1)**2+(x2-x1)**2)/((x-x0)**2+(y-y0)**2))*(y-y0)
        else if ((y2-y1==0.0)) then
            if((y2>=y0).and.(y1>=y0)) then
                af_vector(i,6) = 0.0
                af_vector(i,7) = abs(x2-x1)
            else if((y2<=y0).and.(y1<=y0))then
                af_vector(i,6) = 0.0
                af_vector(i,7) = -1.0*abs(x2-x1)
            end if
        else if(x2-x1==0.0) then
            if((x2>=x0).and.(x1>=x0)) then
                af_vector(i,6) = abs(y2-y1)
                af_vector(i,7) = 0.0
            else if((x2<=x0).and.(x1<=x0))then
                af_vector(i,6) = -1.0*abs(y2-y1)
                af_vector(i,7) = 0.0
            end if
        end if
    end do
    !end conputing af_vector seen by centroids
    !########################################################################################################
    !optional code to view af_vector
    !do i=1,size(af_vector(:,1))
    !    print*,af_vector(i,:)
    !end do
    !########################################################################################################
    !r from cell to face centroid
    allocate(del_r(size(faces(:,1)),5))
    do i = 1,size(faces(:,1))
        del_r(i,1) = i
        if((faces(i,4)/=0).and.(faces(i,5)/=0)) then
            del_r(i,2) = face_centroids(i,2)-centroids(faces(i,4),2)
            del_r(i,3) = face_centroids(i,3)-centroids(faces(i,4),3)
            del_r(i,4) = face_centroids(i,2)-centroids(faces(i,5),2)
            del_r(i,5) = face_centroids(i,3)-centroids(faces(i,5),3)
        else if(faces(i,4)==0)then
            del_r(i,2) = 0.0
            del_r(i,3) = 0.0
            del_r(i,4) = face_centroids(i,2)-centroids(faces(i,5),2)
            del_r(i,5) = face_centroids(i,3)-centroids(faces(i,5),3)
        else if(faces(i,5)==0)then
            del_r(i,2) = face_centroids(i,2)-centroids(faces(i,4),2)
            del_r(i,3) = face_centroids(i,3)-centroids(faces(i,4),3)
            del_r(i,4) = 0.0
            del_r(i,5) = 0.0
        end if
    end do
    !end del_r calculation
    !#########################################################################################################
    !optional code to view del_r
    !do i=1,size(del_r(:,1))
    !    print*,del_r(i,:)
    !end do
    !########################################################################################################
    !e_zhi, r_zhi calculation as seen from every cell
    allocate(e_zhi_cell(ele_count,7))
    allocate(r_zhi_cell(ele_count,7))
    do i = 1,ele_count
        e_zhi_cell(i,1) = i
        r_zhi_cell(i,1) = i
        if(elements(i,6)/=0) then
            e_zhi_cell(i,2) = (centroids(elements(i,6),2)-centroids(elements(i,1),2))/delta_zhi(elements(i,9),2)
            e_zhi_cell(i,3) = (centroids(elements(i,6),3)-centroids(elements(i,1),3))/delta_zhi(elements(i,9),2)
            r_zhi_cell(i,2) = (centroids(elements(i,6),2)-centroids(elements(i,1),2))
            r_zhi_cell(i,3) = (centroids(elements(i,6),3)-centroids(elements(i,1),3))
        else
            e_zhi_cell(i,2) = (face_centroids(elements(i,9),2)-centroids(elements(i,1),2))/delta_zhi(elements(i,9),2)
            e_zhi_cell(i,3) = (face_centroids(elements(i,9),3)-centroids(elements(i,1),3))/delta_zhi(elements(i,9),2)
            r_zhi_cell(i,2) = (face_centroids(elements(i,9),2)-centroids(elements(i,1),2))
            r_zhi_cell(i,3) = (face_centroids(elements(i,9),3)-centroids(elements(i,1),3))
        end if
        if(elements(i,7)/=0) then
            e_zhi_cell(i,4) = (centroids(elements(i,7),2)-centroids(elements(i,1),2))/delta_zhi(elements(i,10),2)
            e_zhi_cell(i,5) = (centroids(elements(i,7),3)-centroids(elements(i,1),3))/delta_zhi(elements(i,10),2)
            r_zhi_cell(i,4) = (centroids(elements(i,7),2)-centroids(elements(i,1),2))
            r_zhi_cell(i,5) = (centroids(elements(i,7),3)-centroids(elements(i,1),3))
        else
            e_zhi_cell(i,4) = (face_centroids(elements(i,10),2)-centroids(elements(i,1),2))/delta_zhi(elements(i,10),2)
            e_zhi_cell(i,5) = (face_centroids(elements(i,10),3)-centroids(elements(i,1),3))/delta_zhi(elements(i,10),2)
            r_zhi_cell(i,4) = (face_centroids(elements(i,10),2)-centroids(elements(i,1),2))
            r_zhi_cell(i,5) = (face_centroids(elements(i,10),3)-centroids(elements(i,1),3))
        end if
        if(elements(i,8)/=0) then
            e_zhi_cell(i,6) = (centroids(elements(i,8),2)-centroids(elements(i,1),2))/delta_zhi(elements(i,11),2)
            e_zhi_cell(i,7) = (centroids(elements(i,8),3)-centroids(elements(i,1),3))/delta_zhi(elements(i,11),2)
            r_zhi_cell(i,6) = (centroids(elements(i,8),2)-centroids(elements(i,1),2))
            r_zhi_cell(i,7) = (centroids(elements(i,8),3)-centroids(elements(i,1),3))
        else
            e_zhi_cell(i,6) = (face_centroids(elements(i,11),2)-centroids(elements(i,1),2))/delta_zhi(elements(i,11),2)
            e_zhi_cell(i,7) = (face_centroids(elements(i,11),3)-centroids(elements(i,1),3))/delta_zhi(elements(i,11),2)
            r_zhi_cell(i,6) = (face_centroids(elements(i,11),2)-centroids(elements(i,1),2))
            r_zhi_cell(i,7) = (face_centroids(elements(i,11),3)-centroids(elements(i,1),3))
        end if
    end do
    !end e_zhi,r_zhi_cell calculation for every cell
    !############################################################################################################
    !optional code to view e_zhi_cell
    !do i=1,size(e_zhi_cell(:,1))
    !    print*,e_zhi_cell(i,:)
    !end do
    !optional code to view r_zhi_cell
    !do i=1,size(r_zhi_cell(:,1))
    !    print*,r_zhi_cell(i,:)
    !end do
    !#########################################################################################################
    !main solver
    !initialize phi values
    allocate(phi(ele_count,2))
    allocate(phi_dash(ele_count,2))
    allocate(phi_dash_gs(ele_count,2))
    allocate(phi_star(size(faces(:,1)),2))
    phi = 0.0
    
    !set_index
    do i=1,ele_count
        phi(i,1) = i
        phi_dash(i,1) = i
        phi_dash_gs(i,1) = i
    end do
    
    !phi initialization end
    !##################################################################################################
    !optional code to view phi
    !do i=1,ele_count
    !    print*,phi(i,:)
    !end do
    !#################################################################################################
    !calculate gma
    allocate(gma(size(faces(:,1)),2))
    do i=1,size(faces(:,1))
        gma(i,1) = i
        gma(i,2) = 10.0
    end do
    !end gma calculation
    !#####################################################################################################
    !optional code to view gma
    !do i=1,size(gma(:,1))
    !    print*,gma(i,:)
    !end do
    !#####################################################################################################
    !anb_f calculation, defined as attribute of every cell
    allocate(anb_f(ele_count,4))
    do i=1,ele_count
        anb_f(i,1) = i
        anb_f(i,2) = gma(elements(i,9),2)*af_dot_af(elements(i,9),2)/(delta_zhi(elements(i,9),2)*af_dot_e_zhi(elements(i,9),2))
        anb_f(i,3) = gma(elements(i,10),2)*af_dot_af(elements(i,10),2)/(delta_zhi(elements(i,10),2)*af_dot_e_zhi(elements(i,10),2))
        anb_f(i,4) = gma(elements(i,11),2)*af_dot_af(elements(i,11),2)/(delta_zhi(elements(i,11),2)*af_dot_e_zhi(elements(i,11),2))
    end do
    !end anb_f calculation
    !############################################################################################################
    !optional code to view anb_f
    !do i=1,size(anb_f(:,1))
    !    print*,anb_f(i,:)
    !end do
    !#########################################################################################################
    !allocations
    allocate(phi_f(size(faces(:,1)),2))
    allocate(grad_phi(ele_count,3))
    allocate(grad_phi_dash(ele_count,3))
    allocate(grad_phi_f(size(faces(:,1)),3))
    allocate(sgrad(ele_count,4))
    allocate(sc(ele_count,2))
    allocate(sp(ele_count,2))
    allocate(ap(ele_count,2))
    allocate(b(ele_count,2))
    !###########################################################################################################
    !set index
    do i=1,ele_count
        sc(i,1) = i
        sp(i,1) = i
        ap(i,1) = i
        b(i,1) = i
    end do
    !set index
    do i=1,ele_count
        grad_phi(i,1) = i
        grad_phi_dash(i,1) = i
    end do
    !########################################################################################################
    !least squares allocation
    allocate(dummy(2,2),dummy2(2,2),dummy3(2,1),m(3,2),mt(2,3),d(3,1),d_dash(3,1))
    !end least squares allocation
    !########################################################################################################
    
    !dirichlet solver
    !initialize boundary value
    !inner boundary
    phi_f =0.0
    !outer boundary
    do i=1,size(faces(:,1))
        if (faces(i,6)==2) then
            phi_f(i,2) =0.0
        end if
    end do
    
    !neuman solver
    !neuman boundary allocation
    !make sure you satisfy energy conservation if you do not want funny results
    neuman = .true.
    if (neuman) then
    allocate(q_flux(size(faces(:,1)),2))
    !inner boundary
    q_flux = -500.0
    !outer boundary
    do i=1,size(faces(:,1))
        if (faces(i,6)==2) then
            q_flux(i,2) =-500.0
        end if
    end do
    !fixed point must be inside the surface boundary not allowed
    x_q = 0.50
    y_q = 0.50
    !find the centroid closest to this fixed point
    fx=1
    do i=1,ele_count
        if (sqrt((centroids(i,2)-x_q)**2+(centroids(i,3)-y_q)**2)<sqrt((centroids(fx,2)-x_q)**2+(centroids(fx,3)-y_q)**2)) then
            fx = i
        end if
    end do
    !temperature of the fixed point
    phi(fx,2) = 1000.0
    end if
    
    !inner dirichlet outer neuman
    !set your desired inner boundary value to dirichlet BC above
    !make neuman true, make both boundary fluxes the same flux you want at your boundary, 
    !caution ## both inlet and outlet dirichlet cannot be set true here
    inner_dirichlet = .false.
    
    !outer dirichlet inner neuman
    !set your desired outer boundary value to dirichlet BC above
    !make neuman true, make both boundary fluxes the same flux you want at your boundary
    !caution ## both inlet and outlet dirichlet cannot be set true here
    outer_dirichlet = .false.
    
    
    
    !gradient initialization
    grad_phi_f = 0.0
    grad_phi = 0.0
    grad_phi_dash = 0.0
    
    !relaxation
    alpha = 1.0
    
    if (1) then
    !loop start
    count = 0
    !source term loop
    do while(.true.)
        !initialize error for source term
        error = 0.0
        !count to keep track of the iteration
        count = count+1
        print*,count,'outer'
        !backup array
        phi_dash = phi
        !declare source term
        sc(:,2) = 1000.0
        sp(:,2) = 0.0
        ap(:,2) = anb_f(:,2)+anb_f(:,3)+anb_f(:,4)-sp(:,2)*delv(:,2)
        
        if(0)then
        !code to calculate gradient green gauss cell based
        !phi_f is phi value seen at the face
        !update inner face values using linear approximations
        do i=1,size(faces(:,1))
            phi_f(i,1) = i
            if((faces(i,4)/=0).and.(faces(i,5)/=0)) then
                phi_f(i,2) = 0.5*(phi(faces(i,4),2)+phi(faces(i,5),2))
            end if
        end do
        !end computing first phase of cell gradient
        !########################################################################################################
        !optional code to view phi_f
        !do i=1,size(phi_f(:,1))
        !    print*,phi_f(i,:)
        !end do
        !########################################################################################################
        count_grad = 0
        k=0
        !while loop to converge on grad_phi
        do while(.true.)
            !k=k+1
            count_grad = count_grad+1
            print*,count_grad,'grad'
            !backup
            grad_phi_dash = grad_phi
        !#########################################################################################################
        !calculate cell gradient based on the initial phi_f values
            do i=1,ele_count
                grad_phi(i,1)=i
                grad_phi(i,2)=(1.0/delv(i,2))*(phi_f(elements(i,9),2)*af_vector(i,2)+phi_f(elements(i,10),2)*af_vector(i,4)+phi_f(elements(i,11),2)*af_vector(i,6))
                grad_phi(i,3)=(1.0/delv(i,2))*(phi_f(elements(i,9),2)*af_vector(i,3)+phi_f(elements(i,10),2)*af_vector(i,5)+phi_f(elements(i,11),2)*af_vector(i,7))
            end do
        !end first calculation of grad_phi
        !#########################################################################################################
        !optional code to view grad_phi
        !do i=1,size(grad_phi(:,1))
        !    print*,grad_phi(i,:)
        !end do
        !########################################################################################################
        !update phi_f
            do i=1,size(phi_f(:,1))
                if((faces(i,4)/=0).and.(faces(i,5)/=0)) then
                    phi_f(i,2) = (1.0/(norm2(del_r(i,2:3))+norm2(del_r(i,4:5))))*(norm2(del_r(i,4:5))*(phi(faces(i,4),2)+dot_product(grad_phi(faces(i,4),2:3),del_r(i,2:3)))+norm2(del_r(i,2:3))*(phi(faces(i,5),2)+dot_product(grad_phi(faces(i,5),2:3),del_r(i,4:5))))
                end if
            end do
        !end updating phi_f
        !##########################################################################################################
        !optional code to view phi_f
        !do i=1,size(phi_f(:,1))
        !    print*,phi_f(i,:)
        !end do
        !#########################################################################################################
            !print*,sqrt(sum((grad_phi(:,2:3)-grad_phi_dash(:,2:3))**2))
            if (sqrt(sum((grad_phi(:,2:3)-grad_phi_dash(:,2:3))**2))<1e-6) exit
        end do
        !converged to grad_phi
        !##########################################################################################################
        end if
        
        if(0) then
            !gradient calculation using modified green gauss
            !step1 face value calculation
            do i=1,size(faces(:,1))
                phi_f(i,1) = i
                if((faces(i,4)/=0).and.(faces(i,5)/=0)) then
                    phi_f(i,2) = (1.0/(norm2(del_r(i,2:3))+norm2(del_r(i,4:5))))*(norm2(del_r(i,4:5))*phi(faces(i,4),2)+norm2(del_r(i,2:3))*phi(faces(i,5),2))
                end if
            end do
            !step2 cell gradient calculation
             do i=1,ele_count
                grad_phi(i,1)=i
                grad_phi(i,2)=(1.0/delv(i,2))*(phi_f(elements(i,9),2)*af_vector(i,2)+phi_f(elements(i,10),2)*af_vector(i,4)+phi_f(elements(i,11),2)*af_vector(i,6))
                grad_phi(i,3)=(1.0/delv(i,2))*(phi_f(elements(i,9),2)*af_vector(i,3)+phi_f(elements(i,10),2)*af_vector(i,5)+phi_f(elements(i,11),2)*af_vector(i,7))
             end do
        end if
                    
        
        if(1) then
        !gradient calculation using least squares method
        !step1 face value calculation
            do i=1,size(faces(:,1))
                phi_f(i,1) = i
                if((faces(i,4)/=0).and.(faces(i,5)/=0)) then
                    phi_f(i,2) = (1.0/(norm2(del_r(i,2:3))+norm2(del_r(i,4:5))))*(norm2(del_r(i,4:5))*phi(faces(i,4),2)+norm2(del_r(i,2:3))*phi(faces(i,5),2))
                end if
            end do
        do i = 1,ele_count
            !assign matrix m
            m(1,:) = r_zhi_cell(i,2:3)
            m(2,:) = r_zhi_cell(i,4:5)
            m(3,:) = r_zhi_cell(i,6:7)
            !assign matrix d
            d_dash(:,1) = phi_f(elements(i,9:11),2)-phi(elements(i,1),2)
            if(elements(i,6)/=0) then
                d_dash(1,1) = phi(elements(i,6),2)-phi(elements(i,1),2)
            end if
            if(elements(i,7)/=0) then
                d_dash(2,1) = phi(elements(i,7),2)-phi(elements(i,1),2)
            end if
            if(elements(i,8)/=0) then
                d_dash(3,1) = phi(elements(i,8),2)-phi(elements(i,1),2)
            end if
            
            !calculate transpose and inverse
            mt = transpose(m)
            dummy2=matmul(mt,m)
            call inverse(dummy2,dummy,2)
            !find gradient
            dummy3 = matmul(matmul(dummy,mt),d_dash)
            grad_phi(i,2:3) = dummy3(:,1)
        end do
        end if
            
        !finding the value of grad_phi at inner faces, non-uniform
        do i = 1,size(faces(:,1))
            grad_phi_f(i,1) = i
            if((faces(i,4)/=0).and.(faces(i,5)/=0)) then
                grad_phi_f(i,2:3) = (1.0/(norm2(del_r(i,2:3))+norm2(del_r(i,4:5))))*(norm2(del_r(i,4:5))*grad_phi(faces(i,4),2:3)+norm2(del_r(i,2:3))*grad_phi(faces(i,5),2:3))
            else if(faces(i,4)==0)then
                grad_phi_f(i,2:3) = grad_phi(faces(i,5),2:3)
            else if (faces(i,5)==0)then
                grad_phi_f(i,2:3) = grad_phi(faces(i,4),2:3)
            end if
        end do
        !end finding the value of grad_phi at inner faces
        !############################################################################################################
        !optional code to view grad_phi_f
        !do i=1,size(grad_phi_f(:,1))
        !    print*,grad_phi_f(i,:)
        !end do
        !#########################################################################################################
        !secondary gradient calculation viewed as the property of the cell with corresponding neighbors
        do i=1,ele_count
            sgrad(i,1) = i
            sgrad(i,2) = gma(elements(i,9),2)*dot_product(grad_phi_f(elements(i,9),2:3),af_vector(i,2:3))-anb_f(i,2)*dot_product(grad_phi_f(elements(i,9),2:3),e_zhi_cell(i,2:3))*delta_zhi(elements(i,9),2)
            sgrad(i,3) = gma(elements(i,10),2)*dot_product(grad_phi_f(elements(i,10),2:3),af_vector(i,4:5))-anb_f(i,3)*dot_product(grad_phi_f(elements(i,10),2:3),e_zhi_cell(i,4:5))*delta_zhi(elements(i,10),2)
            sgrad(i,4) = gma(elements(i,11),2)*dot_product(grad_phi_f(elements(i,11),2:3),af_vector(i,6:7))-anb_f(i,4)*dot_product(grad_phi_f(elements(i,11),2:3),e_zhi_cell(i,6:7))*delta_zhi(elements(i,11),2)
        end do
        !end secondary gradient calculation
        !########################################################################################################
        !optional code to view sgrad
        !do i=1,size(sgrad(:,1))
        !    print*,sgrad(i,:)
        !end do
        !#########################################################################################################
        !consolidated source term
        b(:,2) = sc(:,2)*delv(:,2)+sgrad(:,2)+sgrad(:,3)+sgrad(:,4)
        !inner count start
        count_in = 0
        !inner GS loop
        do while(.true.)
            !initialize error for gs iteration
            error_in = 0.0
            !count to keep track of the iteration
            count_in = count_in +1
            print*,count_in, 'inner'
            !backup array
            phi_dash_gs = phi
            phi_star = phi_f
            !update values using gauss seidel iterations
            do i = 1,ele_count
                if((neuman).and.(i==fx).and.(inner_dirichlet/=.true.).and.(outer_dirichlet/=.true.)) cycle
                d(:,1) = phi_f(elements(i,9:11),2)
                if(elements(i,6)/=0) then
                    d(1,1) = phi(elements(i,6),2)
                end if
                if(elements(i,7)/=0) then
                    d(2,1) = phi(elements(i,7),2)
                end if
                if(elements(i,8)/=0) then
                    d(3,1) = phi(elements(i,8),2)
                end if
                phi(i,2) = (1.0/ap(i,2))*(anb_f(i,2)*d(1,1)+anb_f(i,3)*d(2,1)+anb_f(i,4)*d(3,1)+b(i,2))           
            end do
            !update boundary values only for the neuman iterations
            if(neuman.and.(inner_dirichlet/=.true.).and.(outer_dirichlet/=.true.)) then    
                do i=1,size(faces(:,1))
                    if(faces(i,5)==0) then
                        phi_f(i,2) = (q_flux(i,2)+(gma(i,2)/delta_zhi(i,2))*phi(faces(i,4),2))/(gma(i,2)/delta_zhi(i,2))
                    end if
                end do
            end if
            !update boundary values only for inner dirichlet condition
            if(neuman.and.(inner_dirichlet==.true.).and.(outer_dirichlet/=.true.)) then
                do i=1,size(faces(:,1))
                    if((faces(i,5)==0).and.(faces(i,6)==2)) then
                        phi_f(i,2) = (q_flux(i,2)+(gma(i,2)/delta_zhi(i,2))*phi(faces(i,4),2))/(gma(i,2)/delta_zhi(i,2))
                    end if
                end do
            end if
            !update boundary values only for outer dirichlet condition
            if(neuman.and.(inner_dirichlet/=.true.).and.(outer_dirichlet==.true.)) then
                do i=1,size(faces(:,1))
                    if((faces(i,5)==0).and.(faces(i,6)/=2)) then
                        phi_f(i,2) = (q_flux(i,2)+(gma(i,2)/delta_zhi(i,2))*phi(faces(i,4),2))/(gma(i,2)/delta_zhi(i,2))
                    end if
                end do
            end if
            !relaxation
            if(neuman.and.(inner_dirichlet/=.true.).and.(outer_dirichlet/=.true.)) then
                phi(:,2) = alpha*phi(:,2)+(1.0-alpha)*phi_dash_gs(:,2)
                phi_f(:,2) = alpha*phi_f(:,2)+(1.0-alpha)*phi_star(:,2)
            end if
            !error
            error_in = sqrt(sum((phi(:,2)-phi_dash_gs(:,2))**2))
            if(error_in<1e-6) exit
        end do
        !end inner GS loop
        error = sqrt(sum((phi(:,2)-phi_dash(:,2))**2))
        if(error<1e-6) exit
    end do
    !end outer loop
    
    end if
    !########################################################################################################################
    !optional code to view phi
    !do i=1,ele_count
    !    print*,phi(i,:)
    !end do
    !############################################################################################################################
    !check contribution of the secondary gradient
    !do i = 1,ele_count
    !    print*, i,sgrad(i,2)+sgrad(i,3)+sgrad(i,4)
    !end do
    !################################################################################################################################
    !optional code to view grad_phi
    !do i=1,size(grad_phi(:,1))
    !    print*,grad_phi(i,:)
    !end do
    !###################################################################################################################################
    !file output
    !write phi values
    open(unit= 4,file='phi.txt',access='sequential',action='write')
    do i = 1,ele_count
        write(4,*) centroids(i,2),centroids(i,3),phi(i,2)
    end do
    close(4)
    !end file write operation
    !###################################################################################################################################
    !take sample output on a line horizontal at the middle of the domain
    allocate(line_fake(ele_count,3))
    y_avg = 0.5*(minval(centroids(:,3))+maxval(centroids(:,3)))
    tolerance_avg = sum(delta_eta(:,2))/size(delta_eta(:,1))
    k=0
    !print*,y_avg,tolerance_avg
    do i = 1,ele_count
        if((centroids(i,3)<=(y_avg+tolerance_avg)).and.(centroids(i,3)>=(y_avg-tolerance_avg))) then
            k=k+1
            line_fake(k,1)=centroids(i,2)
            line_fake(k,2)=y_avg
            line_fake(k,3)=phi(i,2)
        end if
    end do
    allocate(line(k,3))
    do i = 1,k
        line(i,:) = line_fake(i,:)
    end do
    open(unit= 5,file='line.txt',access='sequential',action='write')
    do i = 1,k
        write(5,*) line(i,:)
    end do
    close(5)
            
    
    
    
    deallocate(nodes,elements,faces_fake,faces,centroids,face_centroids,phi,gma,delta_eta,delta_zhi,af_dot_af,af,e_zhi,af_dot_e_zhi,phi_f,delv,af_vector,grad_phi,del_r,grad_phi_dash,grad_phi_f,anb_f,e_zhi_cell,sgrad,phi_dash,sc,sp,ap,phi_dash_gs,b,r_zhi_cell,dummy,m,mt,d,dummy2,dummy3,d_dash,hull_dummy,hull,line_fake,line,phi_star)
    if(neuman) deallocate(q_flux)
end program reading
    