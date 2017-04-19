!    integer,allocatable::bcyes(:,:),bcnode(:,:),bcvalue(:,:)
!    integer,allocatable::numbc_d(:,:),pointbc_d(:,:),numbc_v(:,:),pointbc_v(:,:),numbc_f(:,:),pointbc_f(:,:)

!***************************************************************
!                Applied_Displacement                          *
!          功能：   施加位移荷载                               *
!***************************************************************    
subroutine Applied_Displacement()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer::i,j
    if(bcyes(iload,1)==1)then
        do idime=1,ndime
            if (numbc_d(iload,idime)>0)then
                do i=1,numbc_d(iload,idime)
                    j=bcnode(iload,pointbc_d(iload,idime)+i-1)
                    disp(j,idime)=bcvalue(iload,pointbc_d(iload,idime)+i-1)*loadg(iload,iincs)

                end do
            end if
        end do
    end if
    
end subroutine Applied_Displacement

!***************************************************************
!                   Applied_Velocity                           *
!          功能：    施加速度荷载                              *
!***************************************************************   
subroutine Applied_Velocity()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer::i,j
    if(bcyes(iload,2)==1)then
        do idime=1,ndime
            if (numbc_v(iload,idime)>0)then
                do i=1,numbc_v(iload,idime)
                    j=bcnode(iload,pointbc_v(iload,idime)+i-1)
                    velhalfold(j,idime)=bcvalue(iload,pointbc_v(iload,idime)+i-1)*loadg(iload,iincs)
                end do
            end if
        end do
    end if
    
end subroutine Applied_Velocity      
!***************************************************************
!             Applied_Force                                    *
!          功能： 施加等效体积密度荷载                         *
!***************************************************************    
subroutine Applied_Force()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer::i,j
    integer loadid
    if(bcyes(iload,3)==1)then
        do idime=1,ndime
            if(numbc_f(iload,idime)>0)then
                do i=1,numbc_f(iload,idime)
                    j=bcnode(iload,pointbc_f(iload,idime)+i-1)
                    bforce(j,idime)=bcvalue(iload,pointbc_f(iload,idime)+i-1)*loadg(iload,iincs)

                end do
            end if
        end do
    end if
    !施加历时荷载
    !if(iload>1)then
    !    do loadid=1:iload-1
    !        if (bcyes(loadid,3)==1)then
    !            do i=1,numbc_f(iload,idime)
    !                j=bcnode(loadid,pointbc_f(loadid,idime)+i-1)
    !                bforce(j,idime)=bforce(j,idime)+bcvalue(loadid,pointbc_f(loadid,idime)+i-1)*loadg(loadid,solve(loadid)%nincs)
    !            end do
    !        end if
    !    end do
    ! end if
end subroutine Applied_Force
    

    
  