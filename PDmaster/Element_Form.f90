!***************************************************************
!                   element_form                               *
!          功能：                                              *
!***************************************************************   
subroutine element_form()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer::judge_result,fjudge_crack,i,j,k
    numfam=0
    pointfam=0
    nodefam=0
    idistfam=0.0d0
    efam=0.0d0
    edfam=0.0d0
    
    cnumfam=0
    cpointfam=0
    cnodefam=0
    cidistfam=0.0d0
    delta=(n_delta+0.01d0)*dx
!***确定近场范围物质点    
    do ni = 1,totnode     
       if(ni==1)then
           pointfam(ni)=1
           cpointfam(ni)=1 
       else
           pointfam(ni)=pointfam(ni-1)+numfam(ni-1)
           cpointfam(ni)=cpointfam(ni-1)+cnumfam(ni-1)
       end if
       if(nodeyes(ni).eq.0)cycle
	   do nj = 1,totnode
           if(nodeyes(nj).eq.0)cycle
           !计算物质点间距
           call sub_idist()
           !判断是否在近场范围内
           if(ni.ne.nj)then
               if(idist<=delta)then 
                   
                   cnumfam(ni)=cnumfam(ni)+1
                   cnodefam(cpointfam(ni)+cnumfam(ni)-1)=nj
                   cidistfam(cpointfam(ni)+cnumfam(ni)-1)=idist
                   
                   judge_result=fjudge_crack(coord(ni,1:2),coord(nj,1:2))
                   judge_result=1
                   if(judge_result==1)then
                      numfam(ni)=numfam(ni)+1
                      nodefam(pointfam(ni)+numfam(ni)-1)=nj
                      idistfam(pointfam(ni)+numfam(ni)-1)=idist
                   end if
               end if
           end if
       end do
    end do      
end subroutine element_form
    
subroutine celement_form()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer::judge_result,fjudge_crack,i,j,k 
    real*8 delx
    cnumfam=0
    cpointfam=0
    cnodefam=0
    cidistfam=0.0d0
!***确定近场范围物质点    
    do ni = 1,totnode     
       if(ni==1)then
           cpointfam(ni)=1 
       else
           cpointfam(ni)=cpointfam(ni-1)+cnumfam(ni-1)
       end if
       if(nodeyes(ni).eq.0)cycle
	   do nj = 1,totnode
           if(nodeyes(nj).eq.0)cycle
           !计算物质点间距
          nlength=0.0d0
          do idime=1,ndime
             delx=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx*delx
          end do
          nlength = dsqrt(nlength)
           
           !判断是否在近场范围内
           if(ni.ne.nj)then
               if(nlength<=delta)then 
                   cnumfam(ni)=cnumfam(ni)+1
                   cnodefam(cpointfam(ni)+cnumfam(ni)-1)=nj
                   cidistfam(cpointfam(ni)+cnumfam(ni)-1)=idist
               end if
           end if
       end do
    end do      
end subroutine celement_form    
    
    
    
!***************************************************************
!                   sub_idist                                  *
!          功能：                                              *
!***************************************************************        
subroutine sub_idist()
    use scalars
    use serialArrays
    use globalParameters
    implicit none  
!***计算物质点间距
    if (ndime==1)then
       idist = dabs(coord(nj,1)-coord(ni,1))   
    elseif(ndime==2)then
       idist = dsqrt((coord(nj,1)-coord(ni,1))**2 + (coord(nj,2)-coord(ni,2))**2)
    else
       idist = dsqrt((coord(nj,1)-coord(ni,1))**2 + (coord(nj,2)-coord(ni,2))**2+(coord(nj,3)-coord(ni,3))**2)
    endif
!***
end subroutine sub_idist
    
    
!***************************************************************
!                   fjudge_crack                               *
!          功能：                                              *
!***************************************************************       
function fjudge_crack(coord1,coord2) result(judge_result)
    use scalars
    use serialArrays  
    use globalParameters
    implicit none 
    real(8)::coord1(1,2),coord2(1,2)
    real(8)::min_x,min_y,max_x,max_y,l,di,dj
    real(8)::delta1,delta2,a0x,a0y,b0x,b0y,c0x,c0y,ax,ay,bx,by,cx,cy
    integer::judge_result
    judge_result=1
    do icrack=1,ncrack
            !判断ipart，jpart是否在线段上  
            min_y=min(cracks(icrack,2),cracks(icrack,4))
            max_y=max(cracks(icrack,2),cracks(icrack,4))
            min_x=min(cracks(icrack,1),cracks(icrack,3))
            max_x=max(cracks(icrack,1),cracks(icrack,3))
            l=sqrt((cracks(icrack,1)-cracks(icrack,3))**2+(cracks(icrack,2)-cracks(icrack,4))**2)
            !i点到线段的距离
            di=abs((coord1(1,1)-cracks(icrack,1))*(coord1(1,2)-cracks(icrack,4))-(coord1(1,2)-cracks(icrack,2))*(coord1(1,1)-cracks(icrack,3)))/l
            if(di<=0.01*dx.and.(coord1(1,1)>=min_x.and.max_x>=coord1(1,1)).and.(coord1(1,2)>=min_y.and.max_y>=coord1(1,2))) then
               judge_result=0
               exit
            end if
            !j点到线段的距离
            dj=abs((coord2(1,1)-cracks(icrack,1))*(coord2(1,2)-cracks(icrack,4))-(coord2(1,2)-cracks(icrack,2))*(coord2(1,1)-cracks(icrack,3)))/l
            if(dj<=0.01*dx.and.(coord2(1,1)>=min_x.and.max_x>=coord2(1,1)).and.(coord2(1,2)>=min_y.and.max_y>=coord2(1,2))) then
               judge_result=0
               exit
            end if
               
            !判断i-j连线是否与裂纹线相交
            a0x=coord2(1,1)-coord1(1,1)
	        a0y=coord2(1,2)-coord1(1,2)
	        b0x=cracks(icrack,1)-coord1(1,1)
	        b0y=cracks(icrack,2)-coord1(1,2)
	        c0x=cracks(icrack,3)-coord1(1,1)
	        c0y=cracks(icrack,4)-coord1(1,2)
	        ax=cracks(icrack,3)-cracks(icrack,1)
            ay=cracks(icrack,4)-cracks(icrack,2)
            bx=coord1(1,1)-cracks(icrack,1)
	        by=coord1(1,2)-cracks(icrack,2)
            cx=coord2(1,1)-cracks(icrack,1)
	        cy=coord2(1,2)-cracks(icrack,2)
            delta1=(a0x*b0y-a0y*b0x)*(a0x*c0y-a0y*c0x)
            delta2=(ax*by-ay*bx)*(ax*cy-ay*cx)
            if(delta1<0.and.delta2<0) then
	           judge_result=0
            endif
    enddo
end function fjudge_crack    