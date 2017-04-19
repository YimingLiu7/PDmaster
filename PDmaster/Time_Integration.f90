!***************************************************************
!                    Time_Integration                          *
!          功能： 施加等效体积密度荷载                         *
!***************************************************************    
subroutine Time_Integration()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 dens,cv,dt
    cv=solve(iload)%cv
    dt=solve(iload)%dt
!***计算位移
    do idime=1,ndime
       do ni=1,totnode
          disp(ni,idime)=disp(ni,idime)+velhalfold(ni,idime)*dt
       end do
    end do
!***施加位移荷载
    call Applied_Displacement      
!***计算体力
    call internal_force
!***施加体力荷载
    call applied_force()
!***计算接触力
!    call contect_force
!***计算合力force   
    do idime=1,ndime
       do ni=1,totnode
          force(ni,idime)=pforce(ni,idime)+bforce(ni,idime)!+cforce(ni,idime)+gforce(ni,idime)
       end do
    end do
!***计算速度velhalf
    do ni=1,totnode 
       dens=mats(matsnode(ni))%dens
       do idime=1,ndime
          !velhalf(ni,idime)=(force(ni,idime)+(dens/dt-cv/2)*velhalfold(ni,idime))/(dens/dt+cv/2)
          velhalf(ni,idime)=(2*dt*force(ni,idime)+(2*dens-cv*dt)*velhalfold(ni,idime))/(2*dens+cv*dt)
          vel(ni,idime)=(velhalf(ni,idime)+velhalfold(ni,idime))/2       
          acc(ni,idime)=(velhalf(ni,idime)-velhalfold(ni,idime))/dt
          velhalfold(ni,idime)=velhalf(ni,idime)
       end do
    end do
!***施加速度荷载velhalf
    call Applied_Velocity()      

end subroutine Time_Integration
    
    
    
    
!***************************************************************
!                  check_convergence                           *
!          功能： 施加等效体积密度荷载                         *
!***************************************************************    
subroutine check_convergence()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 fsq,bfsq,ptoler,toler
    toler=solve(iload)%toler
    do idime=1,ndime
        do ni=1,totnode 
            if(nodeyes(ni).eq.0)cycle
            !位移加载点1，速度加载点2不计算
            if (typenode(iload,ni)==1.or.typenode(iload,ni)==2)cycle 
            
	  	    fsq=fsq+force(ni,idime)*force(ni,idime)
            bfsq=bfsq+bforce(ni,idime)*bforce(ni,idime)+gforce(ni,idime)*gforce(ni,idime)
        end do
    end do
	ptoler=sqrt(fsq)/sqrt(bfsq)
    if(ptoler<toler) check_c=1
    write(7,'(1x,i2,1x,i3,1x,i4,1x,e12.5,1x,e12.5,1x,e12.5)')iload,iincs,istep,ptoler,fsq,bfsq
end subroutine check_convergence      
       
!***************************************************************
!                  check_damage                                *
!          功能： 施加等效体积密度荷载                         *
!***************************************************************    
subroutine check_damage()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 s,st,sc,fac
    integer i,j,numbond
!***判断键损伤
    do ni=1,totnode 
       !判断物质点ni是否被剔除
       if (nodeyes(ni)==0)cycle
       !判断物质点ni是否为边界点或边界领域
       if(solve(iload)%loadtype==1)then
           if (typenode(iload,ni).ne.0)cycle 
       endif
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          !判断是否有加载域
          if(solve(iload)%loadtype==1)then
             !加载点123以及加载域-1标记的点不破坏
             if (typenode(iload,nj).ne.0)cycle
          endif
          s=efam(pointfam(ni)+j-1)/idistfam(pointfam(ni)+j-1)
          
          !设置材料界面系数
          fac=1.0d0
          if(matsnode(ni).ne.matsnode(nj)) fac=1.2d0
          
          st=min((mats(matsnode(ni))%st),(mats(matsnode(nj))%st))*fac
          sc=min((mats(matsnode(ni))%sc),(mats(matsnode(nj))%sc))*fac
          if (s>st.or.s<-sc)then
             nodefam(pointfam(ni)+j-1)=0
          end if
       end do
    end do
    
!***计算物质点损伤    
    do ni=1,totnode 
       if (nodeyes(ni)==0)cycle
       numbond=0
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj.ne.0) numbond=numbond+1
       end do
       damage(ni)=1.0-dble(numbond)/dble(numfam(ni))
    end do
end subroutine check_damage       