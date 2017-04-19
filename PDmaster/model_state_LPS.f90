!***************************************************************
!                    model_state_LPS                           *
!          功能：                                              *
!***************************************************************         
subroutine model_state_LPS()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 kmod,gmod,v,k,alpha,gama,t,delx(3)
    integer i,j
!!***测试
!    do idime=1,ndime
!        do ni=1,totnode
!           disp(ni,idime)=0.001*coord(ni,idime)
!        end do
!    end do    
   
!***计算物质点dilatation
    theta=0.0d0
    do ni=1,totnode 
       if (nodeyes(ni)==0)cycle
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          !计算ni和nj距离
          if (nj==0)cycle 
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !读取ni和nj初始距离
          idist= idistfam(pointfam(ni)+j-1)
          !计算变形量
          e=nlength-idist
          efam(pointfam(ni)+j-1)=e
          !计算体积修正系数,并修正长度
          call sub_v_fac()
          call sub_omega()
          theta(ni)=theta(ni)+omega*idist*e*vol(nj)*v_fac
       end do
       if (ndime==3) then
           theta(ni)=3.0d0/wvolume(ni)*theta(ni)        !3D
       elseif (ndime==2)then
          v=mats(matsnode(ni))%v 
          
          if (planetype=='stress')then
             theta(ni)=(4.0d0*v-2.0d0)/(v-1.0d0)/wvolume(ni)*theta(ni) !2D,stress
          elseif(planetype=='strain')then
             theta(ni)=2.0d0/wvolume(ni)*theta(ni)   !2D,strain
          else
              write(*,*)planetype,'平面问题类型出错！！！'
              stop
          endif
             
       else
          !1D
       end if
    end do
    
    
!***初始化
    pforce=0.0d0
!***计算物质点内力密度
    do ni=1,totnode 
       if (nodeyes(ni)==0)cycle
       
       !计算模型参数
       kmod=mats(matsnode(ni))%kmod
       gmod=mats(matsnode(ni))%gmod
       v=mats(matsnode(ni))%v 
       !计算微模量
       if (ndime==3)then
           !3D
          k=kmod
          alpha=15.0d0*gmod
          gama=3.0d0
       elseif(ndime==2)then
          if (planetype=='stress')then
             !2D stress
             k=kmod+gmod/9.0d0*((1.0d0+v)/(2.0d0*v-1.0d0))**2 
             alpha=8.0d0*gmod/wvolume(ni)
             gama=(4.0d0*v-2.0d0)/(v-1.0d0)
          elseif(planetype=='strain')then
             !2D strain
             k=kmod+gmod/9.0d0
             alpha=8.0d0*gmod/wvolume(ni)
             gama=2.0d0
          else
              write(*,*)planetype,'平面问题类型出错！！！'
              stop
          endif
       else 
            !1D
       endif 
       
       !直接读取模型参数
       !k=mods(matsnode(ni))%k
       !alpha=mods(matsnode(ni))%alpha
       !gama=mods(matsnode(ni))%gama
       
       !遍历ni领域所有物质点
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle      
          !计算现时构形物质点间距
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength= nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          
          idist= idistfam(pointfam(ni)+j-1)
          !计算拉伸标量
          e=nlength-idist 
          !计算拉伸标量偏量部分
          !ed=e-theta(ni)*idist/3
          call sub_v_fac() 
          call sub_omega()
          !计算力标量状态
          !简化前
          !t=(gama*k/wvolume(ni)-alpha/3*(1-gama/3))*theta(ni)*omega*idist+alpha*omega*ed
          !简化后
          t=(gama*k/wvolume(ni)-alpha/3*(2-gama/3))*theta(ni)*omega*idist+alpha*omega*e
          
          !if(istep==1) write(8,*)ni,nj,t,e,k,gama,wvolume(ni),theta(ni)
          
          do idime=1,ndime
             pforce(ni,idime)=pforce(ni,idime)+t*delx(idime)/nlength*vol(nj)*v_fac
             pforce(nj,idime)=pforce(nj,idime)-t*delx(idime)/nlength*vol(ni)*v_fac
          end do
       end do
    end do
end subroutine model_state_LPS 