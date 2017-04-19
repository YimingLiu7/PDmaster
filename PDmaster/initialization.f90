!***************************************************************
!                   initialization                             *
!          功能：变量初始化赋值                                *
!***************************************************************      
subroutine initialization()
    use scalars
    use serialArrays
    use globalParameters
    implicit none 
    integer::i,j
    real(8)::emod,v,kmod,gmod,sum,maxdt,ceff
    !近场范围
    force=0.0d0
    pforce=0.0d0
    pforceold=0.0d0
    bforce=0.0d0
    cforce=0.0d0
    gforce=0.0d0
    disp=0.0d0
    vel=0.0d0
    velhalf=0.0d0
    velhalfold=0.0d0
    acc=0.0d0
    damage=0.0d0
!***对速度边界条件初始化
    if(bcyes(1,2)==1)then
        do idime=1,ndime
            do i=1,numbc_v(iload,idime)
                j=bcnode(iload,pointbc_v(iload,idime)+i-1)
                velhalfold(j,idime)=bcvalue(iload,pointbc_v(iload,idime)+i-1)
            end do
        end do
    end if
!!***施加自重荷载
!    do ni=1,totnode
!       gforce(ni,2)=-g*mats(matsnode(ni))%dens
!    end do
    
!***计算weighted volume
    if(model=='LPS')then
        do ni=1,totnode
           wvolume(ni)=0.0d0
           if (nodeyes(ni)==0)cycle
           do j=1,numfam(ni)
              nj=nodefam(pointfam(ni)+j-1)
              idist=idistfam(pointfam(ni)+j-1)
              call sub_v_fac() 
              call sub_omega()
              wvolume(ni)=wvolume(ni)+omega*idist*idist*vol(nj)*v_fac
           end do
        end do
    end if
    
!***计算最大时间步长
        !maxdt=100
        !do ni=1,totnode
        !   if (nodeyes(ni)==0)cycle
        !   sum=0
        !   ceff=9*mats(matsnode(ni))%emod/(pi*delta**3)
        !   do j=1,numfam(ni)
        !      nj=nodefam(pointfam(ni)+j-1)
        !      idist=idistfam(pointfam(ni)+j-1)
        !      call sub_v_fac() 
        !      sum=sum+ceff/idist*vol(nj)*v_fac
        !   end do
        !   sum=sqrt(2*mats(matsnode(ni))%dens/sum)
        !   if(sum<maxdt)maxdt=sum
        !end do
        !write(8,*)maxdt
    
    
    
    
    
    
!!***LPS材料的微模量  
!    if(model=='LPS')then
!        do i=1,nmats
!           kmod=mats(i)%kmod
!           gmod=mats(i)%gmod
!           v=mats(i)%v 
!           !计算微模量
!           if (ndime==3)then
!               !3D
!              mods(i)%k=kmod
!              mods(i)%alpha=15.0d0*gmod
!              mods(i)%gama=3.0d0
!           elseif(ndime==2)then
!              !2D stress
!              mods(i)%k=kmod+gmod/9.0d0*((1.0d0+v)/(2.0d0*v-1.0d0))**2 
!              mods(i)%alpha=8.0d0*gmod/wvolume(ni)
!              mods(i)%gama=(4.0d0*v-2.0d0)/(v-1.0d0)
!              !!2D strain
!              !mods(i)%k=kmod+gmod/9.0d0
!              !mods(i)%alpha=8.0d0*gmod/wvolume(ni)
!              !mods(i)%gama=2.0d0
!           else 
!                !1D
!           endif
!        end do
!    end if
!    
!!***PMB材料的微模量  
!    if(model=='PMB')then
!        do i=1,nmats
!           emod=mats(i)%emod
!           !v=mats(i)%v   
!           if (ndime==3) then !c=18*k/(pi*delta**4)
!               v=1/4
!               mods(i)%c=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
!           elseif (ndime==2) then   !c=12*k/(pi*delta**3)
!               v=1/3
!               mods(i)%c=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
!               !v=1/4
!               !mods(i)%c=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
!           else    !c=2E/delta^2
!               mods(i)%c=2*emod/(pi*delta**2)
!           end if 
!        end do
!    end if
!    
!!***CDG材料的微模量  
!    if(model=='CDG')then
!        do i=1,nmats
!           emod=mats(i)%emod
!           v=mats(i)%v 
!           if (ndime==3) then 
!               !3D
!               mods(i)%c=6*emod/(pi*delta**4*(1-2*v))
!               mods(i)%d=3*emod*(1-4*v)/(pi*delta*delta*(1-2*v)*(1+v))
!           elseif (ndime==2) then
!               !Plane stress
!               mods(i)%c=6*emod/(pi*delta**3*(1-v))
!               mods(i)%d=emod*(1-3*v)/(6*pi*delta*(1-v*v))
!               !!Plane strain
!               !mods(i)%c=6*emod/(pi*delta**3*(1-2*v)*(1+v))
!               !mods(i)%d=emod*(1-4*v)/(6*pi*delta*(1-2*v)*(1+v))
!           else   
!               mods(i)%c=2*emod/(pi*delta**2)
!               mods(i)%d=0.0d0
!           end if
!        end do
!    end if
    
end subroutine initialization
!***************************************************************
!                  sub_omega                                   *
!          功能：计算影响函数                                  *
!***************************************************************       
subroutine sub_omega()
    use scalars
    use serialArrays
    use globalParameters
    implicit none  
!***Gaussian influence function    
    !omega=exp(-idist*idist/(delta*delta))
!***
    omega=1.0d0
end subroutine sub_omega 
!***************************************************************
!                   sub_v_fac                                  *
!          功能：计算体积修正系数                              *
!***************************************************************      
subroutine sub_v_fac()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    if(idist<=delta-dx/2)then
        v_fac=1.0d0
    elseif(idist<=delta)then
        v_fac=(delta+dx/2-idist)/dx
        !idist=idist-dx*(1-v_fac)/2
    elseif(idist<=delta+dx/2)then
        v_fac=(delta+dx/2-idist)/dx
        !idist=delta-dx*v_fac/2
    else
        v_fac=0.0d0
    end if
    !v_fac=1.0d0
end subroutine sub_v_fac  
    
!***************************************************************
!                   sub_surface_fac                            *
!          功能：计算体积修正系数                              *
!***************************************************************          
subroutine sub_surface_fac ()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,Wel
    do idime=1,ndime
        do ni=1,totnode
           disp(ni,idime)=0.001*coord(ni,idime)
        end do
    end do
    
    if (model=='PMB')then
       call sub_Strain_Energy_PMB
    !elseif(model=='LPS')then
       call model_state_LPS()
    !elseif(model=='CDG')then
       !call model_bond_CDG()
    else
       write(*,*)model,'材料模型类型出错！！！'
       stop
    end if  
    do ni=1,totnode 
    !计算ni点应变能密度积分值及修正系数
       emod=mats(matsnode(ni))%emod
       v=mats(matsnode(ni))%v   
       if (ndime==3) then 
           Wel=emod*0.001*0.001
       elseif (ndime==2) then  
           if (planetype=='stress')then
               Wel=emod*0.001*0.001/(1-v)
           elseif(planetype=='strain')then
               Wel=emod*0.001*0.001
            else
               write(*,*)planetype,'平面问题类型出错！！！'
               stop
            endif
       else   
           Wel=emod*0.001*0.001
       end if
       s_fac(ni)=Wel/Strain_Energy(ni)
    end do
end subroutine sub_surface_fac  


