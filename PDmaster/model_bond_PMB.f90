!***************************************************************
!                    model_bond_PMB                            *
!          功能：                                              *
!***************************************************************    
subroutine model_bond_PMB
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,c,ci,cj,delx(3),pf
    integer i,j
!***初始化
    pforce=0.0d0
!***计算ni点的内力密度
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       
       !计算ni点的微模量,注意PMB泊松比为固定值
       emod=mats(matsnode(ni))%emod
       !v=mats(matsnode(ni))%v   
       if (ndime==3) then !c=18*k/(pi*delta**4)
           !v=1/4
           !ci=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
           
           !简化后
           ci=12*emod/(pi*delta**4)
       elseif (ndime==2) then   !c=12*k/(pi*delta**3)
           if (planetype=='stress')then
               !v=1/3
               !ci=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
               !简化后
               ci=9*emod/(pi*delta**3)!Plane stress
            elseif(planetype=='strain')then
               !!v=1/4
               !ci=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
               !简化后
               ci=48/5*emod/(pi*delta**3)!Plane strain
            else
               write(*,*)planetype,'平面问题类型出错！！！'
               stop
            endif
       else    !c=2E/delta^2
           ci=2*emod/(pi*delta**2)
       end if 
       
       !读取ni点的微模量
       !ci=mods(matsnode(ni))%c

       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          
          !计算nj点的微模量，PMB模型中泊松比为固定值
          emod=mats(matsnode(nj))%emod
          !v=mats(matsnode(nj))%v          
          if (ndime==3) then !c=18*k/(pi*delta**4)
              !v=1/4
              !cj=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
              cj=12*emod/(pi*delta**4)
          elseif (ndime==2) then   !c=12*k/(pi*delta**3)
              if (planetype=='stress')then
                 !v=1/3
                 !cj=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
                 !简化后
                 cj=9*emod/(pi*delta**3)!Plane stress
               elseif(planetype=='strain')then
                 !v=1/4
                 !cj=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
                 !简化后
                 cj=48/5*emod/(pi*delta**3)!Plane strain
               else
                  write(*,*)planetype,'平面问题类型出错！！！'
                  stop
                endif
          else    !c=2E/delta^2
              cj=2*emod/(pi*delta**2)
          end if  
          
          !读取nj点的微模量
          !cj=mods(matsnode(nj))%c
          
          !当键两端材料不同时，微模量取串联值
          c=2/(1/ci+1/cj)
          !计算键长
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !读取初始键长
          idist= idistfam(pointfam(ni)+j-1)
          !计算键伸长量并存储
          e=nlength-idist
          efam(pointfam(ni)+j-1)=e
          !修正i领域内j点积分体积Vj         
          call sub_v_fac()
          pf=c*e/idist*vol(nj)*v_fac
          
          !if(istep==50)write(8,*)ci,cj,e
          !if(istep==50)write(8,*)c,c*e/idist*delx(1)/nlength,c*e/idist*delx(2)/nlength,pf
          
          do idime=1,ndime
             pforce(ni,idime)=pforce(ni,idime)+pf*delx(idime)/nlength
          end do
          
         ! if(istep==50)write(8,*)pforce(ni,1),pforce(ni,2)
          
       end do
    end do
end subroutine model_bond_PMB
    
    
subroutine sub_Strain_Energy_PMB
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,c,ci,cj,delx(3),pf
    integer i,j
!***初始化
    Strain_Energy=0.0d0
!***计算ni点的应变能
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       
       !计算ni点的微模量,注意PMB泊松比为固定值
       emod=mats(matsnode(ni))%emod
       !v=mats(matsnode(ni))%v   
       if (ndime==3) then !c=18*k/(pi*delta**4)
           !v=1/4
           !ci=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
           
           !简化后
           ci=12*emod/(pi*delta**4)
       elseif (ndime==2) then   !c=12*k/(pi*delta**3)
           if (planetype=='stress')then
               !v=1/3
               !ci=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
               !简化后
               ci=9*emod/(pi*delta**3)!Plane stress
            elseif(planetype=='strain')then
               !!v=1/4
               !ci=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
               !简化后
               ci=48/5*emod/(pi*delta**3)!Plane strain
            else
               write(*,*)planetype,'平面问题类型出错！！！'
               stop
            endif
       else    !c=2E/delta^2
           ci=2*emod/(pi*delta**2)
       end if 
       
       !读取ni点的微模量
       !ci=mods(matsnode(ni))%c

       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          
          !计算nj点的微模量，PMB模型中泊松比为固定值
          emod=mats(matsnode(nj))%emod
          !v=mats(matsnode(nj))%v          
          if (ndime==3) then !c=18*k/(pi*delta**4)
              !v=1/4
              !cj=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
              cj=12*emod/(pi*delta**4)
          elseif (ndime==2) then   !c=12*k/(pi*delta**3)
              if (planetype=='stress')then
                 !v=1/3
                 !cj=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
                 !简化后
                 cj=9*emod/(pi*delta**3)!Plane stress
               elseif(planetype=='strain')then
                 !v=1/4
                 !cj=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
                 !简化后
                 cj=48/5*emod/(pi*delta**3)!Plane strain
               else
                  write(*,*)planetype,'平面问题类型出错！！！'
                  stop
                endif
          else    !c=2E/delta^2
              cj=2*emod/(pi*delta**2)
          end if  
          
          !读取nj点的微模量
          !cj=mods(matsnode(nj))%c
          
          !当键两端材料不同时，微模量取串联值
          c=2/(1/ci+1/cj)
          !计算键长
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !读取初始键长
          idist= idistfam(pointfam(ni)+j-1)
          !计算键伸长量并存储
          e=nlength-idist
          !修正i领域内j点积分体积Vj         
          call sub_v_fac()
          
          Strain_Energy(ni)=Strain_Energy(ni)+1/4*c*e*e/idist*vol(nj)*v_fac
       end do
    end do
end subroutine sub_Strain_Energy_PMB