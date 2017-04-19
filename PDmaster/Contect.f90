!***************************************************************
!                  contect_force                               *
!          功能： 施加等效体积密度荷载                         *
!***************************************************************    
subroutine contect_force()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,c,delx(3),cf,dij
    integer i,j
!***初始化
    cforce=0.0d0
!***计算ni点的内力密度
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       !计算ni点的微模量,注意PMB泊松比为固定值
       emod=mats(matsnode(ni))%emod
       !v=mats(matsnode(ni))%v   
       if (ndime==3) then !c=18*k/(pi*delta**4)
           c=9*emod/(pi*delta**4)
       elseif (ndime==2) then   !c=12*k/(pi*delta**3)
           if (planetype=='stress')then
               c=9*emod/(pi*delta**3)!Plane stress
            elseif(planetype=='strain')then
               c=48/5*emod/(pi*delta**3)!Plane strain
            else
               write(*,*)planetype,'平面问题类型出错！！！'
               stop
            endif
       else    !c=2E/delta^2
           c=2*emod/(pi*delta**2)
       end if 
       
       !读取ni点的微模量
       !c=mods(matsnode(ni))%c

       do j=1,numfam(ni)
          nj=cnodefam(cpointfam(ni)+j-1)
          if (nj==0)cycle
          !计算键长
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !读取初始键长
          idist= cidistfam(cpointfam(ni)+j-1)
          !计算键伸长量并存储
          dij=min((0.9*idist),(1.35*dx))
          !修正i领域内j点积分体积Vj  
          if(nlength<dij)then
              call sub_v_fac()
              cf=c*(dij-nlength)/delta*vol(nj)*v_fac
              do idime=1,ndime
                 cforce(ni,idime)=cforce(ni,idime)-cf*delx(idime)/nlength
                 cforce(nj,idime)=cforce(nj,idime)+cf*delx(idime)/nlength
              end do
          end if
          
       end do
    end do
end subroutine contect_force