!***************************************************************
!                    model_bond_CDG                            *
!          功能：                                              *
!***************************************************************    
subroutine model_bond_CDG
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,c,ci,cj,d,di,dj,delx(3),l,m,n,k(3,3),dis(3)
    integer i,j
!!!!***测试
!    do idime=1,ndime
!        do ni=1,totnode
!           disp(ni,idime)=0.001*coord(ni,idime)
!        end do
!    end do 
    
!***初始化
    pforce=0.0d0
!***计算ni点的内力密度
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       
       !计算ni点的微模量
       emod=mats(matsnode(ni))%emod
       v=mats(matsnode(ni))%v   
       if (ndime==3) then 
           !3D
           ci=6*emod/(pi*delta**4*(1-2*v))
           di=3*emod*(1-4*v)/(pi*delta*delta*(1-2*v)*(1+v))
       elseif (ndime==2) then
          if (planetype=='stress')then
             !Plane stress
             ci=6*emod/(pi*delta**3*(1-v))
             di=emod*(1-3*v)/(6*pi*delta*(1-v*v))
          elseif(planetype=='strain')then
             !Plane strain
             ci=6*emod/(pi*delta**3*(1-2*v)*(1+v))
             di=emod*(1-4*v)/(6*pi*delta*(1-2*v)*(1+v))
          else
              write(*,*)planetype,'平面问题类型出错！！！'
              stop
          endif
       else   
           ci=2*emod/(pi*delta**2)
           di=0.0d0
       end if 
       
       !直接读取ni点微模量
       !ci=mods(matsnode(ni))%c
       !di=mods(matsnode(ni))%d
       
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          
          !计算nj点的微模量
           emod=mats(matsnode(nj))%emod
           v=mats(matsnode(nj))%v   
           if (ndime==3) then 
              !3D
              cj=6*emod/(pi*delta**4*(1-2*v))
              dj=3*emod*(1-4*v)/(pi*delta*delta*(1-2*v)*(1+v))
           elseif (ndime==2) then
              if (planetype=='stress')then
                 !Plane stress
                 cj=6*emod/(pi*delta**3*(1-v))
                 dj=emod*(1-3*v)/(6*pi*delta*(1-v*v))
              elseif(planetype=='strain')then
                 !Plane strain
                 cj=6*emod/(pi*delta**3*(1-2*v)*(1+v))
                 dj=emod*(1-4*v)/(6*pi*delta*(1-2*v)*(1+v))
              else
                  write(*,*)planetype,'平面问题类型出错！！！'
                  stop
              endif
           else   
              cj=2*emod/(pi*delta**2)
              dj=0.0d0
           end if 
          !直接读取nj点微模量
          !cj=mods(matsnode(nj))%c
          !dj=mods(matsnode(nj))%d
          
          !当键两端材料不同时，微模量取串联值
          c=2/(1/ci+1/cj)
          d=2/(1/di+1/dj)
          !计算键长
          nlength=0.0d0
          delx=0.0d0
          dis=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             dis(idime)=disp(nj,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !读取初始键长
          idist= idistfam(pointfam(ni)+j-1)
          !计算键伸长量并存储
          e=nlength-idist
          efam(pointfam(ni)+j-1)=e
          !与坐标轴夹角余弦
          l=delx(1)/nlength
          m=delx(2)/nlength
          n=delx(3)/nlength
          if(ndime<3)n=1
          if(ndime<2)m=1
          k(1,1)=l*l*c/idist+12*l*l*m*m*d/((l*l+n*n)*idist**3)+12*n*n*d/((l*l+n*n)*idist**3)
          k(1,2)=l*m*c/idist-12*l*m*d/idist**3
          k(1,3)=l*n*c/idist+12*l*m*m*n*d/((l*l+n*n)*idist**3)+12*l*n*d/((l*l+n*n)*idist**3)
          k(2,1)=k(1,2)
          k(2,2)=m*m*c/idist+12*(l*l+n*n)*d/idist**3
          k(2,3)=m*n*c/idist-12*m*n*d/idist**3
          k(3,1)=k(1,3)
          k(3,2)=k(2,3)
          k(3,3)=n*n*c/idist+12*m*m*n*n*d/((l*l+n*n)*idist**3)+12*l*l*d/((l*l+n*n)*idist**3)
          
          !修正i领域内j点积分体积Vj
          call sub_v_fac()
          !计算合力
          do idime=1,ndime
             do i=1,ndime
                pforce(ni,idime)=pforce(ni,idime)+k(idime,i)*dis(i)*vol(nj)*v_fac
             end do
          end do
          
          !if(istep==50) write(8,*)pforce(ni,1),pforce(ni,2)
          
       end do
    end do
end subroutine model_bond_CDG