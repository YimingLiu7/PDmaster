!***************************************************************
!                    model_state_LPS                           *
!          ���ܣ�                                              *
!***************************************************************         
subroutine model_state_LPS()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 kmod,gmod,v,k,alpha,gama,t,delx(3)
    integer i,j
!!***����
!    do idime=1,ndime
!        do ni=1,totnode
!           disp(ni,idime)=0.001*coord(ni,idime)
!        end do
!    end do    
   
!***�������ʵ�dilatation
    theta=0.0d0
    do ni=1,totnode 
       if (nodeyes(ni)==0)cycle
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          !����ni��nj����
          if (nj==0)cycle 
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !��ȡni��nj��ʼ����
          idist= idistfam(pointfam(ni)+j-1)
          !���������
          e=nlength-idist
          efam(pointfam(ni)+j-1)=e
          !�����������ϵ��,����������
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
              write(*,*)planetype,'ƽ���������ͳ�������'
              stop
          endif
             
       else
          !1D
       end if
    end do
    
    
!***��ʼ��
    pforce=0.0d0
!***�������ʵ������ܶ�
    do ni=1,totnode 
       if (nodeyes(ni)==0)cycle
       
       !����ģ�Ͳ���
       kmod=mats(matsnode(ni))%kmod
       gmod=mats(matsnode(ni))%gmod
       v=mats(matsnode(ni))%v 
       !����΢ģ��
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
              write(*,*)planetype,'ƽ���������ͳ�������'
              stop
          endif
       else 
            !1D
       endif 
       
       !ֱ�Ӷ�ȡģ�Ͳ���
       !k=mods(matsnode(ni))%k
       !alpha=mods(matsnode(ni))%alpha
       !gama=mods(matsnode(ni))%gama
       
       !����ni�����������ʵ�
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle      
          !������ʱ�������ʵ���
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength= nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          
          idist= idistfam(pointfam(ni)+j-1)
          !�����������
          e=nlength-idist 
          !�����������ƫ������
          !ed=e-theta(ni)*idist/3
          call sub_v_fac() 
          call sub_omega()
          !����������״̬
          !��ǰ
          !t=(gama*k/wvolume(ni)-alpha/3*(1-gama/3))*theta(ni)*omega*idist+alpha*omega*ed
          !�򻯺�
          t=(gama*k/wvolume(ni)-alpha/3*(2-gama/3))*theta(ni)*omega*idist+alpha*omega*e
          
          !if(istep==1) write(8,*)ni,nj,t,e,k,gama,wvolume(ni),theta(ni)
          
          do idime=1,ndime
             pforce(ni,idime)=pforce(ni,idime)+t*delx(idime)/nlength*vol(nj)*v_fac
             pforce(nj,idime)=pforce(nj,idime)-t*delx(idime)/nlength*vol(ni)*v_fac
          end do
       end do
    end do
end subroutine model_state_LPS 