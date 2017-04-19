!***************************************************************
!                  contect_force                               *
!          ���ܣ� ʩ�ӵ�Ч����ܶȺ���                         *
!***************************************************************    
subroutine contect_force()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,c,delx(3),cf,dij
    integer i,j
!***��ʼ��
    cforce=0.0d0
!***����ni��������ܶ�
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       !����ni���΢ģ��,ע��PMB���ɱ�Ϊ�̶�ֵ
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
               write(*,*)planetype,'ƽ���������ͳ�������'
               stop
            endif
       else    !c=2E/delta^2
           c=2*emod/(pi*delta**2)
       end if 
       
       !��ȡni���΢ģ��
       !c=mods(matsnode(ni))%c

       do j=1,numfam(ni)
          nj=cnodefam(cpointfam(ni)+j-1)
          if (nj==0)cycle
          !�������
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !��ȡ��ʼ����
          idist= cidistfam(cpointfam(ni)+j-1)
          !������쳤�����洢
          dij=min((0.9*idist),(1.35*dx))
          !����i������j��������Vj  
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