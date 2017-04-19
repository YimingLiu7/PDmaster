!***************************************************************
!                    model_bond_PMB                            *
!          ���ܣ�                                              *
!***************************************************************    
subroutine model_bond_PMB
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 emod,v,c,ci,cj,delx(3),pf
    integer i,j
!***��ʼ��
    pforce=0.0d0
!***����ni��������ܶ�
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       
       !����ni���΢ģ��,ע��PMB���ɱ�Ϊ�̶�ֵ
       emod=mats(matsnode(ni))%emod
       !v=mats(matsnode(ni))%v   
       if (ndime==3) then !c=18*k/(pi*delta**4)
           !v=1/4
           !ci=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
           
           !�򻯺�
           ci=12*emod/(pi*delta**4)
       elseif (ndime==2) then   !c=12*k/(pi*delta**3)
           if (planetype=='stress')then
               !v=1/3
               !ci=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
               !�򻯺�
               ci=9*emod/(pi*delta**3)!Plane stress
            elseif(planetype=='strain')then
               !!v=1/4
               !ci=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
               !�򻯺�
               ci=48/5*emod/(pi*delta**3)!Plane strain
            else
               write(*,*)planetype,'ƽ���������ͳ�������'
               stop
            endif
       else    !c=2E/delta^2
           ci=2*emod/(pi*delta**2)
       end if 
       
       !��ȡni���΢ģ��
       !ci=mods(matsnode(ni))%c

       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          
          !����nj���΢ģ����PMBģ���в��ɱ�Ϊ�̶�ֵ
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
                 !�򻯺�
                 cj=9*emod/(pi*delta**3)!Plane stress
               elseif(planetype=='strain')then
                 !v=1/4
                 !cj=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
                 !�򻯺�
                 cj=48/5*emod/(pi*delta**3)!Plane strain
               else
                  write(*,*)planetype,'ƽ���������ͳ�������'
                  stop
                endif
          else    !c=2E/delta^2
              cj=2*emod/(pi*delta**2)
          end if  
          
          !��ȡnj���΢ģ��
          !cj=mods(matsnode(nj))%c
          
          !�������˲��ϲ�ͬʱ��΢ģ��ȡ����ֵ
          c=2/(1/ci+1/cj)
          !�������
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !��ȡ��ʼ����
          idist= idistfam(pointfam(ni)+j-1)
          !������쳤�����洢
          e=nlength-idist
          efam(pointfam(ni)+j-1)=e
          !����i������j��������Vj         
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
!***��ʼ��
    Strain_Energy=0.0d0
!***����ni���Ӧ����
    do ni=1,totnode
       if (nodeyes(ni)==0)cycle
       
       !����ni���΢ģ��,ע��PMB���ɱ�Ϊ�̶�ֵ
       emod=mats(matsnode(ni))%emod
       !v=mats(matsnode(ni))%v   
       if (ndime==3) then !c=18*k/(pi*delta**4)
           !v=1/4
           !ci=18*emod/(pi*delta**4)/(3*(1-2*v))!k=E/(3*(1-2v))
           
           !�򻯺�
           ci=12*emod/(pi*delta**4)
       elseif (ndime==2) then   !c=12*k/(pi*delta**3)
           if (planetype=='stress')then
               !v=1/3
               !ci=12*emod/(pi*delta**3)/(2*(1-v))                     !k=E/(2*(1-v))!Plane stress
               !�򻯺�
               ci=9*emod/(pi*delta**3)!Plane stress
            elseif(planetype=='strain')then
               !!v=1/4
               !ci=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
               !�򻯺�
               ci=48/5*emod/(pi*delta**3)!Plane strain
            else
               write(*,*)planetype,'ƽ���������ͳ�������'
               stop
            endif
       else    !c=2E/delta^2
           ci=2*emod/(pi*delta**2)
       end if 
       
       !��ȡni���΢ģ��
       !ci=mods(matsnode(ni))%c

       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          
          !����nj���΢ģ����PMBģ���в��ɱ�Ϊ�̶�ֵ
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
                 !�򻯺�
                 cj=9*emod/(pi*delta**3)!Plane stress
               elseif(planetype=='strain')then
                 !v=1/4
                 !cj=12*emod/(pi*delta**3)/(2*(1-v-2*v**2))              !k=E/(2*(1-v-2v**2))!Plane strain
                 !�򻯺�
                 cj=48/5*emod/(pi*delta**3)!Plane strain
               else
                  write(*,*)planetype,'ƽ���������ͳ�������'
                  stop
                endif
          else    !c=2E/delta^2
              cj=2*emod/(pi*delta**2)
          end if  
          
          !��ȡnj���΢ģ��
          !cj=mods(matsnode(nj))%c
          
          !�������˲��ϲ�ͬʱ��΢ģ��ȡ����ֵ
          c=2/(1/ci+1/cj)
          !�������
          nlength=0.0d0
          do idime=1,ndime
             delx(idime)=coord(nj,idime)+disp(nj,idime)-coord(ni,idime)-disp(ni,idime)
             nlength=nlength+delx(idime)*delx(idime)
          end do
          nlength = dsqrt(nlength)
          !��ȡ��ʼ����
          idist= idistfam(pointfam(ni)+j-1)
          !������쳤�����洢
          e=nlength-idist
          !����i������j��������Vj         
          call sub_v_fac()
          
          Strain_Energy(ni)=Strain_Energy(ni)+1/4*c*e*e/idist*vol(nj)*v_fac
       end do
    end do
end subroutine sub_Strain_Energy_PMB