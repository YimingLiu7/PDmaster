!***************************************************************
!                    Time_Integration                          *
!          ���ܣ� ʩ�ӵ�Ч����ܶȺ���                         *
!***************************************************************    
subroutine Time_Integration()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 dens,cv,dt
    cv=solve(iload)%cv
    dt=solve(iload)%dt
!***����λ��
    do idime=1,ndime
       do ni=1,totnode
          disp(ni,idime)=disp(ni,idime)+velhalfold(ni,idime)*dt
       end do
    end do
!***ʩ��λ�ƺ���
    call Applied_Displacement      
!***��������
    call internal_force
!***ʩ����������
    call applied_force()
!***����Ӵ���
!    call contect_force
!***�������force   
    do idime=1,ndime
       do ni=1,totnode
          force(ni,idime)=pforce(ni,idime)+bforce(ni,idime)!+cforce(ni,idime)+gforce(ni,idime)
       end do
    end do
!***�����ٶ�velhalf
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
!***ʩ���ٶȺ���velhalf
    call Applied_Velocity()      

end subroutine Time_Integration
    
    
    
    
!***************************************************************
!                  check_convergence                           *
!          ���ܣ� ʩ�ӵ�Ч����ܶȺ���                         *
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
            !λ�Ƽ��ص�1���ٶȼ��ص�2������
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
!          ���ܣ� ʩ�ӵ�Ч����ܶȺ���                         *
!***************************************************************    
subroutine check_damage()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real*8 s,st,sc,fac
    integer i,j,numbond
!***�жϼ�����
    do ni=1,totnode 
       !�ж����ʵ�ni�Ƿ��޳�
       if (nodeyes(ni)==0)cycle
       !�ж����ʵ�ni�Ƿ�Ϊ�߽���߽�����
       if(solve(iload)%loadtype==1)then
           if (typenode(iload,ni).ne.0)cycle 
       endif
       do j=1,numfam(ni)
          nj=nodefam(pointfam(ni)+j-1)
          if (nj==0)cycle
          !�ж��Ƿ��м�����
          if(solve(iload)%loadtype==1)then
             !���ص�123�Լ�������-1��ǵĵ㲻�ƻ�
             if (typenode(iload,nj).ne.0)cycle
          endif
          s=efam(pointfam(ni)+j-1)/idistfam(pointfam(ni)+j-1)
          
          !���ò��Ͻ���ϵ��
          fac=1.0d0
          if(matsnode(ni).ne.matsnode(nj)) fac=1.2d0
          
          st=min((mats(matsnode(ni))%st),(mats(matsnode(nj))%st))*fac
          sc=min((mats(matsnode(ni))%sc),(mats(matsnode(nj))%sc))*fac
          if (s>st.or.s<-sc)then
             nodefam(pointfam(ni)+j-1)=0
          end if
       end do
    end do
    
!***�������ʵ�����    
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