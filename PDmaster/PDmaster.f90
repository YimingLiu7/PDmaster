program master
    use scalars
    use serialArrays  
    use globalParameters
    use CMD_Progress
    implicit none
    character(len=20)::str1
    type( CLS_CMD_Progress ) ::Progress 
!    write(Progress_name, '(1X,A,I2,A,I2,A)' ) "The Progress of ",iincs,'/',nincs,': '
!    Progress % Prefix = AdjustL(Trim( Progress_name) )!���ý���������
!***��������ͨ��
    open(2,file='load_grading.dat',status='old')
    open(4,file='mesh.dat',status='old')
    open(5,file='input.dat',status='old')
    open(6,file='output.dat')
    open(7,file='out_c.dat')
    open(8,file='outtest.dat')
!***�������ݶ����ӳ���
    call data_input() 
    write(*,*)'end data_input'
!***ʶ��������ʵ�
    call element_form
    write(*,*)'end element_form'
!***������ʼ����ֵ
    call initialization()
    write(*,*)'end initialization'
!***������
    ii=-1
    do iload=1,nload
        call Progress % Set( N = solve(iload)%maxstep, L = 25 )!// ��������ʾ����25
        write(str1,'(1X,I3)') solve(iload)%nincs
        do iincs=1,solve(iload)%nincs
            
            write(Progress_name, '(1X,A,I2,A,I3,A,A,A)' ) "The Progress of ",iload,':',iincs,'/',AdjustL(Trim(str1)),': '
            Progress % Prefix = AdjustL(Trim(Progress_name))!���ý���������
            
            check_c=0
            do istep=1,solve(iload)%maxstep
                
                !ʶ��Ӵ������÷�Χ
                !if(mod(istep,100)==0) call celement_form()   
                
                call Time_Integration()
                !�������
                if(solve(iload)%ncrit==1)then
                   call check_damage()
                end if
                !���������
                call check_convergence()
                call data_output()
                call data_outtest()  !//���������Ϣ
                
                call Progress % Put( istep , CMD_PROGRESS_ABSOLUTE ) !// �۲��������
                if(check_c==1) exit
            end do
            !������������
            if (solve(iload)%ncrit==2)then
                call check_damage()
            end if
            
            !������������Ҵﵽ�����������˳�����
            if (maxdamage>0.3.and.istep>=solve(iload)%maxstep) exit
            
        end do
        
    end do
    write(*,*)
    write(*,*)'end PDmaster'
end program master  