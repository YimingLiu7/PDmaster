program master
    use scalars
    use serialArrays  
    use globalParameters
    use CMD_Progress
    implicit none
    character(len=20)::str1
    type( CLS_CMD_Progress ) ::Progress 
!    write(Progress_name, '(1X,A,I2,A,I2,A)' ) "The Progress of ",iincs,'/',nincs,': '
!    Progress % Prefix = AdjustL(Trim( Progress_name) )!设置进度条名称
!***开辟数据通道
    open(2,file='load_grading.dat',status='old')
    open(4,file='mesh.dat',status='old')
    open(5,file='input.dat',status='old')
    open(6,file='output.dat')
    open(7,file='out_c.dat')
    open(8,file='outtest.dat')
!***调用数据读入子程序
    call data_input() 
    write(*,*)'end data_input'
!***识别近场物质点
    call element_form
    write(*,*)'end element_form'
!***变量初始化赋值
    call initialization()
    write(*,*)'end initialization'
!***求解计算
    ii=-1
    do iload=1,nload
        call Progress % Set( N = solve(iload)%maxstep, L = 25 )!// 进度条显示长度25
        write(str1,'(1X,I3)') solve(iload)%nincs
        do iincs=1,solve(iload)%nincs
            
            write(Progress_name, '(1X,A,I2,A,I3,A,A,A)' ) "The Progress of ",iload,':',iincs,'/',AdjustL(Trim(str1)),': '
            Progress % Prefix = AdjustL(Trim(Progress_name))!设置进度条名称
            
            check_c=0
            do istep=1,solve(iload)%maxstep
                
                !识别接触力作用范围
                !if(mod(istep,100)==0) call celement_form()   
                
                call Time_Integration()
                !检测损伤
                if(solve(iload)%ncrit==1)then
                   call check_damage()
                end if
                !检测收敛性
                call check_convergence()
                call data_output()
                call data_outtest()  !//输出测试信息
                
                call Progress % Put( istep , CMD_PROGRESS_ABSOLUTE ) !// 观察迭代过程
                if(check_c==1) exit
            end do
            !收敛后检测损伤
            if (solve(iload)%ncrit==2)then
                call check_damage()
            end if
            
            !如果出现损伤且达到最大迭代步则退出计算
            if (maxdamage>0.3.and.istep>=solve(iload)%maxstep) exit
            
        end do
        
    end do
    write(*,*)
    write(*,*)'end PDmaster'
end program master  