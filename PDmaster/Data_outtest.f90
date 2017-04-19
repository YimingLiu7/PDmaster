!***************************************************************
!               data_outtest                                   *
!          功能： 控制变量输出                                 *
!***************************************************************
subroutine data_outtest()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer::i,j,blockid
    !real*8 sum
    !write(8,*) istep,disp(59141,1),disp(59141,2)
    !sum=0.0
    !blockid=1 
    !do i=1,numblock(blockid)
    !   ni=nodeblock(pointblock(blockid)+i-1)
    !   sum=sum+pforce(ni,1)*vol(ni)
    !end do
    !write(8,*)istep,sum
    
end subroutine data_outtest