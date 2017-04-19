!***************************************************************
!                         data_input                           *         
!             功能：从数据文件读入数据                         *
!***************************************************************
    subroutine data_input()
    use scalars
    use serialArrays  
    use globalParameters
    implicit none
    character(len=255)::a
    integer::i,j,k,m,n,num,ilist,inode,maxincs
    integer,allocatable::matsblock(:)
    type boundary
         integer numlist
         integer list(12)
         integer  xyz(12)
         real*8  value(12)
    end type boundary
    type(boundary),allocatable::bc_d(:),bc_v(:),bc_f(:)
    type crack_type
         real*8 x
         real*8 y
         real*8 theta
         real*8 l
         real*8 d
    end type crack_type
    type(crack_type),allocatable::crack(:)
    
!***读取主控信息  
    read(5,*)a
    !读取物质点总数，结构总块数，合集物质点总数，集合总数,荷载步
    read(5,*)totnode,nblock,totlist,nlist,ncrack,nload
    allocate(coord(totnode,3),block_id(totnode),vol(totnode),damage(totnode),nodeyes(totnode))
    allocate(numfam(totnode),pointfam(totnode),nodefam(totnode*maxfam),efam(totnode*maxfam),edfam(totnode*maxfam),idistfam(totnode*maxfam))
    allocate(cnumfam(totnode),cpointfam(totnode),cnodefam(totnode*maxfam),cidistfam(totnode*maxfam))
    allocate(force(totnode,3),pforce(totnode,3),pforceold(totnode,3),bforce(totnode,3),cforce(totnode,3),gforce(totnode,3))
    allocate(disp(totnode,3),vel(totnode,3),velhalf(totnode,3),velhalfold(totnode,3),acc(totnode,3))
    allocate(wvolume(totnode),theta(totnode))
    allocate(typenode(nload,totnode))
    allocate(s_fac(totnode),Strain_Energy(totnode))
!***读取模型信息    
    read(5,*)a
    !读取物质点间距，近场物质点倍数，维数，材料模型
    read(5,*)dx,n_delta,ndime,model
    read(5,*)planetype  !读取平面模型类型
!***读取网格信息  
    allocate(numlist(nlist),pointlist(nlist),nodelist(totlist),faclist(totlist))
    read(4,*)a
    !读取坐标、分块及体积
    do i=1,totnode
        read(4,*) coord(i,1:3),block_id(i),vol(i)
    end do
    
    read(4,*)a
    !读取各集合物质点数目
    read(4,*)numlist(1:nlist)
    pointlist(1)=1
    !计算各集合起始点位置
    do i=2,nlist
        pointlist(i)=pointlist(i-1)+numlist(i-1)
    end do
    read(4,*)a
    !顺序读取集合所有物质点号
    do i=1,totlist
        read(4,*) nodelist(i),faclist(i)
    end do
!***读取裂纹信息
    if(ncrack>0)then
        open(3,file='crackpoint.dat',status='old')
        allocate(crack(ncrack))
        allocate(cracks(ncrack,5))
        read(3,*) a
        !输入裂纹xy平面中心点坐标、角度、长度、宽度
        do i=1,ncrack
            read(3,*) crack(i)
        end do
    end if
    do i=1,ncrack
        cracks(i,1)=crack(i)%x-crack(i)%l*dcos(crack(i)%theta*pi/180)/2
        cracks(i,2)=crack(i)%y-crack(i)%l*dsin(crack(i)%theta*pi/180)/2
        cracks(i,3)=crack(i)%x+crack(i)%l*dcos(crack(i)%theta*pi/180)/2
        cracks(i,4)=crack(i)%y+crack(i)%l*dsin(crack(i)%theta*pi/180)/2
        cracks(i,5)=crack(i)%d
    end do
!***读取分块信息
    allocate(numblock(nblock),pointblock(nblock),nodeblock(totnode))
    k=0
    do i=1,nblock
        numblock(i)=0
        do j=1,totnode
            if(block_id(j)==i)then
                k=k+1
                numblock(i)=numblock(i)+1
                nodeblock(k)=j
            end if
        end do
    end do
    pointblock(1)=1
    do i=2,nblock
        pointblock(i)=pointblock(i-1)+numblock(i-1)
    end do
    
!***读取材料信息    
    read(5,*)a
    read(5,*)nmats
    !读取各材料参数
    allocate(mats(nmats),mods(nmats))
    do i=1,nmats
        read(5,*)mats(i)
    end do
    !为材料块指定材料
    allocate(matsblock(nblock))
    read(5,*) matsblock(1:nblock)
    !为物质点附加材料信息
    allocate(matsnode(totnode))
    do i=1,nblock
        do j=1,numblock(i)
            matsnode(nodeblock(pointblock(i)+j-1))=matsblock(i)
        end do
    end do
   

    allocate(bc_d(nload),bc_v(nload),bc_f(nload))
    allocate(bcyes(nload,3),bcnode(nload,totlist*3),bcvalue(nload,totlist*3))
    allocate(numbc_d(nload,3),pointbc_d(nload,3),numbc_v(nload,3),pointbc_v(nload,3),numbc_f(nload,3),pointbc_f(nload,3))
    allocate(solve(nload),mout(nload))
    !初始化赋值
    bcyes=0
    bcnode=0
    bcvalue=0.0d0
    numbc_d=0
    pointbc_d=0
    numbc_v=0
    pointbc_v=0
    numbc_f=0
    pointbc_f=0
    typenode=0
!***读取各加载步荷载求解及输出信息   
    do iload=1,nload
       read(5,*)a
!***读取边界信息
       read(5,*)a
       !读取初始信息
       read(5,*)bc_d(iload)%numlist,bc_v(iload)%numlist,bc_f(iload)%numlist
       read(5,*)a  
       if (bc_d(iload)%numlist/=0)then
           read(5,*)bc_d(iload)%list(1:bc_d(iload)%numlist)
           read(5,*)bc_d(iload)%xyz(1:bc_d(iload)%numlist)
           read(5,*)bc_d(iload)%value(1:bc_d(iload)%numlist)
       end if
       read(5,*)a
       if (bc_v(iload)%numlist/=0)then
           read(5,*)bc_v(iload)%list(1:bc_v(iload)%numlist)
           read(5,*)bc_v(iload)%xyz(1:bc_v(iload)%numlist)
           read(5,*)bc_v(iload)%value(1:bc_v(iload)%numlist)
       end if
       read(5,*)a
       if (bc_f(iload)%numlist/=0)then
           read(5,*)bc_f(iload)%list(1:bc_f(iload)%numlist)
           read(5,*)bc_f(iload)%xyz(1:bc_f(iload)%numlist)
           read(5,*)bc_f(iload)%value(1:bc_f(iload)%numlist)
       end if
    !生成边界全局信息
        n=0
        !生成力等效边界条件信息
        if (bc_f(iload)%numlist/=0)then
            bcyes(iload,3)=1
            do i=1,3
                m=0
                do j=1,bc_f(iload)%numlist
                    if (bc_f(iload)%xyz(j)==i)then
                        ilist=bc_f(iload)%list(j)
                        num=numlist(ilist)
                        do k=1,num
                           m=m+1
                           n=n+1
                           inode=nodelist(pointlist(ilist)+k-1)
                           bcnode(iload,n)=inode
                           typenode(iload,inode)=3
                           bcvalue(iload,n)=bc_f(iload)%value(j)*faclist(pointlist(ilist)+k-1)
                        end do
                    end if
                end do
                numbc_f(iload,i)=m
                 if(m.ne.0) pointbc_f(iload,i)=n-m+1
            end do
        end if
        !生成位移边界条件信息
        if (bc_d(iload)%numlist/=0)then
            bcyes(iload,1)=1
            do i=1,3
                m=0
                do j=1,bc_d(iload)%numlist
                    if (bc_d(iload)%xyz(j)==i)then
                        ilist=bc_d(iload)%list(j)
                        num=numlist(ilist)
                        do k=1,num
                           m=m+1
                           n=n+1
                           inode=nodelist(pointlist(ilist)+k-1)
                           bcnode(iload,n)=inode
                           typenode(iload,inode)=1
                           bcvalue(iload,n)=bc_d(iload)%value(j)*faclist(pointlist(ilist)+k-1)
                        end do
                    end if
                end do
                numbc_d(iload,i)=m
                if(m.ne.0) pointbc_d(iload,i)=n-m+1
            end do
        end if
        !生成速度边界条件信息
        if (bc_v(iload)%numlist/=0)then
            bcyes(iload,2)=1
            do i=1,3
                m=0
                do j=1,bc_v(iload)%numlist
                    if (bc_v(iload)%xyz(j)==i)then
                        ilist=bc_v(iload)%list(j)
                        num=numlist(ilist)
                        do k=1,num
                           m=m+1
                           n=n+1
                           inode=nodelist(pointlist(ilist)+k-1)
                           bcnode(iload,n)=inode
                           typenode(iload,inode)=2
                           bcvalue(iload,n)=bc_v(iload)%value(j)*faclist(pointlist(ilist)+k-1)
                        end do
                    end if
                end do
                numbc_v(iload,i)=m
                 if(m.ne.0)  pointbc_v(iload,i)=n-m+1
            end do
        end if

    
!***读取求解信息
        read(5,*)a 
        read(5,*)a 
    !读取加载步求解信息
        read(5,*) solve(iload) 
        
!***读取输出信息
        read(5,*)a 
    !读取加载步求解信息
        read(5,*)a 
        !输出块材料块
        read(5,*)mout(iload)%numblocks
        read(5,*)a 
        !输出块材料块编号
        read(5,*)mout(iload)%blocks(1:mout(iload)%numblocks)
        read(5,*)a 
        !输出频率
        read(5,*)mout(iload)%frequency
    end do
!***读取加载级
    
    maxincs=0
    do iload=1,nload
        if(maxincs<solve(iload)%nincs)then
            maxincs=solve(iload)%nincs
        end if
    end do  
    allocate(loadg(nload,maxincs))
    loadg=1
    do iload=1,nload
        do iincs=1,solve(iload)%nincs
            read(2,*)loadg(nload,iincs)
        end do
    end do
    
!***标记裂纹位置物质，生成有宽度裂纹    
    call node_delete() 
    
!***确定输出信息
    allocate(totout(nload),outnode(nload,totnode))
    do iload=1,nload
        totout(iload)=0
        do i=1,mout(iload)%numblocks
           j=mout(iload)%blocks(i)
           do m=1,numblock(j) 
               ni=nodeblock(pointblock(j)-1+m)
               if(nodeyes(ni).eq.0)cycle
               totout(iload)=totout(iload)+1
               outnode(iload,totout(iload))=ni
           end do
        end do
    end do
!***结束数据读取及预处理    
end subroutine data_input

    
!***************************************************************
!                      node_delete                             *
!          功能：    删除预制裂纹处物质点                      *
!***************************************************************
subroutine node_delete()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    real(8)::x_coord,y_coord
    real::x1,y1,x2,y2,x3,y3,d,d0,DistanceLine
    nodeyes=1
    do ni=1,totnode
    !标记有裂纹宽度范围内的物质点
		if(ncrack>0) then
            do icrack=1,ncrack
		        x1=cracks(icrack,1)
		        y1=cracks(icrack,2)
		        x2=cracks(icrack,3)
		        y2=cracks(icrack,4)
		        x3=coord(ni,1)
		        y3=coord(ni,2)
		        d=DistanceLine(x1,y1,x2,y2,x3,y3)
		        d0=cracks(icrack,5)
		        if(d<d0) nodeyes(ni)=0
            end do
        end if 
    enddo
	return
end subroutine node_delete
!*****************************************************    
    
!***************************************************************
!                  Distance                                    *
!   功能：求以点(x1,y1)和(x2,y2)为端点的线段的距离             *
!***************************************************************
!求两点间距离
    function Distance(x1,y1,x2,y2)	result(d)
	   implicit none
	   real::x1,y1,x2,y2,d
	   d=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
    end  function distance  
!***************************************************************
!                    DistanceLine                              *
!  功能：求点(x3,y3)到以(x1,y1)和(x2,y2)为端点的线段的距离     *
!*************************************************************** 
	 function DistanceLine(x1,y1,x2,y2,x3,y3) result(d)
	   implicit none
	   real::x1,y1,x2,y2,x3,y3,Distance,ab,bb,d
	   ab=(x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)
	   bb=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
	   if(ab<=0)then
	     d=Distance(x1,y1,x3,y3)
	   else if(ab>=bb)then
	     d=Distance(x2,y2,x3,y3)
	   else
	     d=abs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)/Distance(x1,y1,x2,y2)
	   endif
	 end function DistanceLine    