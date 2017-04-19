!******************************************************************
!     module scalars
!******************************************************************
    module scalars
    integer totnode,nblock,totlist,nlist,nload,ncrack,icrack,ni,nj
    integer ndime,idime,nmats
    integer iload,istep,iincs,ii
    character(len=3) model
    character(len=6) planetype
    !planetype:平面问题类型
    real*8 dx,n_delta,delta,v_fac,idist,nlength
!***check参数
    integer check_c
    !check_c:收敛结果：1收敛；0未收敛
!***state参数
    real*8 omega,e,ed
    !omega:influence function
    real*8 maxdamage
    end module scalars
    
!******************************************************************
!     module serialArrays
!******************************************************************
    module serialArrays


!***物质点相关信息   
    real(8),allocatable::coord(:,:),vol(:),damage(:)
    real(8),allocatable::disp(:,:),vel(:,:),velhalf(:,:),velhalfold(:,:),acc(:,:)
    real(8),allocatable::force(:,:),pforce(:,:),pforceold(:,:),bforce(:,:),cforce(:,:),gforce(:,:)
    integer,allocatable::block_id(:),matsnode(:),nodeyes(:)
    integer,allocatable::numfam(:),pointfam(:),nodefam(:),cnumfam(:),cpointfam(:),cnodefam(:)
    real(8),allocatable::efam(:),edfam(:),idistfam(:),cidistfam(:)
    real(8),allocatable::wvolume(:),theta(:)
    real(8),allocatable::s_fac(:),Strain_Energy(:)
    !wvolume(:):weighted volume
    !theta(:):dilatation
!***物质点集合相关信息
    integer,allocatable::numlist(:),pointlist(:),nodelist(:)
    real(8),allocatable::faclist(:)
!***材料相关变量 
    type materials
!        character(len=10)::m_type
         real*8 emod    !Elastic modulus
         real*8 v       !Poisson's ratio
         real*8 kmod    !Bulk Modulus
         real*8 gmod    !Shear Modulus  
         real*8 dens    !Density
         real*8 st      !Critical Stretch拉
         real*8 sc      !Critical Stretch压
    end type materials
    type(materials),allocatable::mats(:)
    type materials_mod
         real*8 c    
         real*8 d       
         real*8 k    
         real*8 alpha            
         real*8 gama   
    end type materials_mod
    type(materials_mod),allocatable::mods(:)
!***裂纹相关信息
    real(8),allocatable::cracks(:,:)
!***分块信息 
    integer,allocatable::numblock(:),nodeblock(:),pointblock(:)
    
!***边界相关变量
    integer,allocatable::bcyes(:,:),bcnode(:,:)
    real,allocatable::bcvalue(:,:)
    integer,allocatable::numbc_d(:,:),pointbc_d(:,:),numbc_v(:,:),pointbc_v(:,:),numbc_f(:,:),pointbc_f(:,:)
    integer,allocatable::typenode(:,:)
    !typenode：物质点类型，1位移边界点，2速度边界点，3力边界点，0自由点，-1，边界领域点
!***求解相关信息
    type solve_type
         character(len=3)::stype
         real*8  dt       !迭代时间步长
         real*8  cv       !阻尼系数
         real*8  toler    !收容差
         integer maxstep  !最大迭代步
         integer nincs    !增量步
         integer ncrit    !破坏准则
         integer loadtype!加载方式控制是否有加载域
    end type solve_type
    type(solve_type),allocatable::solve(:)
    
    real(8),allocatable::loadg(:,:)!荷载系数列阵
!***输出相关信息
    type output 
         integer::numblocks  !输出块材料块数
         integer,dimension(10)::blocks !材料块编号
         integer frequency   !输出频率      
    end type output
    type(output),allocatable::mout(:)
    integer,allocatable::totout(:),outnode(:,:)
    end module serialArrays
    
!******************************************************************
!     module serialArrays
!******************************************************************   
    module globalParameters
    integer,parameter::maxfam=100
    !maxfam: Maximum number of material points inside a horizon of a material point
    real(8),parameter::pi=dacos(-1.0d0),g=9.81
    !pi: Pi
    !g: Acceleration of gravity(m/s^2)

    end module globalParameters	
    