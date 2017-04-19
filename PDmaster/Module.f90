!******************************************************************
!     module scalars
!******************************************************************
    module scalars
    integer totnode,nblock,totlist,nlist,nload,ncrack,icrack,ni,nj
    integer ndime,idime,nmats
    integer iload,istep,iincs,ii
    character(len=3) model
    character(len=6) planetype
    !planetype:ƽ����������
    real*8 dx,n_delta,delta,v_fac,idist,nlength
!***check����
    integer check_c
    !check_c:���������1������0δ����
!***state����
    real*8 omega,e,ed
    !omega:influence function
    real*8 maxdamage
    end module scalars
    
!******************************************************************
!     module serialArrays
!******************************************************************
    module serialArrays


!***���ʵ������Ϣ   
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
!***���ʵ㼯�������Ϣ
    integer,allocatable::numlist(:),pointlist(:),nodelist(:)
    real(8),allocatable::faclist(:)
!***������ر��� 
    type materials
!        character(len=10)::m_type
         real*8 emod    !Elastic modulus
         real*8 v       !Poisson's ratio
         real*8 kmod    !Bulk Modulus
         real*8 gmod    !Shear Modulus  
         real*8 dens    !Density
         real*8 st      !Critical Stretch��
         real*8 sc      !Critical Stretchѹ
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
!***���������Ϣ
    real(8),allocatable::cracks(:,:)
!***�ֿ���Ϣ 
    integer,allocatable::numblock(:),nodeblock(:),pointblock(:)
    
!***�߽���ر���
    integer,allocatable::bcyes(:,:),bcnode(:,:)
    real,allocatable::bcvalue(:,:)
    integer,allocatable::numbc_d(:,:),pointbc_d(:,:),numbc_v(:,:),pointbc_v(:,:),numbc_f(:,:),pointbc_f(:,:)
    integer,allocatable::typenode(:,:)
    !typenode�����ʵ����ͣ�1λ�Ʊ߽�㣬2�ٶȱ߽�㣬3���߽�㣬0���ɵ㣬-1���߽������
!***��������Ϣ
    type solve_type
         character(len=3)::stype
         real*8  dt       !����ʱ�䲽��
         real*8  cv       !����ϵ��
         real*8  toler    !���ݲ�
         integer maxstep  !��������
         integer nincs    !������
         integer ncrit    !�ƻ�׼��
         integer loadtype!���ط�ʽ�����Ƿ��м�����
    end type solve_type
    type(solve_type),allocatable::solve(:)
    
    real(8),allocatable::loadg(:,:)!����ϵ������
!***��������Ϣ
    type output 
         integer::numblocks  !�������Ͽ���
         integer,dimension(10)::blocks !���Ͽ���
         integer frequency   !���Ƶ��      
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
    