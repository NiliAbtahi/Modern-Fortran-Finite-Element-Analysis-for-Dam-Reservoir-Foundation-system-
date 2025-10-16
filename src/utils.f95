!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!


!!
!!
!!  precision kind
!!
!!  real64 = 64bit real or 8byte or double precision
!!  int8   = 8bit  integer takes value from -128 to 127     : all = 256   = 2^8
!!  int16  = 16bit integer takes value from -32768 to 32767 : all = 65536 = 2^16
!!  int32  = 32bit integer takes value from -2147483649 to 2147483648   : all = 4294967296 = 2^32
!!
!!  default integer in fortran is 4bytes or 32 bits
!!
!!  element_type = 1  : Q8
!!  element_type = 2  : Brick
!!

module  mod_utils
   
   
   use iso_fortran_env, only: real32, real64, real128 , output_unit
   
   implicit none
   
   integer, parameter , public :: r_kind     =  real64
   integer, parameter , public :: qp_kind    =  real128
   integer, parameter , public :: stdout     =  output_unit
   
   
   
   

!!====================================================================!!
!! define different kind of elements here. To every element an array
!! of size 4 is assigned whose elements are specified as follows:
!!
!! 1. num_main_nodes
!! 2. num_node_in_element
!! 3. num_element_around_node
!! 4. num_face
!!====================================================================!!


  integer , parameter , dimension( 2 , 4) , public :: elements_base_info =   reshape(   &
  [ 4 , 8  , 4  , 4  ,   &                                                 ! row 1 is for Q8 element or elem_id = 1
 &  8 , 20 , 8  , 6  ] , shape( elements_base_info ) , order = [ 2 , 1 ] ) ! row 2 is for brick element or elem_id = 2 
           
           
              
  
!!=========================================================================================================!!
 
   integer , dimension(1:2,1:8 )  , target  , public :: Q8_node_parent = reshape(  &
   & (/ 1 , 1 , -1 , -1 , 1 , 0 , -1 , 0 , &
   &   -1 , 1 ,  1 , -1 , 0 , 1 ,  0 ,-1 /) , (/2, 8/) , order = (/ 2 , 1 /) )
   
   
   
   integer , dimension(1:3,1:20)  , target  , public :: brick_node_parent = reshape( &
   & (/ 1 , 1 , -1 , -1 , 1 , 1 , -1 , -1 , 1 , 0 , -1 , 0 , 1 , 1 , -1 , -1 , 1 , 0 , -1 , 0 , &
   &   -1 , 1 ,  1 , -1 ,-1 , 1 ,  1 , -1 , 0 , 1 ,  0 ,-1 ,-1 , 1 ,  1 , -1 , 0 , 1 ,  0 ,-1 , &
   &   -1 ,-1 , -1 , -1 , 1 , 1 ,  1 ,  1 ,-1 ,-1 , -1 ,-1 , 0 , 0 ,  0 ,  0 , 1 , 1 ,  1 , 1   &
   & /) , (/3,20 /) , order = (/ 2 , 1 /) )
   
   
   

         
!!=====================================================================!!
!!   class for physical properties of material   
             
    type , public :: t_material     
           real(kind = r_kind) , public :: elastic_module               !E
           real(kind = r_kind) , public :: poisson_ratio                !POI
           real(kind = r_kind) , public :: special_weight               !GAMMA 
           real(kind = r_kind) , public :: n_factor                     !N = definition?
           real(kind = r_kind) , public :: sigma0_y                     !SIGY0 :plastic threshold: non-linear analysis
           real(kind = r_kind) , public :: h_prime                      !H' :non-linear analysis
    end type t_material
        

  
!!=====================================================================!!
!! class for boundary loads or condition, i.e. forces acting on faces: 
             
    type , public :: t_boundary_load   
           integer                  , public :: k_type                  ! KTYPE 
           real(    kind = r_kind ) , public :: strength                ! PR
           integer                  , public :: face_id                 ! NFACE
           integer                  , public :: direction               ! LDIR 
           real(    kind = r_kind ) , public :: reference_height        ! XIREFF
    end type t_boundary_load

      
!!=====================================================================!!

    type, public :: t_model_description
        integer, public :: accustic_bulk_model_id = 0 
        integer, public :: elastic_bulk_model_id = 0 
        integer, public :: reservior_inf_boundary_model_id = 0          ! 0 for Somerfled BC, 1 for HW model, 2 for GN model
        integer, public :: foundation_inf_boundary_model_id = 0         ! 0 for Lysmer BC,  1 for HW model
        integer, public :: order_of_reservior_higher_order_bc(2)  = (/ 0 , 0 /)  ! num_function_for_incident_wave , evanescent_wav
        integer, public :: order_of_foundation_higher_order_bc(2) = (/ 0 , 0 /) ! num_function_for_incident_wave , evanescent_wav
        
        real( kind = r_kind ), allocatable , public :: reservior_incident_coefficient(:)
        real( kind = r_kind ), allocatable , public :: reservior_evanescent_coefficient(:)
        real( kind = r_kind ), allocatable , public :: foundation_incident_coefficient(:)
        real( kind = r_kind ), allocatable , public :: foundation_evanescent_coefficient(:)
        
        contains 
           
               procedure , public  , pass :: initialize_model_description
               procedure , public  , pass :: free_model_description
           
    end type t_model_description

!!=====================================================================!!
    
    type, public :: t_enforced_disp
        integer, public :: num_nodes = 0                               !! NDIS in input file
        integer, allocatable, public:: node_ids(:)                     !! node ids
        real( kind = r_kind ), allocatable , public :: disp_vector(:,:)
        contains 
           
               procedure , public  , pass :: initialize_enforced_disp
               procedure , public  , pass :: free_enforced_disp
           
    end type t_enforced_disp
    
    
!!=====================================================================! 
!!  Virtual spring for fixing foundation boundary nodes in static step of time domain analysis

    type, public :: t_boundary_spring
        integer, public :: num_fixed_nodes = 0                               !! NDIS in input file
        integer, allocatable, public:: fixed_node_ids(:)                     !! node ids
        integer, allocatable, public:: k_type(:)                        !! k-type = 0 (at all steps), 1(static step), 2(dynamic step)
        real( kind = r_kind ), allocatable , public :: disp_vector(:,:)
        real( kind = r_kind ), allocatable , public :: spring_stiffness(:,:)
        contains 
           
               procedure , public  , pass :: initialize_boundary_spring
               procedure , public  , pass :: free_boundary_spring
           
    end type t_boundary_spring
    
    
!!=====================================================================! 
!!  nodal loads

    type, public :: t_nodal_load
        integer, public :: num_point_load = 0                           !! 
        integer, allocatable, public:: node_ids(:)                     !! node ids
        real( kind = r_kind ), allocatable , public :: load_vector(:,:)
        contains 
           
               procedure , public  , pass :: initialize_nodal_load
               procedure , public  , pass :: free_nodal_load
           
    end type t_nodal_load


!!=====================================================================! 
!!  load condition
    
    type, public :: t_load_case
        integer , public :: num_rows = 0  ;
        real( kind = r_kind ) :: tmin ;
        real( kind = r_kind ) :: tmax ;
        real( kind = r_kind ), public  :: load_coeffs(7) ;
        real( kind = r_kind ), allocatable , public :: time_load(:,:)   !! (time , load)
        contains 
                procedure , public  , pass :: initialize_load_case
                procedure , public  , pass :: free_load_case
    end  type t_load_case
    


!!====================================================================!!
!! DRM Layer

    type, public :: t_drm_layer
        integer  :: layer_id = 0 ;
        integer  :: num_element_in_layer = 0 ;
        integer  :: num_levels           = 0 ; ! levels at which EQ was deconvoluted (for drm_ndyn = 2)
        integer  :: num_node_in_layer    = 0 ;
        integer, allocatable  :: element_index(:)
        integer, allocatable  :: internal_face(:,:)
        integer, allocatable  :: paired_node_index(:,:)
        integer(kind = r_kind) , allocatable  :: levels(:)  ! levels at which EQ was deconvoluted (for drm_ndyn = 2)
        integer :: num_rows  = 0 ;                          !! Num Rows in DRM-NDYN = 2
        integer :: num_cols  = 0 ;                          !! Num cols in DRM-NDYN = 2
        contains 
                procedure , public  , pass :: initialize_drm_layer
                procedure , public  , pass :: free_drm_layer
    end type t_drm_layer

      
!!=====================================================================!!
  

   public :: get_nodes_on_element_face
   
   
   public :: local_coordinate_on_face

  
   public :: is_member
   
   
   
  contains
  

    

!!====================================================================!!
!!
!!  procedures of t_drm_layer
!!
!!
!!=====================================================================!!

 subroutine initialize_drm_layer( this )
    class ( t_drm_layer ) , intent( inout) :: this
    integer::  k
     
    if( this%num_element_in_layer > 0 ) then
        allocate( this%element_index( this%num_element_in_layer ) , stat = k )
        if( k /= 0) stop 'Failed in allocating memory in t_drm_layer class in utils module.'
        
        allocate( this%internal_face( this%num_element_in_layer ,2 ) , stat = k )
        if( k /= 0) stop 'Failed in allocating memory in t_drm_layer class in utils module.'
        
        
        if( this%num_node_in_layer > 0 ) then
            allocate( this%paired_node_index( this%num_node_in_layer , 2) , stat = k )
            if( k /= 0) stop 'Failed in allocating memory in t_drm_layer class in utils module.'
        end if 
    end if 
     
  end subroutine initialize_drm_layer
  
  
  
  subroutine free_drm_layer( this )
    class ( t_drm_layer ) , intent( inout) :: this
    integer::  k = 0
     
    if( allocated(this%element_index) ) deallocate(this%element_index ,stat=k)
    if( k /= 0) stop 'Failed in deallocating memory in t_drm_layer class in utils module.'
    
    if( allocated(this%internal_face) ) deallocate(this%internal_face ,stat=k)
    if( k /= 0) stop 'Failed in deallocating memory in t_drm_layer class in utils module.' 
     
    
    if( allocated( this%paired_node_index)) deallocate( this%paired_node_index , stat = k )
    if( k /= 0) stop 'Failed in deallocating memory in t_drm_layer class in utils module.'
     
    if( allocated( this%levels)) deallocate( this%levels , stat = k )
    if( k /= 0) stop 'Failed in deallocating memory in t_drm_layer class in utils module.'
      
     
  end subroutine free_drm_layer
   
    
  
!!====================================================================!!
!!
!!  procedures of t_model_description class
!!  boundary_id = 1 : reservior
!!  boundary_id = 2 : foundation
!!
!!=====================================================================!!

  subroutine initialize_model_description( this , boundary_id )
     class ( t_model_description ) , intent( inout) :: this
     integer , intent( in ) :: boundary_id 
     integer:: alloc_stat , k
     
     
     !! reservior truncation boundary: 0 for Somerfel, 1 for HW , 2 for GN
     if( ( this%reservior_inf_boundary_model_id > 0 ) .and. ( boundary_id .eq. 1 ) ) then   
         
         !! check the order(N,M) of the model: both GN and HW are allowed to take O(0,0)
         k = this%order_of_reservior_higher_order_bc(1) + 1             ! \{ a_j \} vector beginning from j = 0 , 1 , ... , N
             
         if( this%reservior_inf_boundary_model_id .eq. 2 ) then         ! GN model: there is no a_0
             k = k -1 
         end if
         
         if( k > 0 ) then
             allocate( this%reservior_incident_coefficient(k) , stat = alloc_stat )
             if( alloc_stat /= 0) stop 'Failed in allocating memory in t_model_description class in utils module.'
         else
             this%order_of_reservior_higher_order_bc(1) = 0 
         end if
         
         !! the evanescent part: in both models it can be zero
         k = this%order_of_reservior_higher_order_bc(2)                 ! \{b_j \} vector beginning from j = 1 , M
         
         if( k > 0 ) then
             allocate( this%reservior_evanescent_coefficient(k) , stat = alloc_stat )
             if( alloc_stat /= 0) stop 'Failed in allocating memory in t_model_description class in utils module.'
         else
             this%order_of_reservior_higher_order_bc(2) = 0 
         end if 
             
     end if
     
     
     
     !! foundation truncation boundary: 0 for Lysmer, 1 for HW
    if( ( this%foundation_inf_boundary_model_id > 0 ) .and. ( boundary_id .eq. 2 ) ) then      
         !! check the order(N,M) of the model: HW is allowed to take O(0,0)
         k = this%order_of_foundation_higher_order_bc(1) + 1            ! \{ a_j \} vector beginning from j = 0 , 1 , ... , N
             
         if( ( this%foundation_inf_boundary_model_id > 1 ) .or.  &
           & ( this%foundation_inf_boundary_model_id < 0 ) ) then       ! Invalid model id
    stop 'Invalid model id at foundation truncation boundary, choose either Lysmer (id =0) or HW (id=1)' 
         end if
         
         
         allocate( this%foundation_incident_coefficient(k) , stat = alloc_stat )
         if( alloc_stat /= 0) stop 'Failed in allocating memory in t_model_description class in utils module.'
         
         
         !! the evanescent part: in both models it can be zero
         k = this%order_of_foundation_higher_order_bc(2)                ! \{b_j \} vector beginning from j = 1 , M
         
         if( k > 0 ) then
             allocate( this%foundation_evanescent_coefficient(k) , stat = alloc_stat )
             if( alloc_stat /= 0) stop 'Failed in allocating memory in t_model_description class in utils module.'
         else
             this%order_of_foundation_higher_order_bc(2) = 0 
         end if 
             
    end if
    
              
    end subroutine initialize_model_description
   
   
   

    subroutine free_model_description( this )
        class ( t_model_description ) , intent( inout) :: this
        integer :: state = 0 , state1 = 0
        if(allocated(this%reservior_incident_coefficient)) deallocate(this%reservior_incident_coefficient,stat=state)
        state1 = state1 + state 
        if(allocated(this%reservior_evanescent_coefficient)) deallocate(this%reservior_evanescent_coefficient ,stat=state)
        state1 = state1 + state 
        if(allocated(this%foundation_incident_coefficient)) deallocate(this%foundation_incident_coefficient,stat=state)
        state1 = state1 + state 
        if(allocated(this%foundation_evanescent_coefficient)) deallocate(this%foundation_evanescent_coefficient,stat=state)
        state1 = state1 + state      
   
        if( state1 /= 0 ) stop 'Failed in deallocating memory in t_model_description class in utils module'
   
    end subroutine free_model_description
  
  
!!=====================================================================!!
!!
!! initialize enforced displacement type
!!
!! 
!!====================================================================!!

    subroutine initialize_enforced_disp(this)
        class ( t_enforced_disp ) , intent( inout) :: this
        integer:: alloc_stat = 0 ;
        
        if(allocated(this%node_ids)) deallocate(this%node_ids,stat=alloc_stat)
        if(allocated(this%disp_vector)) deallocate(this%disp_vector,stat=alloc_stat)
        if( alloc_stat /= 0 ) stop 'Failed in deallocating memory in t_enforced_dispn class in utils module'
        
        if( this%num_nodes > 0 ) then
            allocate( this%node_ids( this%num_nodes ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_enforced_dispn class in utils module.'
            
            allocate( this%disp_vector( this%num_nodes , 3 ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_enforced_dispn class in utils module.'
        end if
    
    
    end subroutine initialize_enforced_disp
  
  
!!=====================================================================!!
!!
!! free enforced displacement type
!!
!! 
!!====================================================================!!

    subroutine free_enforced_disp(this)
        class ( t_enforced_disp ) , intent( inout) :: this
        integer :: state = 0 
    
        if(allocated(this%node_ids)) then
           deallocate(this%node_ids,stat=state)
           if( state /= 0 ) stop 'Failed in deallocating memory in t_enforced_dispn class in utils module'
        end if
        
        if(allocated(this%disp_vector)) then
            deallocate(this%disp_vector,stat=state)
            if( state /= 0 ) stop 'Failed in deallocating memory in t_enforced_dispn class in utils module'
        end if
        
    end subroutine free_enforced_disp
    
    
!!=====================================================================!!
!!
!! initialize boundary spring type
!!
!! 
!!====================================================================!!

    subroutine initialize_boundary_spring(this)
        class ( t_boundary_spring ) , intent( inout) :: this
        integer:: alloc_stat = 0 ;
        
        if(allocated(this%fixed_node_ids))   deallocate(this%fixed_node_ids,stat=alloc_stat)
        if(allocated(this%k_type) )          deallocate(this%k_type,stat=alloc_stat)
        if(allocated(this%disp_vector))      deallocate(this%disp_vector,stat=alloc_stat)
        if(allocated(this%spring_stiffness)) deallocate(this%spring_stiffness,stat=alloc_stat)
        
        if( alloc_stat /= 0 ) stop 'Failed in deallocating memory in t_boundary_spring class in utils module'
        
        if( this%num_fixed_nodes > 0 ) then
            allocate( this%fixed_node_ids( this%num_fixed_nodes ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_boundary_spring class in utils module.'
            
            allocate( this%k_type( this%num_fixed_nodes ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_boundary_spring class in utils module.'
            
            
            allocate( this%disp_vector( this%num_fixed_nodes , 3 ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_boundary_spring class in utils module.'
            
            allocate( this%spring_stiffness( this%num_fixed_nodes , 3 ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_boundary_spring class in utils module.'
        end if
    
    
    end subroutine initialize_boundary_spring
  
    

!!=====================================================================!!
!!
!! free boundary spring type
!!
!! 
!!====================================================================!!

    subroutine free_boundary_spring(this)
        class ( t_boundary_spring ) , intent( inout) :: this
        integer :: state = 0 , state1 = 0
    
        if(allocated(this%fixed_node_ids))  then
           deallocate(this%fixed_node_ids,stat=state)
           if( state /= 0 ) stop 'Failed in deallocating memory in t_boundary_spring class in utils module'
        end if
        
        
        if(allocated(this%k_type ))   then
            deallocate(this%k_type, stat=state)
            if( state /= 0 ) stop 'Failed in deallocating memory in t_boundary_spring class in utils module'
        end if
        
        if(allocated(this%disp_vector))  then
            deallocate(this%disp_vector,stat=state)
            if( state /= 0 ) stop 'Failed in deallocating memory in t_boundary_spring class in utils module'
        end if
        
        if(allocated(this%spring_stiffness)) then
            deallocate(this%spring_stiffness,stat=state)
            if( state /= 0 ) stop 'Failed in deallocating memory in t_boundary_spring class in utils module'
        end if
        
        
    end subroutine free_boundary_spring
    
    
!!=====================================================================!!
!!
!! initialize nodal loads
!!
!! 
!!====================================================================!!

    subroutine initialize_nodal_load( this )
        class ( t_nodal_load ) , intent( inout) :: this
        integer:: alloc_stat = 0 ;
        
        if(allocated(this%node_ids))    deallocate(this%node_ids   , stat=alloc_stat )
        if(allocated(this%load_vector)) deallocate(this%load_vector, stat=alloc_stat )
        
        if( alloc_stat /= 0 ) stop 'Failed in deallocating memory in t_nodal_load class in utils module'
        
        if( this%num_point_load > 0 ) then
            allocate( this%node_ids( this%num_point_load ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_nodal_load class in utils module.'
            
            allocate( this%load_vector( this%num_point_load , 3 ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_nodal_load class in utils module.'
        end if
    end subroutine initialize_nodal_load
      
      
!!=====================================================================!!
!!
!!  free nodal loads
!! 
!!====================================================================!!

    subroutine free_nodal_load( this )
        class ( t_nodal_load ) , intent( inout ) :: this
        integer :: state = 0 , state1 = 0
    
        if(allocated(this%node_ids)) then
            deallocate( this%node_ids , stat = state )
            if( state /= 0 ) stop 'Failed in deallocating memory in t_nodal_load class in utils module'
        end if
       
       
        if(allocated(this%load_vector)) then
           deallocate( this%load_vector , stat = state )
           if( state1 /= 0 ) stop 'Failed in deallocating memory in t_nodal_load class in utils module'
        end if
        
        
    end subroutine free_nodal_load
      
    

!!=====================================================================!!
!!
!! initialize load conditions
!!
!! 
!!====================================================================!!
    
    

    subroutine initialize_load_case( this )
        class ( t_load_case ) , intent( inout) :: this
        integer:: alloc_stat = 0 ;
        
        if(allocated(this%time_load )) deallocate(this%time_load   , stat=alloc_stat )
        if( alloc_stat /= 0 ) stop 'Failed in deallocating memory in t_load_condition class in utils module'
        
        if( this%num_rows > 0 ) then
            allocate( this%time_load( this%num_rows , 8 ) , stat = alloc_stat )
            if( alloc_stat /= 0) stop 'Failed in allocating memory in t_load_condition class in utils module.'
            this%time_load(:,:) = 0.0_r_kind ;
        end if
    end subroutine initialize_load_case
      
    
!!=====================================================================!!
!!
!!  free load conditions
!! 
!!====================================================================!!

    subroutine free_load_case( this )
        class ( t_load_case ) , intent( inout ) :: this
        integer :: state 
        
        if(allocated(this%time_load )) then
            deallocate( this%time_load , stat = state )
            if( state /= 0 ) stop 'Failed in deallocating memory in t_load_conditiond class in utils module'
        end if 
        
    end subroutine free_load_case
      
     
    
  
!!====================================================================!!
!!
!!  get local coordinate for given face
!!
!!====================================================================!!
  
  subroutine  local_coordinate_on_face( n_face , vector )
     implicit none
     integer , intent( in    ) :: n_face 
     integer , intent( inout ) :: vector(2)
     
     
     select  case ( n_face ) 
              case ( 1  )                   ! [xi , eta , zeta ] = [ 1 , ... , ...]
                   vector = [ 1 , +1 ]      ! xi is dim = 1  and values = +1 so r = [dim , value ] = [ 1 , 1]                 
              case ( 2  )
                   vector = [ 1 , -1 ]      ! xi = -1  -->  [dim , value]  = [ 1 , -1]            
              case ( 3  )            
                   vector = [ 2 , +1 ]      ! eta = +1 -->  [dim , value]  = [ 2 , 1]                 
              case (  4 )     
                   vector = [ 2 , -1 ]      ! eta = -1 -->  [dim , value]  = [ 2 , -1]               
              case ( 5 ) 
                   vector = [ 3 , +1 ]      ! zeta = +1 --> [dim , value] = [ 3 , +1]                 
              case ( 6 ) 
                   vector = [ 3 , -1 ]      ! zeta = -1 --> [dim , value] = [ 3 , -1]   
              case default
                  write(*,*) 'face_id = ' , n_face
                  stop 'Error: face id is out of bound in face_to_node_index in utils module '
      end select
  
  end subroutine
 
 
!!============================================================================!!
!! get nodes that lie on a face of an element
!! input: n_face     : the id of face under consideration
!!        n_dim      : the element type
!!        node_index : the index of nodes that lie on the face
!!
!! return: n     : num of non-zero elements of node_index
!!============================================================================!!

  
   function  get_nodes_on_element_face( n_face , the_node_parent , node_index  )  result ( n_node )  
      implicit none
      integer , intent( in    )  :: n_face 
      integer , dimension(:,:) , intent( in )  :: the_node_parent
      integer , intent( inout )  :: node_index(:)
      integer :: i , r(2) , n_node 
      
      node_index = 0
     
      
      call  local_coordinate_on_face( n_face , r)
          
      
      !-----------------------------------------------------------------!
      n_node = 0
    
      do i = 1 , size( the_node_parent  , 2 )
         if( the_node_parent( r(1) , i ) .eq. r(2) ) then 
             n_node = n_node + 1
             node_index( n_node ) = i 
         end if
      end do 
 
   end function get_nodes_on_element_face
   
   

  
  
  
  
  

!!=====================================================================!!
!!
!!  check if element x belong to a vector up to elements of length +1, i.e. the last element is included. 
!!  The set may be accendingly ordered or not.
!!  
!!====================================================================!!
   
  function  is_member( vector, a_number , length , is_ordered )  result( answer )
      
      implicit none
      integer , dimension(:) , intent( in ) :: vector
      integer , intent( in ) :: a_number
      integer , intent(in ) , optional :: length
      logical , intent(in) , optional :: is_ordered                      ! default is false
      logical :: answer
      
      integer :: i , until                                       
      
      
      answer = .false.
      if( present(length ) ) then
          until = length
      else
         until  = size( vector )
      end if
      
      if( until .eq. 0 ) return  
      
      if( present( is_ordered ) ) then
          if( is_ordered ) then
              if( ( a_number < vector(1) ) .or. &                       ! smaller than first element
                & ( a_number > vector(until) )  ) return                ! larger than last element
          end if
      end if
      
      
      
      do i = 1 , until
         if( a_number .eq. vector(i) ) then
             answer = .true.
             return
         end if
      end do 
      
   
  end function  is_member
  
  
  
  
  
  
  
  
  
  
  
  
  
  
end module mod_utils

