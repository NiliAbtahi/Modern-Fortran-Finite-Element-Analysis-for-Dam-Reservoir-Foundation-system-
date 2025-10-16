!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!


!
!   The base element class that defines basic variables required for different types of elements.
!   Each element is defined in separate module.
!
!   There are only five things that should be specilized for each element type:
!   
!        1. num_main_nodes        
!        2. num_node_in_element
!        3. node parent matrix
!        4. shape vector
!        5. local_d_shape_matrix  (matrix of derivative of shape vector in local coordinates)
!
!    Explanation:
!        1. num_main_nodes 
!           Main nodes are minimum number of nodes that are required for defining an element of desired type.
!                Example 1. for quadrilaterial element Q9 the four first nodes are basic (main nodes)
!                Example 2. for traingular element T8 the three first nodes are basic, otherwise a triangle cannot be defined.
!                Eaxmple 3. for brick element H20 or H27 only 6 first nodes are main nodes.
!
!        Remark: main nodes comes first, that is the first n nodes are taken as basic.
!                Main nodes should not be zero, they should refer to index of existing node.
!
!
!        2. num_node_in_element
!           The number of all nodes in an element, that is both main and optional nodes.
!           Example 1. In a Q9 element, the four main nodes are necessary while each of 5th to 9th nodes can be absent.
!           Example 2. In brick H20 element, the first 6 nodes are necessary while each of 7th to 20th nodes can be absent,
!           
!           In example 1, num_node_in_element = 9  even if some of nodes are absent.
!           In example 2, num_node_in_element = 20 even if some of nodes are absent.
!
!
!       3. node parent matrix
!          The node parent matrix is required for every type of nodes.
!          Examples are clear for Q8 and H20 elements.
!
!       4. shape vector
!          To each subclass of base-element (,i.e. each element type) the shape vector should be defined.
!          
!       5. local_d_shape_matrix
!          To each subclass of base-element (,i.e. each element type) the matrix of derivative of shape
!          functions in local coordinate should be defined.
!
!      The five items defined above depend on the choice of element type.
!      Every other property and methods defined in this class is in common for all element types.
!
!
!     By pass feature, the first (hidden) argument of the procedures is the object itself.
!!==========================================================================================================!!




 module mod_element
     
     use mod_utils , only : r_kind
     implicit none
      
     
     type, abstract , public ::  t_element
     
         
         integer , public  :: num_main_nodes                  
         integer , public  :: num_node_in_element 
         integer , public  :: num_element_around_node
         integer , public  :: num_face
                    
         
         integer , pointer , public  :: node_parent(:,:)     
            
         real(   kind = r_kind ) , allocatable , public  :: shape_vector(:)              ! elements of shape matrix in a vector
         real(   kind = r_kind ) , allocatable , public  :: derivative_shape_matrix(:,:) ! deriavtive of shape vector in global coordinate     
         real(   kind = r_kind ) , allocatable , public  :: jacobian_matrix(:,:)         ! jacobian of transformaion between local and global coordinates
         real(   kind = r_kind )               , public  :: determinant_of_jacobian      ! determinant of jacoban matrix
         real(   kind = r_kind ) , allocatable , public  :: normal_at_face(:)
         
         real(kind = r_kind) , allocatable , private  :: local_d_shape_matrix(:,:)    ! deriavtive of shape vector in local coordinate       
         real(kind = r_kind) , allocatable , private  :: inverse_jacobian_matrix(:,:) ! inverse of jacobian matrix
         
        
         
         
         
     contains
        
        procedure , private  , pass :: base_element_constructor         
        procedure , private  , pass :: are_nodes_related
        procedure , private  , pass :: get_jacobian_matrix
        procedure , private  , pass :: get_inverse_jacobian_matrix
        procedure , private  , pass :: add_optional_node_share           
        
        procedure , public  , pass :: destruct_element                         
        procedure , public  , pass :: show_element  
        procedure , public  , pass :: get_normal_at_face  
        procedure , public  , pass :: volume_element_at_face 
        procedure , public  , pass :: get_derivative_shape_function
             
   
        procedure ( element_constructor               ) , public  , deferred  :: construct_element          
        procedure ( vector_of_shape_function          ) , public  , deferred  :: get_shape_function 
        procedure ( matrix_of_derivatve_shape_function) , private , deferred  :: local_derivative_of_shape_function
            
     end type t_element
   
!!=====================================================================!!

     abstract interface 
       subroutine  element_constructor( this )
            
             import t_element
             class( t_element )  , intent( inout ) :: this
             
       end subroutine
     end interface  

!!=====================================================================!!

     abstract interface 
       subroutine  vector_of_shape_function( this , elem_id , xi , eta , zeta )
            
             use mod_utils , only: r_kind
             import t_element
             
             
             class( t_element )    , intent( inout ) :: this
             integer               , intent( in    ) :: elem_id 
             real( kind = r_kind ) , intent( in    ) :: xi
             real( kind = r_kind ) , intent( in    ) :: eta
             real( kind = r_kind ) , intent( in ) , optional :: zeta     ! geometry can be 2D or 3D
             
       end subroutine
     end interface
     

!!=====================================================================!!

     abstract interface
       subroutine  matrix_of_derivatve_shape_function( this ,elem_id , xi , eta , zeta )
       
             use mod_utils , only: r_kind
             import t_element
             
             
             class( t_element )   , intent( inout ) :: this
             integer              , intent( in    ) :: elem_id 
             real( kind = r_kind) , intent( in    ) :: xi
             real( kind = r_kind) , intent( in    ) :: eta
             real( kind = r_kind) , intent( in ) , optional :: zeta     ! geometry can be 2D or 3D
       
       end subroutine
     end interface
     
     
!=======================================================================!


!!=======================================================================
!!
!!  Q8 element
!!
!!=======================================================================  
    
     type , public , extends( t_element ) ::  t_Q8_element
          
          contains
          procedure , public , pass :: construct_element                => Q8_constructor
          procedure , public , pass :: get_shape_function               => Q8_shape_function
          procedure , public , pass :: local_derivative_of_shape_function    => Q8_local_derivative_shape_function
 
     end type t_Q8_element   
    
      
!!=======================================================================
!!
!!  Brick element
!!
!!=======================================================================  
    
     type , public , extends( t_element ) ::  t_brick_element
          
          contains
          procedure , public , pass :: construct_element                => brick_constructor
          procedure , public , pass :: get_shape_function               => brick_shape_function
          procedure , public , pass :: local_derivative_of_shape_function    => brick_local_derivative_shape_function
          
     end type t_brick_element 

     
  



!!=======================================================================
!!
!!  declare other element type here, as two examples above
!!
!!=======================================================================  
  


  contains
    


 
        
!!=======================================================================
!!  ( Generic for every kind of element)
!!
!!
!!  base class : initialize element class
!!
!!=======================================================================
   
  subroutine  base_element_constructor( the_element )
       
       use mod_geometry , only : num_dim
       implicit none   
       class( t_element ) , intent( inout) :: the_element
       
       integer :: state
       
       
       call the_element%destruct_element()                              ! check if the element arrays are null
                          
       
       ! allocate memory for shape vector
       allocate( the_element%shape_vector( the_element%num_node_in_element ) , stat = state)  
       if ( state /= 0)  then
            call the_element%destruct_element()  
            stop 'Error: failed in allocating shape vector array in element class.'  
       end if
       
       ! allocate memory for derivative_shape_matrix in global coordinate system
       allocate( the_element%derivative_shape_matrix( num_dim , the_element%num_node_in_element ) , stat = state)  
       if ( state /= 0)  then
            call the_element%destruct_element()
            stop 'Error: failed in allocating derivative_shape_matrix in element class.'  
       end if
       
       ! allocate memory for derivative_shape_matrix in local coordinate system
       allocate( the_element%local_d_shape_matrix(num_dim , the_element%num_node_in_element ) , stat = state)  
       if ( state /= 0)  then
            call the_element%destruct_element()
            stop 'Error: failed in allocating derivative_shape_matrix in element class.'  
       end if
       
       ! allocate memory for jacobian_matrix
       allocate( the_element%jacobian_matrix(num_dim , num_dim ) , stat = state)  
       if ( state /= 0)  then
            call the_element%destruct_element()
            stop 'Error: failed in allocating Jacobian matrix in element class.'  
       end if
       
       ! allocate memory for inverse of jacobian_matrix
       allocate( the_element%inverse_jacobian_matrix(num_dim , num_dim ) , stat = state)  
       if ( state /= 0)  then
            call the_element%destruct_element()
            stop 'Error: failed in allocating inverse Jacobian matrix in element class.'  
       end if
       
       allocate( the_element%normal_at_face(num_dim ) , stat = state)  
       if ( state /= 0)  then
            call the_element%destruct_element()
            stop 'Error: failed in allocating normal_to_face array in element class.'  
       end if
       
       
             
  end subroutine base_element_constructor
    
      
     
       
!!=======================================================================
!!  ( Generic for every kind of element)
!!
!!
!!  base class: clean up element class
!!
!!=======================================================================
   
  subroutine  destruct_element( the_element )
       
       implicit none
       class( t_element ) , intent( inout) :: the_element
       
       integer :: state , state1 = 0
       
       
       if( associated( the_element%node_parent ) ) then
           nullify( the_element%node_parent )
       end if
       
       
       if( allocated( the_element%shape_vector )) then
           deallocate( the_element%shape_vector , stat = state)
           state1 = state1 + state
       end if
       
       if( allocated( the_element%derivative_shape_matrix ) ) then
           deallocate( the_element%derivative_shape_matrix , stat = state )
           state1 = state1 + state
       end if
       
       if( allocated( the_element%local_d_shape_matrix ) ) then
           deallocate( the_element%local_d_shape_matrix , stat = state )
           state1 = state1 + state
       end if
       
       if( allocated( the_element%jacobian_matrix ) ) then
           deallocate( the_element%jacobian_matrix , stat = state )
           state1 = state1 + state
       end if
       
       if( allocated( the_element%inverse_jacobian_matrix ) ) then
           deallocate( the_element%inverse_jacobian_matrix , stat = state )
           state1 = state1 + state
       end if
       
       if( allocated( the_element%normal_at_face ) ) then
           deallocate( the_element%normal_at_face , stat = state )
           state1 = state1 + state
       end if
       
      
       
       ! check if things are done correctly
       if( state1 /= 0 ) stop 'Error: failed in deallocating arrays in element class.'  
                  
  end subroutine destruct_element

   
 
  
     
!!=======================================================================
!!  ( Generic for every kind of element)
!!
!!
!!  modify shape function and its derivative to include effect of optional nodes
!!
!! Note: for elements other than Q8 and brick the factor 2 in denum of following algorithm may be changed.
!!       If so, it needs to include the suitable factor into the subroutine arguments.
!!=======================================================================

   subroutine  add_optional_node_share( this , elem_id , to_shape_function )
      
      use mod_geometry, only : element_matrix
      implicit none
      
      class( t_element )  , intent( inout ) :: this
      integer             , intent( in    ) :: elem_id
      logical             , intent( in    ) :: to_shape_function        ! true means shape function, false its derivative
      integer  :: i , k
      
      do i = this%num_main_nodes + 1 , this%num_node_in_element
       
         if( element_matrix( elem_id , i ) .eq. 0 ) cycle               ! this node does not exist
          
         ! loop over base nodes                                                  
         do k = 1 , this%num_main_nodes 
                                              
            if( .not. this%are_nodes_related( k , i )) cycle  
            
            ! if shape function is going to be modified
            if( to_shape_function ) then
                this%shape_vector(k) = this%shape_vector(k) - 0.50d0 * this%shape_vector(i)
            else
                this%local_d_shape_matrix(:,k) =  this%local_d_shape_matrix(:,k) &
                                               -  this%local_d_shape_matrix(:,i) / 2.0d0
            end if
            
         end do ! loop over k  
                                        
      end do    ! loop over i
      
   
   end subroutine add_optional_node_share
   
     
!!=======================================================================
!!  ( Generic for every kind of element)
!!
!!
!!  base class : find Jacobian matrix and its inverse
!!
!!=======================================================================
  
 subroutine get_jacobian_matrix(this, elem_id )  
 
     use mod_geometry, only : num_dim , element_matrix , node
     implicit none
     
     class(t_element) , intent( inout ) :: this
     integer          , intent( in    ) :: elem_id 
     
     integer :: col , n 
   
   
 !!  finding Jacobian matrix    
     
     this%jacobian_matrix = 0.0_r_kind
     
     do col = 1 , num_dim  
        do n = 1 , this%num_node_in_element
        if( element_matrix( elem_id , n ) .eq. 0 ) cycle                ! if the node does not exist
        
        this%jacobian_matrix(:, col) = this%jacobian_matrix(: , col)  &
      + this%local_d_shape_matrix(: , n) * node(element_matrix( elem_id , n) , col )
           
        end do ! over n
     end do ! over col
     
  
 end subroutine get_jacobian_matrix
 
 
!!=====================================================================!!
!!  ( Generic for every kind of element)
!!
!!
!!  get inverse Jacobian matrix
!!
!!=====================================================================!!
 
  subroutine  get_inverse_jacobian_matrix( this )
  
     use mod_geometry , only : num_dim
     implicit none
     
     class(t_element) , intent( inout ) :: this
     real( kind = r_kind ) :: det 
 

    !! determinant of matrix and its inverse if num_dim = 2
    if( num_dim .eq. 2 ) then
        det  = this%jacobian_matrix(1,1) * this%jacobian_matrix(2,2) &
             - this%jacobian_matrix(1,2) * this%jacobian_matrix(2,1)
        
        ! check if matrix is non-singular
        if( abs( det ) .le. 1000.0_r_kind * epsilon( det )) stop ' Error: singular Jacobian matrix in element module' 
        
        this%determinant_of_jacobian = det 
        
        this%inverse_jacobian_matrix(1,1) =  this%jacobian_matrix(2,2) / det 
        this%inverse_jacobian_matrix(1,2) = -this%jacobian_matrix(1,2) / det 
        this%inverse_jacobian_matrix(2,1) = -this%jacobian_matrix(2,1) / det 
        this%inverse_jacobian_matrix(2,2) =  this%jacobian_matrix(1,1) / det 
        
        return

    end if
    
    
    !! if code gets here then num_dim = 3, first find adjugate matrix
    
  this%inverse_jacobian_matrix(1,1) = this%jacobian_matrix(2,2) * this%jacobian_matrix(3,3)  &
                                    - this%jacobian_matrix(2,3) * this%jacobian_matrix(3,2)
    
  this%inverse_jacobian_matrix(2,1) = this%jacobian_matrix(2,3) * this%jacobian_matrix(3,1)  &
                                    - this%jacobian_matrix(2,1) * this%jacobian_matrix(3,3)
    
  this%inverse_jacobian_matrix(3,1) = this%jacobian_matrix(2,1) * this%jacobian_matrix(3,2)  &
                                    - this%jacobian_matrix(2,2) * this%jacobian_matrix(3,1)
                                        
  this%inverse_jacobian_matrix(1,2) = this%jacobian_matrix(1,3) * this%jacobian_matrix(3,2)  &
                                    - this%jacobian_matrix(1,2) * this%jacobian_matrix(3,3)
    
  this%inverse_jacobian_matrix(2,2) = this%jacobian_matrix(1,1) * this%jacobian_matrix(3,3)  &
                                    - this%jacobian_matrix(1,3) * this%jacobian_matrix(3,1)
    
  this%inverse_jacobian_matrix(3,2) = this%jacobian_matrix(1,2) * this%jacobian_matrix(3,1)  &
                                    - this%jacobian_matrix(1,1) * this%jacobian_matrix(3,2)  

  this%inverse_jacobian_matrix(1,3) = this%jacobian_matrix(1,2) * this%jacobian_matrix(2,3)  &
                                    - this%jacobian_matrix(1,3) * this%jacobian_matrix(2,2)
    
  this%inverse_jacobian_matrix(2,3) = this%jacobian_matrix(1,3) * this%jacobian_matrix(2,1)  &
                                    - this%jacobian_matrix(1,1) * this%jacobian_matrix(2,3)
                                          
  this%inverse_jacobian_matrix(3,3) = this%jacobian_matrix(1,1) * this%jacobian_matrix(2,2)  &
                                    - this%jacobian_matrix(1,2) * this%jacobian_matrix(2,1)
    
                        
     ! get determinant
     det  =  this%jacobian_matrix(1,1) * this%inverse_jacobian_matrix(1,1)  &
          +  this%jacobian_matrix(1,2) * this%inverse_jacobian_matrix(2,1)  &
          +  this%jacobian_matrix(1,3) * this%inverse_jacobian_matrix(3,1)    
     
     ! check if matrix is singular
     if( abs( det ) .le. 1000.0_r_kind * epsilon( det  )) stop ' Error: singular Jacobian matrix in element module' 
        
     this%determinant_of_jacobian = det 
     
     ! divide adjugate matrix by determinant
     this%inverse_jacobian_matrix = this%inverse_jacobian_matrix / det 
     
     
  
  end subroutine get_inverse_jacobian_matrix
        
!!====================================================================!!
!!  ( Generic for every kind of element)
!!
!!
!! matrix containing derivative of shape functions in global coordinate
!! This function is called to get derivative matrix of shape function.
!!
!!=====================================================================!!

  subroutine get_derivative_shape_function( this , elem_id , xi , eta , zeta)
     use mod_geometry , only : num_dim
     implicit none
     class(t_element) , intent( inout ) :: this
     integer              , intent( in    ) :: elem_id 
     real( kind = r_kind) , intent( in    ) :: xi
     real( kind = r_kind) , intent( in    ) :: eta
     real( kind = r_kind) , intent( in ) , optional :: zeta             ! geometry can be 2D or 3D
       
     
     if( num_dim .eq. 2 ) then
         call this%local_derivative_of_shape_function( elem_id , xi , eta )
     else 
         if( .not. present( zeta ) ) then
             stop ' Not enough argument in get_derivative_shape_function for 3d element: in element module'
         end if
         
         call this%local_derivative_of_shape_function( elem_id , xi , eta , zeta )
     end if
     
     call this%get_jacobian_matrix( elem_id )       !! find Jacobian matrix
      
     call this%get_inverse_jacobian_matrix()        !! find invarse Jacobian matrix
     
     this%derivative_shape_matrix = matmul( this%inverse_jacobian_matrix , &  !! derivative in global coordinate
                                             this%local_d_shape_matrix    )
  
  end subroutine get_derivative_shape_function
  
  
  
  
  
!!=======================================================================
!!
!!  show element information
!!
!!=======================================================================
 
  subroutine  show_element( the_element )
     
     implicit none
     class( t_element ) , intent( in ) :: the_element
     
     integer :: i 
     
     write(*,*) 'Node parent: '
     do i = 1 , size( the_element%node_parent , 1)
           write(*,*)  the_element%node_parent(i,:)  
     end do ! loop over i
     
     
     write( *,*) ' '
     write( *,*) ' '
     write( *,*) ' The shape function vector : '
     write( *,*) ' '
     write( *,*) ' '
     
     do i = 1 , size( the_element%shape_vector )
          write(*,*)  the_element%shape_vector(i) 
     end do
     
     write( *,*) ' '
     write( *,*) ' '
     write( *,*) ' Derivative of shape function in local coordinate ' 
     write( *,*) ' '
     write( *,*) ' '
     do i = 1 , size( the_element%local_d_shape_matrix , 1)
           write(*,*)  the_element%local_d_shape_matrix(i,:) 
     end do ! loop over i
     
     
     write( *,*) ' '
     write( *,*) ' '
     write( *,*) ' Jacobian matrix '
     write( *,*) ' '
     write( *,*) ' '
     do i = 1 , size( the_element%jacobian_matrix , 1)
        write(*,*) the_element%jacobian_matrix(i , :) 
     end do ! loop over i
     
     
     write( *,*) ' '
     write( *,*) ' '
     write( *,*) ' inverse Jacobian matrix '
     write( *,*) ' '
     write( *,*) ' '
     do i = 1 , size( the_element%jacobian_matrix , 1)
        write(*,*) the_element%inverse_jacobian_matrix(i , :) 
     end do ! loop over i
  
     write( *,*) ' '
     write( *,*) ' '
     write( *,*) ' Derivative of shape function in global coordinate ' 
     write( *,*) ' '
     write( *,*) ' '
     do i = 1 , size( the_element%derivative_shape_matrix , 1)
        write(*,*)  the_element%derivative_shape_matrix(i,:)  
     end do ! loop over i
     
     write( *,*) ' '
     write( *,*) ' '
     write( *,*) ' determinant = ' , the_element%determinant_of_jacobian
     
  end subroutine show_element
 
 
 
 
!!=======================================================================
!!  ( Generic for every kind of element)
!!
!!
!!  check if shape function of node id2 affect that of node id1 for arbitrary element type
!!
!!=======================================================================  
 
 logical function are_nodes_related( the_element , id1 , id2)  result (res)
 
     use mod_geometry, only : num_dim
     implicit none
     
     class( t_element  ) , intent(in) :: the_element
     integer             , intent(in) :: id1 , id2
 
      
     logical  :: mask_is_zero(num_dim) , mask_equality(num_dim) 
     integer  :: n
     
     res = .false.
     if( id1 .eq. id2 ) return
       
    
     mask_is_zero = the_element%node_parent( 1: num_dim , id2 ) .eq. 0
           
     ! if there is no zero coordinate then has no effect
     if( .not. any( mask_is_zero )) return
         
     
     ! equality of non-zero coordinates for the two nodes
      mask_equality = the_element%node_parent(: , id1) .eq. the_element%node_parent(: , id2)
      res = .true.
             
      ! check equality of non-zero coordinates
      do n = 1 , num_dim                            
         if( mask_is_zero(n) ) cycle
         res = res .and. mask_equality(n)                               ! multiplication of boolean variables
      end do 
      
          
 end function are_nodes_related

 

 
!!=======================================================================
!!  ( Generic for every kind of element)
!!
!!
!!  get normal to face
!!  direction = 1 : xi   direction
!!  direction = 2 : eta  direction
!!  direction = 3 : zeta direction
!!
!! In 2D case: the coordinate functions are r(xi , eta ) = x(xi , eta) e_x + y( xi , eta ) e_y
!! so in direction 1, xi = fixed , so we get
!! dr/d eta = ( dx/d eta ) d eta e_x + ( dy/d eta ) d eta e_y
!! This is tangent vector. The normal is 
!!      either n =  (dy/ d eta ) e_x - (dx/d eta) e_y,
!!          or n = -(dy/ d eta ) e_x + (dx/d eta) e_y,
!! Similar expression is valid for other direction.
!!
!! Jacobian matrix is also
!!                         J = [ x_{,xi  }  , y_{,xi  } , z_{,xi  }  ;
!!                               x_{,eta }  , y_{,eta } , z_{,eta }  ;
!!                               x_{,zeta}  , y_{,zeta} , z_{,zeta}  ; ]
!!
!!
!! In 3D element, one of local coordinates is fixed, making one row of
!! Jacobian matrix zero. We can find two other vectors on the boundary face.
!! The cross product of these vectors is used to get normal vector at a given point.
!! 
!! Assume that the xi is the fixed vector. In the eta-zeta plane, we define
!! the position vector as
!!      r(xi , eta , zeta ) = x( xi , eta , zeta) e_x + y( xi , eta , zeta) e_y + z( xi , eta , zeta) e_z,
!!
!! The two velocity vector are obained from derivative along eta-fixed and zeta-fixed directions.
!! 
!!     v_1 d eta = ( dr / d eta ) = ( dx / d eta ) e_x + ( dy / d eta ) e_y + ( dz / d eta ) e_z , 
!!     v_2 d zeta = ( dr / d zeta ) = ( dx / d zeta ) e_x + ( dy / d zeta ) e_y + ( dz / d zeta ) e_z ,
!!
!!     These vectors are second and third rows of the Jacobian matrix.
!!     Similarly if boundary face is located at  eta = const. then first and third  rows of Jacobian give the two vector.
!!     Similarly if boundary face is located at zeta = const. then first and second rows of Jacobian give the two vector.
!!
!!     The two vectors obtained in this way lie in the boundary face. 
!!     The normal vector is the cross product of these two vectors, up to size and sign.
!!     The size is fixed by making the vector unit.
!!     The sign is fixed by making is inward or outward.
!!
!!====================================================================== 
  
   subroutine  get_normal_at_face( this , elem_id , dir , facial_nodes , is_outward ) 
       
       use mod_geometry , only : num_dim , element_matrix , node
       use mod_utils    , only : is_member
       implicit none    
             
       class( t_element )   , intent( inout ) :: this
       integer , intent( in    ) :: dir , elem_id
       integer , dimension(:) , intent(in) :: facial_nodes
       logical , intent(in) , optional        :: is_outward             ! default is outward
       
       integer ::  i , k
       real( kind = r_kind ) :: v(3) , temp
    
                               
       
       !! First find normal vector up to size and a minus sign
       if( num_dim .eq. 2 ) then
           if( dir .eq. 1 ) then                                        ! xi is fixed
               this%normal_at_face(1) =  this%jacobian_matrix( 2 , 2 )
               this%normal_at_face(2) = -this%jacobian_matrix( 2 , 1 )
           else                                                         ! eta is fixed
               this%normal_at_face(1) =  this%jacobian_matrix( 1 , 2 )
               this%normal_at_face(2) = -this%jacobian_matrix( 1 , 1 )
           end if
       
       
       else
            if( dir .eq. 1 ) then      ! xi is fixed : cross product of second and third rows of Jacobian matrix.
            
                this%normal_at_face(1) = this%jacobian_matrix( 2 , 2 ) * this%jacobian_matrix( 3 , 3 ) - &
                                       & this%jacobian_matrix( 2 , 3 ) * this%jacobian_matrix( 3 , 2 )
                
                this%normal_at_face(2) = this%jacobian_matrix( 2 , 3 ) * this%jacobian_matrix( 3 , 1 ) - &
                                       & this%jacobian_matrix( 2 , 1 ) * this%jacobian_matrix( 3 , 3 )
                
                this%normal_at_face(3) = this%jacobian_matrix( 2 , 1 ) * this%jacobian_matrix( 3 , 2 ) - &
                                       & this%jacobian_matrix( 2 , 2 ) * this%jacobian_matrix( 3 , 1 )       
                                                       
            elseif ( dir .eq. 2 ) then   ! eta is fixed: cross product of first and third rows of Jacobian matrix.
            
                this%normal_at_face(1) = this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 3 , 3 ) - &
                                       & this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 3 , 2 )
                
                this%normal_at_face(2) = this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 3 , 1 ) - &
                                       & this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 3 , 3 )
                
                this%normal_at_face(3) = this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 3 , 2 ) - &
                                       & this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 3 , 1 )   
            
            else                         ! zeta is fixed: cross product of first and second rows of Jacobian matrix.
            
                this%normal_at_face(1) = this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 2 , 3 ) - &
                                       & this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 2 , 2 )
                
                this%normal_at_face(2) = this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 2 , 1 ) - &
                                       & this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 2 , 3 )
                
                this%normal_at_face(3) = this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 2 , 2 ) - &
                                       & this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 2 , 1 )  
            end if
       
       end if
       
       
       ! to fix size normalize vector to unit vector
       temp = 0.0_r_kind
       do i = 1 , num_dim
          temp = temp + this%normal_at_face(i) ** 2
       end do
       
       this%normal_at_face = this%normal_at_face / sqrt( temp )
       
       ! to fix sign (direction) get a non-facial node of the element.
       ! First get index of non-facial node

       k = 0
       do i = 1 , this%num_node_in_element
          if( element_matrix( elem_id , i ) .eq. 0 ) cycle              ! node is absent
          if( is_member( facial_nodes , i ) ) cycle                     ! node is facial
          k = i
          exit
       end do
       
       
       !! get coordinate of non-facial node
       v(1:num_dim ) = node( element_matrix( elem_id , k ) , : )
       
       !! get a facial node
       k = 0
       do i = 1 , size( facial_nodes )
          if(  element_matrix( elem_id , facial_nodes(i) ) .eq. 0 ) cycle             ! node is absent
          k = i
          exit
       end do 
       
       !! form displacement vector
       v(1:num_dim) = v(1:num_dim) - node( element_matrix( elem_id , facial_nodes(k) ) , : )
       
       !! form inner product of displacement vector and normal
       temp = 0.0_r_kind
       do i = 1 , num_dim
          temp = temp + v(i) * this%normal_at_face(i)
       end do
       
       
       !! if temp > 0 then normal is inward (angle < 90) , if temp < 0 then normal is outward ( angle > 90)
       !! The angle is never 90 by construction, hence temp = 0 is never obtained.
     
       if( present( is_outward ) ) then
       
           if( is_outward .and.  temp > 0 ) then
          
               this%normal_at_face =  -1.0_r_kind * this%normal_at_face 
              
           elseif( (.not. is_outward ) .and.  temp < 0 ) then
          
               this%normal_at_face =  -1.0_r_kind * this%normal_at_face 
              
           end if
          
       else                                                             ! default is outward
                                               
           if( temp > 0 ) then
               this%normal_at_face =  -1.0_r_kind * this%normal_at_face 
           end if
       end if
           

   end subroutine get_normal_at_face
 
 
!!====================================================================!!
!!  ( Generic for every kind of element)
!!
!!
!!  volume element at face
!!  If element is 2D then face is an 1D edge with length element 
!!  dL = sqrt( dx^2 + dy^2 )
!!  Also the change of variable (from global to local) is performed in integral.
!!
!!  If xi = fixed then dx = (dx/d eta ) d eta  and dy = (dy/d eta ) d eta, 
!!                  so dL = sqrt( (dx/d eta ) ^2 + (dy/d eta )^2) d eta
!!        The change of variable is then ( dL / d \eta ) or ( \partial L / \partial \eta ) 
!!        so the factor is = sqrt( (dx/d eta ) ^2 + (dy/d eta )^2)
!!   
!!  If eqta = fixed, then dx = (dx/d xi ) d xi  and dy = (dy/d xi ) d xi, 
!!                     so dL = sqrt( (dx/d xi ) ^2 + (dy/d xi )^2) d xi
!!         The change of variable is then ( dL / d \eta ) or ( \partial L / \partial \eta ) 
!!         so the factor is = sqrt( (dx/d xi ) ^2 + (dy/d xi )^2)
!!   
!!
!! If element is 3D then face is 2D. In isoparameteric coordinate system,
!! the boundary faces lie at constant values of one of coordinates xi , eta , zeta.
!! 
!! To get colume element we use the fact the area of parallelogram of a surface
!! spanned by two vectors is obtained from magnitude of the cross product of the two vectors.
!! So we define position vector as
!!              r(xi , eta , zeta ) = x( xi , eta , zeta ) e_x + y( xi , eta , zeta ) e_y  + z( xi , eta , zeta ) e_z ,
!!
!! Now we define two tangent vectors in the face xi = constant,
!! 
!! v_1 d eta = (dr / d eta)  , v_1 d zeta = (dr / d zeta) 
!! that are second and third rows of the Jacobian matrix evaluated at the point on the boundary face.
!!
!! In terms of local coordinates, we first have Jacobian matrix as
!! 
!!                         J = [ x_{,xi  }  , y_{,xi  } , z_{,xi  }  ;
!!                               x_{,eta }  , y_{,eta } , z_{,eta }  ;
!!                               x_{,zeta}  , y_{,zeta} , z_{,zeta}  ; ]
!!
!! Define , v1 = J(1 , :) , v2 = J(2 , :) , v3 = J(3 , : ) ,
!!
!!
!!  If xi   = const. then dA = \abs( cross_product( v2 , v3 ) ) d \eta  d\zeta
!!  If eta  = const. then dA = \abs( cross_product( v1 , v3 ) ) d \xi   d\zeta
!!  If zeta = const. then dA = \abs( cross_product( v1 , v2 ) ) d \xi   d\eta
!!
!!=====================================================================!!
 
 function  volume_element_at_face(this , elem_id , dir , xi , eta , zeta ) result( vol )
     use mod_geometry , only : num_dim
     implicit none
     
     class(t_element) , intent( inout ) :: this
     integer              , intent( in    ) :: elem_id 
     integer              , intent( in    ) :: dir                      ! xi(dir = 1) or eta( dir = 2) or zeta( dir = 3) 
     real( kind = r_kind) , intent( in    ) :: xi
     real( kind = r_kind) , intent( in    ) :: eta
     real( kind = r_kind) , intent( in ) , optional :: zeta             ! geometry can be 2D or 3D
     real( kind = r_kind ) :: vol , n(3)
       
     
     if( num_dim .eq. 2 ) then
         call this%local_derivative_of_shape_function( elem_id , xi , eta )
     else 
         if( .not. present( zeta ) ) then
             stop ' Not enough argument in get_derivative_shape_function for 3d element: in element module'
         end if
         
         call this%local_derivative_of_shape_function( elem_id , xi , eta , zeta )
     end if
     
     call this%get_jacobian_matrix( elem_id )                           !! find Jacobian matrix
      
      
     if( num_dim .eq. 2 ) then
     
         if( dir .eq. 1 ) then  !! fixed direction is xi
             vol = sqrt( this%jacobian_matrix( 2 , 1 ) ** 2 + this%jacobian_matrix( 2 , 2 ) ** 2 )                      
         else                   !! fixed direction is eta
             vol = sqrt( this%jacobian_matrix( 1 , 1 ) ** 2 + this%jacobian_matrix( 1 , 2 ) ** 2 )     
         end if
         
         return
     end if  
     
     
     !! 3D problem, face is 2D
     
     if( dir .eq. 1 ) then      !! fixed direction is xi  
        
         n(1) = this%jacobian_matrix( 2 , 2 ) * this%jacobian_matrix( 3 , 3 ) - &
              & this%jacobian_matrix( 2 , 3 ) * this%jacobian_matrix( 3 , 2 )
                
         n(2) = this%jacobian_matrix( 2 , 3 ) * this%jacobian_matrix( 3 , 1 ) - &
              & this%jacobian_matrix( 2 , 1 ) * this%jacobian_matrix( 3 , 3 )
                
         n(3) = this%jacobian_matrix( 2 , 1 ) * this%jacobian_matrix( 3 , 2 ) - &
              & this%jacobian_matrix( 2 , 2 ) * this%jacobian_matrix( 3 , 1 )       
                                                       
     elseif ( dir .eq. 2 ) then   ! eta is fixed: cross product of first and third rows of Jacobian matrix.
            
         n(1) = this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 3 , 3 ) - &
              & this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 3 , 2 )
                
         n(2) = this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 3 , 1 ) - &
              & this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 3 , 3 )
                
         n(3) = this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 3 , 2 ) - &
              & this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 3 , 1 )   
            
     else                         ! zeta is fixed: cross product of first and second rows of Jacobian matrix.
            
          n(1) = this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 2 , 3 ) - &
               & this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 2 , 2 )
                
          n(2) = this%jacobian_matrix( 1 , 3 ) * this%jacobian_matrix( 2 , 1 ) - &
               & this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 2 , 3 )
                
          n(3) = this%jacobian_matrix( 1 , 1 ) * this%jacobian_matrix( 2 , 2 ) - &
               & this%jacobian_matrix( 1 , 2 ) * this%jacobian_matrix( 2 , 1 )  
     end if
       

     
     vol = sqrt( n(1) ** 2 +  n(2) ** 2 + n(3) ** 2 )                  ! norm of normal vector
    

 
 end function volume_element_at_face
 
 
 
!!=======================================================================
!!
!!  subroutines related to Q8 element begin gere
!!
!!=======================================================================  
 
!!========================= Q8 constructor ===========================!
 
     subroutine  Q8_constructor( this )
          use mod_utils , only : Q8_node_parent
          implicit none
          class( t_Q8_element) , intent(inout) :: this
         
          this%num_main_nodes           =  4 
          this%num_node_in_element      =  8
          this%num_element_around_node  =  4  ! a node is surrondded by 4 quad elements 
          this%num_face                 =  4
         
          call  this%base_element_constructor()
          
          this%node_parent => Q8_node_parent
              
     
     end subroutine Q8_constructor
     
     
!! =========================== Q8 shape function ======================!

     subroutine  Q8_shape_function( this , elem_id , xi , eta , zeta )
     
          use mod_utils
          use mod_geometry , only : num_dim , element_matrix
          implicit none
  
        
          class( t_Q8_element) , intent( inout ) :: this
          integer              , intent( in    ) :: elem_id 
          real( kind = r_kind) , intent( in    ) :: xi
          real( kind = r_kind) , intent( in    ) :: eta
          real( kind = r_kind) , intent( in ) , optional :: zeta        ! Not importrant n 2d
         
          integer   :: i , k  , a 
          real(kind = r_kind) :: r(2) , term 
        
          r(1:2) = (/ xi , eta /)
        
          this%shape_vector    = 0.0_r_kind                             ! initial values
        
          do i = 1 , this%num_node_in_element                           ! find base part of shape functions      
             
              if( element_matrix( elem_id , i ) .eq. 0 ) cycle          ! check if node exist 
                                               
              this%shape_vector(i)  = 1.0_r_kind  
               
              do k = 1 , num_dim                                          
                 a = this%node_parent( k , i )                          ! (xi_i , eta_i)
               
                 if( a /= 0 ) then
                     term   = 0.50_r_kind * ( 1.0_r_kind + a * r(k) )               
                 else
                     term   =  1.0_r_kind - r(k) ** 2                        ! main term
                 end if 

                 this%shape_vector(i) = this%shape_vector(i) * term  
               
              end do ! loop over k    
          end do    ! loop over i
 
          ! add share of optional nodes to shape function of main nodes
          call this%add_optional_node_share( elem_id , .true. )
       
     end subroutine Q8_shape_function
  

!===================== Q8 derivative of shape function ================!
  
     
 subroutine Q8_local_derivative_shape_function( this , elem_id , xi , eta , zeta )
 
      use mod_geometry , only : num_dim , element_matrix
      implicit none
      
      class( t_Q8_element) , intent( inout ) :: this
      integer              , intent( in    ) :: elem_id 
      real( kind = r_kind) , intent( in    ) :: xi
      real( kind = r_kind) , intent( in    ) :: eta
      real( kind = r_kind) , intent( in ) , optional :: zeta            ! not importrant in 2d
     
      integer   :: i , k , n , a 
      real(kind = r_kind) :: r(2) , term , d_term
                                   
      r(1:2) = (/ xi , eta /)
       
      this%local_d_shape_matrix    = 0.0_r_kind                         ! initial values
                   
      ! find base part of shape functions 
      do i = 1 , this%num_node_in_element               
                             
         if( element_matrix( elem_id , i ) .eq. 0 ) cycle               ! check if node exist  
                             
         this%local_d_shape_matrix(:,i) = 1.0_r_kind
               
         do k = 1 , num_dim
               
            a = this%node_parent( k , i )                               ! (xi_i , eta_i )
         
            if( a /= 0 ) then
                term   = 0.50_r_kind * ( 1.0_r_kind + a * r(k) )        ! main term
                d_term = 0.50_r_kind * a                                ! its derivatve 
            else
                term   =  1.0_r_kind - r(k) ** 2                        ! main term
                d_term = -2.0_r_kind * r(k)                             ! its derivative
            end if 
            
            
            do n = 1 , num_dim  
               if( k .eq. n ) then                                      ! derivative term is multplied                                
                  this%local_d_shape_matrix(n,i) = this%local_d_shape_matrix(n,i) * d_term  
               else
                  this%local_d_shape_matrix(n,i) = this%local_d_shape_matrix(n,i) * term    ! main term is multiplied
               end if
            end do ! loop over n
             
              
         end do ! loop over k    
      end do    ! loop over i
 
      call this%add_optional_node_share( elem_id , .false. )       ! find modification part of derivative shape matrix 
    
     
 end subroutine Q8_local_derivative_shape_function
 
 

 
 
!!=======================================================================
!!
!!  subroutines related to brick element begin here
!!
!!=======================================================================  
 
!!========================= brick constructor ===========================!   
   
  subroutine  brick_constructor(this)
  
      use mod_utils , only : brick_node_parent   
      implicit none
      class( t_brick_element) , intent(inout) :: this
         
      this%num_main_nodes          =  8 
      this%num_node_in_element     =  20
      this%num_element_around_node =  8
      this%num_face                =  6
         
      call  this%base_element_constructor() 
      
      this%node_parent => brick_node_parent
          
  end subroutine brick_constructor
  
  
!!========================= brick shape function ======================!
  
  subroutine  brick_shape_function( this , elem_id , xi , eta , zeta )
     
      use mod_geometry , only : num_dim , element_matrix
      implicit none
        
      class( t_brick_element ) , intent( inout ) :: this
      integer              , intent( in    ) :: elem_id 
      real( kind = r_kind) , intent( in    ) :: xi
      real( kind = r_kind) , intent( in    ) :: eta
      real( kind = r_kind) , intent( in ) , optional :: zeta
         
      integer :: i , k  , a 
      real( kind = r_kind ) :: r(3) , term 
         
      r(1:3) = (/ xi , eta , zeta /)
        
      this%shape_vector    = 0.0_r_kind                                 ! initial values
        
      do i = 1 , this%num_node_in_element                               ! find base part of shape functions          
             
          if( element_matrix( elem_id , i ) .eq. 0 ) cycle              ! check if node exist 
                                               
          this%shape_vector(i)  = 1.0_r_kind  
               
          do k = 1 , num_dim                                          
               a = this%node_parent( k , i )                            ! (xi_i , eta_i , zeta_i )
               
               if( a /= 0 ) then
                   term   = 0.50_r_kind * ( 1.0_r_kind + a * r(k) )               
               else
                   term   =  1.0_r_kind - r(k) ** 2                     ! main term
               end if 

               this%shape_vector(i) = this%shape_vector(i) * term  
               
          end do ! loop over k    
      end do    ! loop over i
      
      call this%add_optional_node_share( elem_id , .true. )             ! find modification part of shape vector  
       
  end subroutine brick_shape_function
  
  
  
!! ========================= brick derivarive of shape function =======! 
    
  subroutine brick_local_derivative_shape_function( this , elem_id , xi , eta , zeta )
 
      use mod_geometry , only : num_dim , element_matrix
      implicit none
      
      class( t_brick_element) , intent( inout ) :: this
      integer               , intent( in ) :: elem_id 
      real( kind = r_kind ) , intent( in ) :: xi
      real( kind = r_kind ) , intent( in ) :: eta
      real( kind = r_kind ) , intent( in ) , optional :: zeta
     
      integer:: i , k , n , a 
      real( kind = r_kind ) :: r(3) , term , d_term
                                   
      r(1:3) = (/ xi , eta , zeta /)
       
      this%local_d_shape_matrix    = 0.0_r_kind                              ! initial values
                   
      do i = 1 , this%num_node_in_element                              ! find base part of shape functions             
                             
         if( element_matrix( elem_id , i ) .eq. 0 ) cycle               ! check if node exist  
                             
         this%local_d_shape_matrix(:,i) = 1.0_r_kind
             
         do k = 1 , num_dim
               
            a = this%node_parent( k , i )                               ! (xi_i , eta_i , zeta_i )
         
            if( a /= 0 ) then
                term   = 0.50_r_kind * ( 1.0_r_kind + a * r(k) )        ! main term
                d_term = 0.50_r_kind * a                                ! its derivatve 
            else
                term   =  1.0_r_kind - r(k) ** 2                        ! main term
                d_term = -2.0_r_kind * r(k)                             ! its derivative
            end if 
            
            
            do n = 1 , num_dim  
               if( k .eq. n ) then                                      ! derivative term is multplied                                
                  this%local_d_shape_matrix(n,i) = this%local_d_shape_matrix(n,i) * d_term  
               else
                  this%local_d_shape_matrix(n,i) = this%local_d_shape_matrix(n,i) * term    ! main term is multiplied
               end if
            end do ! loop over n
             
              
         end do ! loop over k    
      end do    ! loop over i
       
      call this%add_optional_node_share( elem_id , .false. )            ! find modification part of shape vector 
     
      
 end subroutine brick_local_derivative_shape_function
 
 
 

 
 
 end module mod_element
