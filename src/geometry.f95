!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!


!
!
!  Geometry module containing a geometry class with members: node, element, dof, connectivity matrix
!                     1) nodes  , 2) elements , 3) connectivty matrix
!                     
!  Also, shape function, derivative of shape function and Jacobian matrix for transformation and its determinant                   
!  are also calculated as member functions.
!
!                
!

!!  Element_type as listed in the input file: 
!!
!!        1 : brick element with 8 nodes                              (Elastic media)   (modelled by 3D20nodes )
!!        3 : arbitrary interface element with 8 node and more        (Elastic media)   (not modelled)
!!        4 : arbitrary brick with 8 node and more                    (Elastic media)   (modelled by 3D20Nodes )
!!        5 : quadrilaterial with 8 node for plane strain  condition  (Elastic media)   (modelled by 2D8Nodes  )
!!        6 : quadrilaterial with 8 node for plane stress  condition  (Elastic media)   (modelled by 2D8Nodes  )
!!        12: quadrilaterial with 8 node for acustic media (dam)      (Acustic media)   (modelled by 2D8Nodes  )
!!        14: arbitrary brick with 8 node and more for acustic media  (Acustic media)   (modelled by 3D20Nodes )
!!
!!
!!
!!   Remark. It is supposed that number of nodes is smaller than huge( integer(k =2 ) ) = 32767
!!
!!


!!  substructure comes here.
!!  numgroup is a logical array that says if dam exist, 
!!  possibilities
!!        1. dam   2. dam-reservior  3. dam-foundation  4. dam-reservior-foundation
!!  d.o.f is a logical variable.
!!  The shape function is the same for otherside boundary element. The shape function 
!!  is independent of chosen element.
!! only those things that belong to here comes here.
!! physics part includes: material all of them,
!! boundary condition all of them
!! Later boundary condition are used to define each boundary
!! Matrix element is defined to be of 2000 size and then it is checked it is more then it is reallocated 
!! by doubling it. 
!! Node, elements are geometry.
!! dof and substructure, material are physics, boundary conditions are physics
!! Analysis, node to report
!! forces are just in time domain so not be used here. are in different file


!!  Remark : degree of freedom is stored in matrix whose number of columns are num_dof + 1.
!!           1. The first column stores node classes in physics module.
!!
!!
!!
!!   dof matrix : 
!!     Column  1  =  class_id : (1 to 5 regions) defined in physics module
!!     Columns 2 : 1 + num_dim = translation dof
!!     column  1+ num_dof  = pressure dof (if exist ) 
!!===========================================================================================================

module  mod_geometry   !! geometry module

       use  mod_utils , only :  r_kind                                                     
       implicit none
       
       
       !! member variables
       
       integer , public  :: num_dim   = 0                               ! number of dimensions of model
       integer , public  :: num_node  = 0                               ! number of nodes 
       integer , public  :: num_element  = 0                            ! number of elements in FEM
       integer , public  :: num_degree_of_freedom                       ! DOF
       integer , public  :: max_num_node_in_element
       integer , public  :: num_element_around_node
       integer , public  :: max_num_connected_nodes
       
   
       
       ! input data: read-write
       real( kind = r_kind ) , allocatable , public :: node(:,:)        ! nodes coordinates
       integer , allocatable , public :: element_matrix(:,:)            ! nodes coordinates
       integer , allocatable , public :: degree_of_freedom(:,:)         ! dof
       
       integer , allocatable , public :: connectivity_matrix(:,:)       ! nearest nearby nodes to every node
       
       
      

!!--------------------------------------------------------------------!! 
!! module functions and subroutines
       
       public    :: initialize_node_matrix
       public    :: initialize_element_matrix
       public    :: form_connectivity_matrix
       public    :: clean_up_geometry
       
       public    :: show_geometry
       

 
    
 contains
     


      


     
        
!!=======================================================================
!!
!!  initialize node matrix
!!
!!=======================================================================
!! the dof matrix has an additional column, see structure module for description.
!! First column of dof matrix is class_id (see physics module).
!! Beside there are num_dof additional columns, so num_col = num_dof + 1

  subroutine initialize_node_matrix()
       
       implicit none
       integer :: state
    
       if( ( num_node < 1 ) .or. ( num_dim < 2 ) ) then
           stop ' invalid number of nodes and/or dimensions'
       end if
       
       if( num_degree_of_freedom < 2  ) then
           stop ' invalid number of dof'
       end if
         
       
       if( allocated( node ) ) then
           deallocate( node , stat = state  )
           if (state /= 0 ) stop 'Error: failed deallocating node array in geometry module.'
       end if
       
       if( allocated(  degree_of_freedom ) ) then
           deallocate( degree_of_freedom , stat = state )
           if (state /= 0 ) stop 'Error: failed deallocating degree_of_freedom matrix in geometry module.'
       end if
       
       
       allocate( node( num_node , num_dim ) , stat = state)  
       if ( state /= 0)  stop 'Error: failed in allocating node matrix in geometry module.'
       
        
       allocate( degree_of_freedom( num_node , 2 * num_degree_of_freedom + 1 )  , stat = state) 
       
       
       
       if ( state /= 0)   stop 'Error: failed in allocating dof matrix in geometry module.' 
  
       
       degree_of_freedom = 0
         
  end subroutine initialize_node_matrix
  
  
  
!!=======================================================================
!!
!!  initialize element matrix
!!
!!=======================================================================


  subroutine initialize_element_matrix()
       
       implicit none
       integer :: state
    
       if( (num_element < 1 ) .or. ( max_num_node_in_element < 3 )  ) then
           stop ' invalid number of elements and/or nodes in element'
       end if
       
       if( allocated( element_matrix ) ) then
            deallocate( element_matrix , stat = state  )
            if (state /= 0 ) stop 'Error: failed deallocating element matrix in geometry module.'
       end if
       
       
       allocate( element_matrix( num_element , max_num_node_in_element ) , stat = state)  
       if ( state /= 0)  stop 'Error: failed in allocating element matrix in geometry module.'
      
  end subroutine initialize_element_matrix  




!=======================================================================
!!
!!  form connectivity matrix:
!!
!!  This matrix contains index of nearest nodes to every active node.
!!  If a node has no dof, there is no need to get its connectivity matrix.
!!  In Q8 element, every node is surrounded by at most 20 other nodes, so number of columns is 21.
!!  In brick element, every node is surrounded by at most 80 nodes, so number of columns is 81.
!!
!!  This matrix is used to form two dictionaries in physics module.
!!  
!=======================================================================

  subroutine  form_connectivity_matrix()
     use mod_utils , only : elements_base_info
     implicit none
     integer :: i , j , k , r , n_nodes 
     integer , allocatable , dimension(:,:) :: minmax_element_nodes     ! min and max node ids in each element
     integer , allocatable , dimension(:)   :: an_index_set             ! a temporary index set to store node ids
     logical :: is_in_element
     
     
     !!--------------------------------------------------------------------!!

     !! step 0. assign memory to arrays.
     
     
     
     allocate( minmax_element_nodes( num_element , 2 ) , stat = i )
     if( i /= 0 ) stop 'Error: failed allocating minmax_element_node in geometry module.' 
     
     !! get max size for nearby nodes to a given node: the size is worse case situation
     
     if( num_dim .eq. 2 ) then
         n_nodes = elements_base_info(1,2) * elements_base_info(1,3)    ! num_node_in_element * num_element_around_node
     else
         n_nodes = elements_base_info(2,2) * elements_base_info(2,3)    ! num_node_in_element * num_element_around_node
     end if
     
     allocate( connectivity_matrix( num_node , 2 * n_nodes ) , stat = i )  !! factor 2 handles worse case meshing.
    
     if( i /= 0 ) stop 'Error: failed allocating connectivity_matrix in geometry module.' 
     
     !! if num_element is smaller than num_element_around_node then presumed num_nearby_nodes is larger than num_nodes 
     
     allocate( an_index_set( max( num_node , 2 * n_nodes )  )  , stat = i )      
     if( i /= 0 ) stop 'Failed in alloating an_index_set array in geometry module '
     
     
!!--------------------------------------------------------------------!!
     
     !! step 1. find minmax matrix: min and max node ids in each element. usefuld in forming connectivity matrix
     
     do i = 1 , num_element
        minmax_element_nodes(i,2) = maxval( element_matrix(i, : ) )        ! max_val is captured
        minmax_element_nodes(i,1) = minmax_element_nodes(i,2)              ! initial value of min_val
        
        do k = 1 , max_num_node_in_element
           if( element_matrix(i,k) .eq. 0 ) cycle                          ! among nonzero node indices
           minmax_element_nodes(i,1) = min( minmax_element_nodes(i,1) , element_matrix(i,k) )
           
        end do ! loop ober k
        
     end do ! loop over i
     
!!---------------------------------------------------------------------!!
     
     !! step 2. find connectivity matrix.
     
     max_num_connected_nodes = 0 ! initial value
     connectivity_matrix = 0    ! initial value  
     an_index_set  = 0          ! intial value
     
     do i = 1 , num_node
     
        if( degree_of_freedom( i , 1 ) .eq. 0 ) cycle                   ! node has no dof: does not matter to get nearby nodes
        
        n_nodes = 0                                                     !! number of actual nearby nodes
          
        do  j = 1 , num_element                                         !! step 2.1. find elements that contain i-th node                                
        
            if( ( i < minmax_element_nodes(j , 1 ) ) .or. &
                ( i > minmax_element_nodes(j , 2 ) ) ) cycle            ! fast way to get negative answer.                   
            
            
            is_in_element = .false.
            
            do k = 1 , size( element_matrix , 2 )                       ! check if node is in the element.
               if( element_matrix( j , k ) .eq. i ) then
                   is_in_element = .true.
                   exit
               end if
            end do ! over k
            
            
            if( .not. is_in_element ) cycle                             ! node not in element
               
            
            !! step 2.2. store all nodes in the element in an_index_set vector
            
            if( n_nodes + 1 > size( an_index_set ) ) then
                stop ' Number of nodes in nearby elements of a node is larger than presumed value: geometry module'
            end if
            an_index_set( n_nodes + 1 : n_nodes + size( element_matrix , 2 ) ) = element_matrix( j , : )    
            n_nodes = n_nodes +  size( element_matrix , 2 ) 
                
            
        end do ! over j -loop : all elements.
  
       
        
        
        !! step. 2.3. sort node indices in an_index_set vector assendingly 
        
        do k = 1 , n_nodes - 1                                             
            do j = k + 1 , n_nodes 
            
               if( an_index_set( k ) > an_index_set( j ) ) then
                   r  = an_index_set( k ) 
                   an_index_set( k ) = an_index_set( j )
                   an_index_set( j ) = r 
               end if 
               
            end do ! over j
         end do  ! over k
         
         
         !! step. 2.4.  set repeated elements to zero. 
         
         do k = 1 , n_nodes - 1    
         
            if( an_index_set( k )  .eq. 0 ) cycle
            
            do j = k + 1 ,  n_nodes     
               if( an_index_set( j )  .eq. 0 ) cycle  
                                             
               if( an_index_set( k )  .eq.  an_index_set( j ) ) then
                   an_index_set( j ) = 0 
               end if
               
            end do ! over j
         end do  ! over k
         
         !! step 2.5  put non-zero elements in connectivity matrix
         r = 0
         do  k = 1 , n_nodes 
             if( an_index_set( k )  .eq. 0 ) cycle
              r = r + 1
              if( r > size( connectivity_matrix , 2 ) ) then
                  stop 'Number of nodes in nearby elements exceeds dimension of connectivity matrix: geometry module'
              end if
              connectivity_matrix( i , r ) = an_index_set( k )
         end do ! over k
         
         max_num_connected_nodes = max( max_num_connected_nodes  , r )
         
         
         !! step. 2.6. set n_node and an_index_set  zero for next iteration
         
         an_index_set( 1:n_nodes ) = 0

     
     end do ! over i-loop : num_nodes 
     
    
     !! step 2.7. free unused memory
     
     deallocate( an_index_set , stat = i )
     if( i /= 0 ) stop ' Error in deallocating an_index_set array in geometry module'
     
     deallocate( minmax_element_nodes , stat = i )
     if( i /= 0 ) stop ' Error in deallocating minmax_element_nodes array in geometry module'
     
  
  end subroutine form_connectivity_matrix

  

!=======================================================================
!
!  deallocate memory for geometry object
!
!=======================================================================

  subroutine   clean_up_geometry()
              
              implicit none
              integer :: state = 0 , state1 = 0
              
              if( allocated(  node ) ) then               
                  deallocate( node , stat = state  )
                  state1 = state
              end if 
                    
              
              if( allocated(  element_matrix ) ) then               
                  deallocate( element_matrix , stat = state  )
                  state1 = state1 + state
              end if 
            
              
              if( allocated( degree_of_freedom ) ) then
                  deallocate( degree_of_freedom , stat = state )
                  state1 =  state1 + state 
              end if  
              
              if( allocated(  connectivity_matrix ) ) then               
                  deallocate( connectivity_matrix  , stat = state  )
                  state1 = state
              end if  
              
              
             if (state1 /= 0 ) stop 'Error: failed deallocating memory in geometry module.'
              
  end subroutine clean_up_geometry     
   
   
   
   
   
!!======================================================================
!!
!!  show geometric information
!!
!!======================================================================

  subroutine  show_geometry()
  
          implicit none
          integer :: i
          
          
          
          write(*,*) ' num_nodes = ' , num_node             
          write(*,*) ' num_dim = '   , num_dim
          write(*,*) ' num_element = '   , num_element
          write(*,*) ' num_nodes in element = '   , max_num_node_in_element

          
         
!          do i = 1 , num_node
!             write(*,*)  i , node(i,:)
!          end do
             
!          write(*,*) ' num_elements = ' , num_element
!          do i = 1 , num_element
!             write(*,*)  i , element_matrix(i,:)
!          end do
          
          write(*,*) ' num_dof+1 = ' , size( degree_of_freedom , 2 )
          write(*,*) ' class_id , dof(1:num_dof) ' 
          do i = 1 , num_node
             write(*,*)  i , degree_of_freedom(i,:)
          end do
          

  end subroutine show_geometry
      
          

end module mod_geometry
