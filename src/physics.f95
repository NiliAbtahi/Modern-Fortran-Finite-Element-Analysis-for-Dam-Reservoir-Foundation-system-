!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!


!!
!!  This module contains physical information of the structure.
!!  The name of variables are enough clear.
!!
!!  Different regions are specified based on following rule.
!!  
!!  Region:
!!        1. dam (region_id = 1 ), foundation (region_id = 2 ) , reservior (region_id = 3 )
!!
!!  Boundaries: 
!!            1. dam_reservior boundary:
!!               1.1. in dam region the ktype of such elements is 2, so       (region = 1 , ktype = 2 )
!!               1.2. in reservior region the ktype of such elements is 3, so (region = 3 , ktype = 3 )
!!
!!            2. reservor-foundation boundary:
!!               2.1. in reservor regon the ktype of such elements is 2, so   (region = 3 , ktype = 2 )
!!               2.2. There is no interaction between foundation-reservior,
!!                    so such elements are not considered here.
!!
!!            3. foundation-infinity boundary:
!!               3.1. in foundation the ktype of such elements is 3, so       (region = 2 , ktype = 3 ) 
!!
!!
!!            4. reservior-infinite boundary:
!!               4.1. in reservior the ktype of such elements is 1, so        (region = 3 , ktype = 1 )
!!
!!---------------------------------------------------------------------!!
!!
!!  Classify nodes:
!!
!!  class 0. nodes with no dof.
!!
!!  class 1. dam_reservior boundary : such nodes do have both translational and pressure dof.
!!
!!  class 2. foundation_infinity boundary: the nodes at truncation boundary of foundation.
!!          
!!  class 3. reservior-infinity boundary: the nodes at truncation boundary of reservior.
!!
!!  class 4. elastic_nodes: all nodes in dam and foundation that are not in class  1 , 2 ( or, bulk nodes )
!!          
!!  class 5. acoustic nodes: all nodes in reservior that are not in class 1 and 3 ( or, bulk nodes )
!!
!!
!!---------------------------------------------------------------------!
!!
!!  
!!  num_state_variables(1)  : num state function for elastic nodes in x-direction          : class_id = 4 
!!  num_state_variables(2)  : num state function for elastic nodes in y-direction          : class_id = 4 
!!  num_state_variables(3)  : num state function for elastic nodes in z-direction          : class_id = 4 
!!  num_state_variables(4)  : num state function for acoustic nodes in bulk reservior      : class_id = 5
!!  num_state_variables(5)  : num state function for nodes in reservior-infinite  boundary : class_id = 3  
!!  num_state_variables(6)  : num state function for nodes in foundation-infinite boundary : class_id = 2 
!!
!! 
!!         
!!===========================================================================================!!

 

    module  mod_physics
      
        use mod_utils , only :  r_kind , t_material , t_boundary_load, t_model_description, t_drm_layer
        use mod_utils , only :  t_enforced_disp , t_boundary_spring , t_nodal_load, t_load_case
        implicit none
      
      

!!--------------------------------------------------------------------!!      
!! variables related to  basic physical quantities
   
   
        integer , public  :: num_substructure                           ! NGROUP
        integer , public  :: num_material                               ! sum of NMAT
        integer , public  :: num_boundary_load                          ! sum of NULD
        integer , public  :: num_load_cases = 0                         ! NLOADC
        integer , public  :: newmark_accel_type = 0                     ! NACCEL
        integer , public  :: num_static_step = 0                        ! NSTEP2
        integer , public  :: analysis_kind = 0                          ! NDYN
        logical , public  :: is_time_domian   
        integer , public  :: earthquake_direction                         ! NGDIR = 1 (x_axis) , 2 (y-axis) , 3 (z-axis)
      
      
        real( kind = r_kind ) , public  :: gravitational_constant         ! GRV
        real( kind = r_kind ) , public  :: alpha_reflection_coefficient   !! alpha coeffient in reservior-foundation.
      
        type( t_material      ) , allocatable , public :: material(:)               ! all MAT 
        type( t_boundary_load ) , allocatable , public :: boundary_load(:)          ! all ULD 
        
        type( t_enforced_disp   ) , public :: enforced_displacement                 ! enforced displacement data 
        type( t_boundary_spring ) , public :: boundary_springs                      ! foundation boundary springs 
        type( t_nodal_load      ) , public :: nodal_loads                           ! nodal nodes
        type( t_load_case  )      , public :: load_cases                            ! all load conditions
      
        integer                 , allocatable , public  :: element_property(:,:)
        integer                 , allocatable , public  :: translational_variable(:)
        real( kind = r_kind )   , allocatable , public  :: rigidity_matrix(:,:)     
      

!!--------------------------------------------------------------------!!
!! variables related to physical model problem


    
        type( t_model_description ) , public :: the_model_description
 
    
   
        integer , public  :: num_state_variables(6)                         ! number of state functions for different node classes 
    
        integer , public  :: element_type(3)
     
        logical , public  :: is_dam_present
        logical , public  :: is_reservior_present
        logical , public  :: is_foundation_present
     
                                                                        ! check if foundation mass-model is complete.
        logical , public :: is_complete_mass_model                      ! true if jinertia(igroup = foundation) = 1 
        integer , public :: dam_jinertia = -1 ;
      
       
      
!!--------------------------------------------------------------------!!
!! variables related to node classification and variable space, dictionary

    integer , protected  :: total_num_variables
    integer , protected  :: num_nodes_in_classes(5)
    integer , public     :: num_max_variables_in_eqns 
    
    integer , allocatable  , dimension(:) , protected :: dictionary                ! first  dictionary
    integer , allocatable  , dimension(:) , protected :: dependency                ! second dictionary

    
!!--------------------------------------------------------------------!!
!! variables for time domain analysis

    real( kind = r_kind ) , allocatable  , dimension(:)   , public  ::  dead_weight_force
    real( kind = r_kind ) , allocatable  , dimension(:)   , public  ::  surface_pressure_force
    real( kind = r_kind ) , allocatable  , dimension(:,:) , public  ::  tension
    real( kind = r_kind ) , allocatable  , dimension(:,:) , public  ::  compression
    

!!--------------------------------------------------------------------!!
!! variables for DRM analysis
    integer, public :: drm_layer_ndyn   = 0     ! 0 for no DRM analysis, 1 for free field data , 2 for noFreeField:Direct Deconvlution data
    type( t_drm_layer ) , public :: drm_layer_info
    
    logical , public  :: is_TR_dam
    real( kind = r_kind ) , allocatable, public :: drm_imported_data(:,:)  !! reading a block of data in DRM analysis
    real( kind = r_kind ) , allocatable, public :: u_uddot_njtf_drm(:,:) !! storing time-domain data of u and uddot njtf
    
!!--------------------------------------------------------------------!!
!! coefficient matrices in global set of equations
    
     
    real( kind = r_kind ) , allocatable  , dimension(:) , target  :: mass          ! a vector to store global mass matrix
    real( kind = r_kind ) , allocatable  , dimension(:) , target  :: stiffness     ! a vector to store global stiffness matrix
    real( kind = r_kind ) , allocatable  , dimension(:) , target  :: damping       ! a vector to store global stiffness matrix
    real( kind = r_kind ) , allocatable  , dimension(:,:)         :: force         ! a matrix to store force
    
    real( kind = r_kind ) , dimension(:) , pointer  :: destination_matrix       ! points to either of mass, stifness, or damping
    
 
!!--------------------------------------------------------------------!!
!! subroutines related to memory assignment
  
      
      public    :: initialize_physics
      public    :: clean_up_physics
      
!!--------------------------------------------------------------------!!
      
      public    :: prepare_physical_data
      public    :: get_variable_location
      private   :: classify_nodes
      private   :: define_variables
      private   :: set_dictionary_of_variables
      private   :: initialize_global_matrices
      
      public    :: add_to_global_matrices
      public    :: add_to_resultant_force
      public    :: find_DRM_internal_face

  contains
  


!!=====================================================================!!
!! initialize element property matrix. It has 4 + max(NULD) columns.

!! column 1 = region_id (1 =dam , 2 = foundation , 3 = reservior )   
!! column 2 = material_id of element (refers to materia(i) class )
!! column 3 = num_Gauss_point for each element
!! column 4 to 3 + max(NULD) refers to elements of boundary_load(id) class
!!
!!======================================================================



 subroutine  initialize_physics( n_bc_max )  ! max_nuld
  
     use  mod_geometry , only : num_element , num_dim , num_degree_of_freedom
     implicit none
     integer , intent( inout ) :: n_bc_max
     integer :: state = 0 , k
     

!!--------------------- substructure info  -----------------------------!!     

      is_dam_present        = .false.                                 
      is_reservior_present  = .false.
      is_foundation_present = .false.
         
     
      ! which substructure does exist ?
      if( num_substructure .eq. 1 ) then                               ! dam-only structure
     
           is_dam_present        = .true.
         
      else if( num_substructure .eq. 2 ) then
            is_dam_present        = .true.                              ! dam is always there  
                
         if( num_dim  .eq. num_degree_of_freedom ) then                 ! there is no reservior: dam-foundation structure
             is_foundation_present = .true.
         else
             is_reservior_present  = .true.                             ! dam-reservior structure
         end if
     
      else                                                             ! dam-foundation-reservior structure
         is_dam_present        = .true.                                 
         is_reservior_present  = .true.
         is_foundation_present = .true.
      end if
      
      is_complete_mass_model = .false. 
      
    
!!-------------------- set physical model ----------------------------!!

    num_state_variables = 1    
      
    !! HW bc model at reservior truncation boundary : 
    !! num_valence: N (phi_i) + 1 (P) + 1 (P_i) + M (Phi_i) + 1 (psi) = N+M+3
    !! if M = 0 then there is no psi, so  N(phi_i) + 1 (P) + 1 (P_i) = N+2
      
        if( the_model_description%reservior_inf_boundary_model_id .eq. 1 ) then
            k = 2 + the_model_description%order_of_reservior_higher_order_bc(1) ;   !phi_i {N} + P + P_i   
          
            if( the_model_description%order_of_reservior_higher_order_bc(2) > 0 ) then  ! if M > 0 add M (phi_j) + 1(psi)
                k = k + 1 + the_model_description%order_of_reservior_higher_order_bc(2) ;
            end if  
                    
            num_state_variables(5) = k 
          
        elseif(  the_model_description%reservior_inf_boundary_model_id .eq. 2 ) then  !! GN bc model
            k = 1 + the_model_description%order_of_reservior_higher_order_bc(1)  &  !phi_i {N} + P + P_i + phi_j(M-1) = N+M-1+2 = N+M+ 1
            & + the_model_description%order_of_reservior_higher_order_bc(2)     ! as phi_{N+M} = 0, one lesser than HW
                          
            num_state_variables(5) = k 
        end if
      
     
        !! HW model for foundation truncation boundary
        !! num_valence at each direction: N(phi_j) + 1 ( phi_0 = u ) + M (phi_{N+j} ) + 1 ( phi_x: continuity)
     
        if( the_model_description%foundation_inf_boundary_model_id .eq. 1 ) then
             k = 1 + the_model_description%order_of_foundation_higher_order_bc(1) ;   !\{ phi_i {N} + u_x \}   
          
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then  ! if M > 0 add M (phi_j) + 1(psi)
                k = k + 1 + the_model_description%order_of_foundation_higher_order_bc(2) ;
            end if  
                    
            num_state_variables(6) = k 
        end if
      
      
        if( num_dim  .eq. 2 ) then
            num_state_variables(3) = 0                                    ! elastic nodes z-direction
        end if
          
           
        if( .not. is_reservior_present ) then
            num_state_variables(4) = 0                                    ! acoustic nodes 
            num_state_variables(5) = 0                                    ! reservor-inf nodes
        end if
          
        if ( .not. is_foundation_present ) then
            num_state_variables(6) = 0                                  ! foundation-inf nodes
        end if    
       

!!-------------------- check input validity ---------------------------!!     

        if( ( num_substructure < 1 ) .or. (num_substructure > 3 ) ) then
            stop 'Invalid number of substructures (NGROUP) '
        end if
     
        if( num_material < 1  ) then
            stop 'Error: invalid num material(NMAT) '
        end if

!!------------------------ allocate memory -----------------------------!! 


        state = 0 ;
        if( allocated(  material ) )       deallocate( material , stat = state )
        if ( state /= 0)  stop 'Error: failed in deallocating element_property in physics module.' 
        if( allocated(  boundary_load ) )  deallocate( boundary_load , stat = state )
        if ( state /= 0)  stop 'Error: failed in deallocating element_property in physics module.' 
        if( allocated(  element_property)) deallocate( element_property , stat = state )
        if ( state /= 0)  stop 'Error: failed in deallocating element_property in physics module.' 
        
    
     
        allocate( element_property( num_element , 3 + n_bc_max ) , stat = state )
        if ( state /= 0)  stop 'Error: failed in allocating element_property in physics module.' 
     
     
        allocate( material( num_material ) , stat = state )
        if ( state /= 0)  stop 'Error: failed in allocating material array in physics module.'    
     
     
        if( n_bc_max > 0 ) then
            allocate( boundary_load( num_boundary_load ) , stat = state )
            if ( state /= 0)  stop 'Error: failed in allocating loads on boundary array in physics module.' 
        end if
     
        ! dam is always present so rigidity_matrix is required in any problem
        if( num_dim .eq. 2 ) then
            allocate( rigidity_matrix( 3 , 3 ) , stat = state )
        else
            allocate( rigidity_matrix( 6 , 6 ) , stat = state )
        end if
        if ( state /= 0)  stop 'Error: failed in allocating rigidity_matrix in physics module.' 
     
     
   
        element_property = 0
        element_type = 0
      
      
    end subroutine  initialize_physics
  
 
!!====================================================================!!
!!
!!   prepare physical data for fem analysis in 5 steps:
!!
!!   1. classification of nodes into 6 classes.
!!   2. forming variable space for the problem
!!   3. making suitable dictionaries for non-zero elements of coefficient matrices
!!   5. assigning suitable memory for coefficient matrices.
!! 
!!====================================================================!!
   
   subroutine  prepare_physical_data()

      implicit none
      
      call  classify_nodes()
      
      call  define_variables()
      
      call  set_dictionary_of_variables()
      
      call  initialize_global_matrices()
     
   
   end subroutine prepare_physical_data 
  
  
  
!!====================================================================!!
!!
!!   classify nodes:
!!
!!      class 0. nodes with no dof.
!! 
!!      class 1. dam_reservior boundary : such nodes do have both translational and pressure dof.
!!
!!      class 2. foundation_infinity boundary: nodes on face in elements where  ktype = 3 and region =2
!!
!!      class 3. reservior-infinity boundary: nodes on face in elements where  ktype = 1 and region =3
!!
!!      class 4. elastic_nodes: all nodes in dam and foundation except those of class 2 ( bulk nodes )
!!
!!      class 5. acoustic nodes: all nodes in reservior except those of class 3 ( bulk nodes )
!! 
!!
!!====================================================================!!
    
  subroutine  classify_nodes()
  
     use  mod_geometry , only :  num_node , num_dim , element_matrix , degree_of_freedom
     use  mod_utils    , only :  get_nodes_on_element_face , Q8_node_parent , brick_node_parent
     implicit none
     
     logical :: is_dof(2)                                               ! [ translation_dof , pressure_dof ]
     integer :: i ,  j , k
     integer :: n_face , node_id , ktype,  num_boundary_node , region_id
     
     integer :: node_index(10)                                          ! 3 for Q8 face and 8 for brick face
     
     
     num_nodes_in_classes = 0                                           ! default value is zero
     
     
!!--------------------------------------------------------------------!!     
!! classify nodes based on dof matrix 
     
     do  i = 1 , num_node 
     
         is_dof  = .false.                                              ! default value for different kind of dof
    
         do k = 1 , num_dim
            if( degree_of_freedom( i , 2 * k ) .eq. 1 ) then            ! there is translational dof for node i
                is_dof(1) = .true. 
                exit
            end if
         end do
         
         if( is_reservior_present ) then
             is_dof(2) = degree_of_freedom( i , 2 * num_dim + 2 ) .eq. 1  ! there is pressure dof for node i
         end if
         
         
         if( is_dof(1) .and. is_dof(2) ) then                           ! node is in dam-reservior boundary, class_id = 1 ,
         
             degree_of_freedom( i , 1 )   =  1                           
             num_nodes_in_classes(1) =  num_nodes_in_classes(1) + 1     ! class id = 1
             
         else if( is_dof(1)  .and. ( .not. is_dof(2) ) )  then          ! node is in elastic media: both bulk and boundary
         
             degree_of_freedom( i , 1 )   =  4                          
             num_nodes_in_classes(4) =  num_nodes_in_classes(4) + 1     ! class id = 4
             
         else if( (.not. is_dof(1) ) .and.  is_dof(2)  )  then          ! node is in acoustic media: both bulk and boundary
             
             degree_of_freedom( i , 1 )   =  5                          
             num_nodes_in_classes(5) =  num_nodes_in_classes(5) + 1     ! class id = 5
             
         end if
     end do
     
!!---------------------------------------------------------------------!!  
!!   
!! use property matrix to exclude boundary nodes in semi-infinite regions. 
!!
!! for nodes in foundation-infinity boundary: class_id = 2 and region_id = 2 so class_id = region_id
!! for nodes in reservior-infinity boundary : class_id = 3 and region_id = 3 so class_id = region_id
!! The equality is used to simplify following code.
!!
!!--------------------------------------------------------------------!!
                      
     if( size( element_property  , 2 ) .eq. 3 ) return                  ! there is no boundary condition.
     
     
     do i = 1 , size( element_property , 1 )                            ! loop over all elements
     
        if( element_property( i , 1 ) .eq.  1 ) cycle                   ! if region_id = 1, nodes lie in dam
        
        region_id = element_property( i , 1 )       
        
        do k = 1 , size( element_property  , 2 ) - 3                    ! loop over all boundary conditions
         
           if( element_property( i , 3 + k  ) .eq. 0 ) cycle                    ! there is no boundary condition
           
           ktype = boundary_load( element_property( i , 3 + k  ) )%k_type       ! get k_type
           
           if( ( region_id .eq.  2 ) .and. ( ktype /=  3 ) ) cycle              ! not in foundation-infinity boundary  
           if( ( region_id .eq.  3 ) .and. ( ktype /=  1 ) ) cycle              ! not in reservior-infinity boundary
           
           n_face = boundary_load( element_property( i , 3 + k  ) )%face_id     ! get boundary face
           
           if( n_face .eq. 0 ) cycle                                            ! not reasonable
           
           if( num_dim .eq. 2 ) then
               num_boundary_node  =  get_nodes_on_element_face( n_face , Q8_node_parent , node_index )    ! num of boundary nodes
           else
               num_boundary_node  =  get_nodes_on_element_face( n_face , brick_node_parent , node_index )  ! num of boundary nodes
           end if
           
           
           do j = 1 , num_boundary_node                                 ! loop over all non-zero nodes on n_face
           
              node_id = element_matrix( i , node_index(j) )
              
              if(  node_id  .eq. 0 ) cycle                              ! node is absent in element definition
              
              if( degree_of_freedom( node_id , 1 ) .eq. 0 ) cycle       ! node has no dof
              
              if( degree_of_freedom( node_id , 1 )  .eq. region_id ) cycle       ! node has already been classified
              
              degree_of_freedom( node_id , 1 )  = region_id                      ! class_id = region_id
              
              num_nodes_in_classes( region_id + 2 ) =  num_nodes_in_classes( region_id + 2 ) - 1   ! reduce node share from bulk
              
              num_nodes_in_classes( region_id ) =  num_nodes_in_classes(region_id ) + 1            ! add node share to boundary class
               
           end do ! over j
  
        end do ! over k
        
     end do ! over i
     
     
     
     
  end subroutine classify_nodes


  
!!====================================================================!!
!!
!!  Here the set of all variables are defined for the analysis.
!!  The dof matrix is adopted to variables indices
!!
!!  num_state_variables(1)  : num state function for elastic nodes in x-direction          : class_id = 4 
!!  num_state_variables(2)  : num state function for elastic nodes in y-direction          : class_id = 4 
!!  num_state_variables(3)  : num state function for elastic nodes in z-direction          : class_id = 4 
!!  num_state_variables(4)  : num state function for acoustic nodes in bulk reservior      : class_id = 5
!!  num_state_variables(5)  : num state function for nodes in reservior-infinite  boundary : class_id = 3  
!!  num_state_variables(6)  : num state function for nodes in foundation-infinite boundary : class_id = 2 
!!
!!====================================================================!!
  
  subroutine define_variables()
  
    use mod_geometry, only : degree_of_freedom , num_dim , num_node 
    implicit none
    
    integer :: i ,  k 
    
    
    total_num_variables = 0                                             ! initial value
    
    do i  = 1 , num_node
    
       if( degree_of_freedom(i , 1 ) .eq. 0 ) cycle                     ! node has no dof
       
       do k = 1 , num_dim                                               ! loop over translational dof
       
          if( degree_of_freedom(i , 2 * k ) .eq. 0 ) cycle              ! node has no translational dof in k-th direction
          
          degree_of_freedom(i , 2 * k     ) = total_num_variables + 1  
          
          if( degree_of_freedom( i , 1 ) .eq. 2 ) then                   ! node is in foundation-inf boundary
          
              degree_of_freedom( i , 2 * k + 1 ) = total_num_variables + num_state_variables(6)
          else
              degree_of_freedom( i , 2 * k + 1 ) = total_num_variables + num_state_variables(k)  ! node is in bulk
          end if
          
          total_num_variables =  degree_of_freedom( i , 2 * k + 1 )       ! update total number of variables
       
       end do ! over k
       
       
       if( .not. is_reservior_present ) cycle                           ! there is no pressure dof
       
       k = 2 * num_dim + 2                                              ! location where pressure dof exists
       
       if( degree_of_freedom(i , k ) .eq. 0 ) cycle                     ! node has no pressure dof 
       
       degree_of_freedom( i , k ) = total_num_variables + 1
       
       if( degree_of_freedom( i , 1 ) .eq. 3 ) then                     ! node is in reservior-inf boundary
          
           degree_of_freedom( i , k + 1 ) = total_num_variables + num_state_variables(5)
       else
           degree_of_freedom( i , k + 1 ) = total_num_variables + num_state_variables(4)  ! node is in accoustic bulk
       end if
          
          total_num_variables =   degree_of_freedom( i ,  k + 1 )         ! update total number of variables
       
    end do ! over i
    
   
     
  
  end subroutine define_variables
  
 
 
!!====================================================================!!
!!
!!  set two dictionaries for variables. They are called dictionary and dependency vectors.
!!
!!====================================================================!!
  
  subroutine  set_dictionary_of_variables()
     
     use mod_utils , only : r_kind , elements_base_info
     use mod_geometry
     implicit  none
     
     integer :: i , j , k , r,  num_dependent_variable , id_beg, id_end , tot_id , n_nodes
     
     
     !! step 0. call form connectivity matrix from geometry module
     
     call  form_connectivity_matrix()
     
     
     !! step 1. form dictionary of variables
     
     if( allocated(  dictionary ) ) then
         deallocate( dictionary , stat = i )
         if( i /= 0 ) stop ' Error in deallocating dictionary array in physics module'
     end if
     
     allocate( dictionary( total_num_variables + 1 ) , stat = i )
     if( i /= 0 ) stop ' Error in allocating dictionary array in physics module'
     
     
     num_max_variables_in_eqns = 0
     dictionary = 0
     
     
     
     do i = 1 , num_node
        
        if( degree_of_freedom(i , 1 ) .eq. 0 ) cycle                    !! the node has no dof
        
        num_dependent_variable = 0 
        
        do  j = 1 , max_num_connected_nodes                             !! evaluated in geometry module
        
            if( connectivity_matrix( i , j ) .eq. 0 ) exit              !! remaining elements in this row is zero
            
            !! step 1.1. get all state variables describing j-th node  in the current row of matrix
            
            n_nodes = connectivity_matrix( i , j )                      !! used as temporary variable
            
            if( degree_of_freedom( n_nodes , 1 ) .eq. 0 ) cycle         ! nearby node has no dof
            
            do k = 1 , num_degree_of_freedom
               if( degree_of_freedom( n_nodes , 2 * k ) .eq. 0 ) cycle  ! node has no dof for k-th variable
            
            num_dependent_variable = num_dependent_variable +  &
         &  degree_of_freedom( n_nodes , 2 * k + 1 ) - degree_of_freedom( n_nodes , 2 * k ) + 1
                    
            end do ! over k
            
        end do ! over j
        
        !! step 1.2. every state variable of node i is coupled with num_dependent_variable other variables
        
        num_max_variables_in_eqns = max( num_max_variables_in_eqns , num_dependent_variable  ) 
        
        !!  update dictionary
         
        do j = 1 , num_degree_of_freedom
           if( degree_of_freedom( i , 2 * j ) .eq. 0 ) cycle            ! node has no dof for j-th variable
            
           do k = degree_of_freedom( i , 2 * j ) , degree_of_freedom( i , 2 * j + 1 ) ! loop over all variables of node i
              dictionary( k + 1 ) = dictionary( k ) + num_dependent_variable  
           end do ! over r
        end do ! over j
        
     
     end do ! over nodes
     
     
    
!!!--------------------------------------------------------------------!!
     
     !! step 2. form second dictionary: dependency vector
     

     i = dictionary( total_num_variables + 1 )                          !! size of dependency vector
     
  
     if( allocated(  dependency ) ) then
         deallocate( dependency , stat = k )
         if( k /= 0 ) stop ' Error in deallocating dependency array in physics module'
     end if
     
     allocate( dependency( i ) , stat = k )
     if( k /= 0 ) stop ' Error in allocating dependency array in physics module'
     
     
     dependency = 0
     
     tot_id = 0
     
     do i = 1 , num_node
        
        if( degree_of_freedom(i , 1 ) .eq. 0 ) cycle                    !! the node has no dof
        
        id_beg = tot_id                                                 ! begin index for i-th node state variable
             
        do  j = 1 , max_num_connected_nodes  
        
            if( connectivity_matrix( i , j ) .eq. 0 ) exit              !! remaining elements in this row is zero
            
            !! step 2.1. get all state variables describing j-th node  in the current row of matrix
            
            n_nodes = connectivity_matrix( i , j )                      !! used as temporary variable
            
            if( degree_of_freedom( n_nodes , 1 ) .eq. 0 ) cycle         ! nearby node has no dof
            
            do k = 1 , num_degree_of_freedom
               if( degree_of_freedom( n_nodes , 2 * k ) .eq. 0 ) cycle  ! node has no dof for k-th variable
            
               do r = degree_of_freedom( n_nodes , 2 * k  ) , degree_of_freedom( n_nodes , 2 * k + 1 ) 
                  tot_id = tot_id + 1 
                  dependency( tot_id ) = r 
                
               end do 
               
            end do ! over k
            
        end do ! over j
        
        id_end = tot_id                                                 ! end index for i-th node state variable
        
        
        
        !! step 2.2. get number of state variables describing i-th node
        n_nodes = 0 
         
        do j = 1 , num_degree_of_freedom
           if( degree_of_freedom( i , 2 * j ) .eq. 0 ) cycle            ! node has no dof for j-th variable
            
           do k = degree_of_freedom( i , 2 * j ) , degree_of_freedom( i , 2 * j + 1 ) ! loop over all variables of node i
              n_nodes = n_nodes + 1 
           end do ! over r
        end do ! over j
        
        
        !! n_node is >= 1 as the node has some dof
        
        if( n_nodes .eq. 1 ) cycle                                       !! its dependency has already set.
        
        
        !! step 2.3. replicate the set of variables obtained for all other state variables of i-th node.
        
        r = id_end - id_beg                                  ! number of variables coupled to each state var of i-th node
        
        do j = 1 , n_nodes - 1
           dependency( tot_id + 1 : tot_id + r ) = dependency( id_beg + 1 : id_end ) 
           tot_id = tot_id + r 
        end do ! over j
        
        
     end do ! over i 
     
     !! deallocate connectivity matrix
     
     deallocate( connectivity_matrix , stat = k )
     if( k /= 0 ) stop ' Error in deallocating connectivity_matrix in physics module'
     
   
    
!      write(*,*) 'dof '
!      do i = 1 , size( degree_of_freedom , 1 )
!         write(*,*) i , degree_of_freedom(i, :)
!      end do
      
!      write(*,*) 'element matrix '
!      write(*,*) element_matrix(120 , :)
      
!      do i = 1 , size( element_matrix , 2 )
!         n_nodes = element_matrix(120 , i )
         
!         do k = 1 , num_degree_of_freedom
!            if( degree_of_freedom(n_nodes , 2 * k ) .eq. 0 ) cycle
!            write(*,*) '********************************************vars for node = ' , n_nodes
            
!            do  j = degree_of_freedom(n_nodes , 2 * k ) , degree_of_freedom(n_nodes , 2 * k + 1 )
            
!                id_beg = dictionary( j ) + 1
!                id_end = dictionary( j + 1 )
!                write(*,*) '**************vars for var id = ' , j
!                write(*,*) j , dependency(id_beg : id_end )
!            end do 
            
            
!         end do
!      end do

   
     
  end subroutine set_dictionary_of_variables
  
  
  
 
!!====================================================================!!
!!
!! assign memory to global matrices
!!
!!====================================================================!!

subroutine  initialize_global_matrices()
     
    use mod_geometry , only : num_dim
    implicit none
     
    integer :: state 
     
    if( allocated( mass ) ) then
        deallocate( mass , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating mass vector in physics module.' 
    end if
     
    if( allocated( stiffness ) ) then
        deallocate( stiffness , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating stiffness vector in physics module.' 
    end if
     
    if( allocated( damping ) ) then
         deallocate( damping , stat = state )
         if( state /= 0 ) stop 'Error: failed in deallocating damping vector in physics module.' 
    end if
     
    if( allocated( force ) ) then
        deallocate( force , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating force vector in physics module.' 
    end if
     
    allocate( mass( size( dependency ) ) , stat = state )
    if( state /= 0 ) stop 'Error: failed in allocating mass vector in physics module.' 
     
    allocate( stiffness( size( dependency ) ) , stat = state )
    if( state /= 0 ) stop 'Error: failed in allocating stiffness vector in physics module.' 
     
    allocate( damping( size( dependency ) ) , stat = state )
    if( state /= 0 ) stop 'Error: failed in allocating damping vector in physics module.' 
     
    allocate( force( total_num_variables , num_dim ) , stat = state )
    if( state /= 0 ) stop 'Error: failed in allocating force vector in physics module.'
     
    mass      = 0.0_r_kind                                             ! set initial value
    stiffness = 0.0_r_kind                                             ! set initial value
    damping   = 0.0_r_kind                                             ! set initial value
    force     = 0.0_r_kind                                             ! set initial value
     
    if( .not. is_time_domian ) return ;                                !! time domain analysis
    
    if( allocated( dead_weight_force ) ) then
        deallocate( dead_weight_force , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating dead_weight_force vector in physics module.' 
    end if
    
    if( allocated( surface_pressure_force ) ) then
        deallocate( surface_pressure_force , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating surface_pressure_force vector in physics module.' 
    end if
    
     
    allocate( dead_weight_force( size( dictionary ) - 1 ) , stat = state )
    if( state /= 0 ) stop 'Error: failed in allocating dead_weight_force vector in physics module.' 
    
    allocate( surface_pressure_force( size( dictionary ) - 1 ) , stat = state )
    if( state /= 0 ) stop 'Error: failed in allocating surface_pressure_force vector in physics module.' 
    
    
    dead_weight_force(:)      = 0.0_r_kind ;
    surface_pressure_force(:) = 0.0_r_kind ;
    
    
end subroutine initialize_global_matrices


    
!!====================================================================!!
!!
!! add elemental matrices to global matrices
!! input:
!!       local_matrix: is elemental matrices, either mass, stiffness or damping. 
!!                     the size of this matrix is num_row * num_col
!!
!!       id_of_dest_matrix : index identifying destination matrix whose values are:
!!                           1. destination matrix is global mass matrix
!!                           2. destination matrix is global stiffness matrix
!!                           3. destination matrix is global damping matrix
!!
!!       row_index : a vector specifying rows of glabal matrix that new elements should be added.
!!                   each row shows one of variables.
!!  
!!       num_row : number of variables in row_index ( the default size of row_index is num_max_variables_in_eqns ) 
!!
!!       col_index : a vector specifying columns of global matrix that new elements should be added.
!!                   
!!       num_col : number of variables in col_index.
!!
!!   Question 1. How to add elements in global matrices?
!!   Answer 1. If row_index(i) .neq. 0 then we search for col_index(k) .neq. 0 and 
!!             add local_matrix( i , k ) to destination_matrix( row_index(i) , col_index(k) ).
!!             The destination_matrix is either of global mass matrix, global stiffness matrix, or global damping matrix.
!!
!!
!!   Question 2. What are the meaning of row_index and col_index in elemental equations?
!!   Answer 2.   The answer is clear from an example. Assume that the equation goverining the fluid is
!!               G \ddot(p) + H p + B \ddot(r) = -B * J * a_g
!!               where, for Q8 element, G and H are 8 *8 matrces and B is 8*16 matrix.
!!               These matrices are linear maps acting on a domain vector and returning an image vector.
!!               The domain is input of these maps and image vector is output of the maps.
!!
!!               For the first term: input is \ddot(p) , so col_index is index of pressure dofs of nodes in the element.
!!                                   the matix G maps col_index to row_index where is, here, the same index of pressure dof.
!!                                   So for the first term col_index = row_index.
!!
!!               For the second term: input is p , so col_index (input) is pressure of dofs of nodesin element.
!!                                    The matrix H also maps col_index to row_index that is here the same, too.
!! 
!!               For the third term: input is \ddot(r), so col index is index of translational dofs of nodes in the element.
!!                                   output is along pressure dofs. Here, num_rows is 8 pressure dofs of nodes
!!                                   while num_cols is 16 translational dofs of the same nodes.
!!                           
!!               This example shows that row_index and col_index can be the same in many situations.
!!               However, in general this is not true and an example is dam-reservior boundary.
!! 
!!  
!!====================================================================!!

  subroutine  add_to_global_matrices( local_matrix , id_of_dest_matrix , row_index , num_row , col_index , num_col )
    
    implicit none
    real( kind = r_kind ) , dimension(:,:) , intent( in ) :: local_matrix
    integer               , intent( in ) ::  id_of_dest_matrix , num_row , num_col
    integer , dimension(:) , intent( in ) :: row_index , col_index
    

    integer :: i , k , var_location  , prev_var_id 
    real( kind = r_kind ) :: effective_epsilon
    
    
    if( ( num_row < 1 ) .or. ( num_col < 1 ) ) return
    
    effective_epsilon = 10.0_r_kind * epsilon( local_matrix(1,1) )
    
    
    !! set destination matrix
    
    select  case ( id_of_dest_matrix )
       case ( 1 )                                                       !! add to mass matrix
       
          if( .not. associated( destination_matrix , mass ) ) then
              if( associated( destination_matrix )  ) then
                  nullify( destination_matrix )
              end if
              destination_matrix => mass  
          end if
       
       case ( 2 )                                                       !! add to stiffness matrix
       
          if( .not. associated( destination_matrix , stiffness ) ) then
              if( associated( destination_matrix )  ) then
                  nullify( destination_matrix )
              end if
              destination_matrix => stiffness 
          end if
             
       
       case ( 3 )                                                       !! add to damping matrix
       
          if( .not. associated( destination_matrix , damping ) ) then
              if( associated( destination_matrix )  ) then
                  nullify( destination_matrix )
              end if
              destination_matrix => damping  
          end if
       
       case default
           stop 'Error: invalid id for destination matrix in physics module' 
    end select 
        
    
    !! add local matrix to destination matrix
    
    do  i = 1 , num_row                                                 ! loop over rows of local matrix 
        
       if( row_index( i ) .eq.  0 ) cycle                               ! i-th variable is not involved in the matrix
       
       prev_var_id  = -2                                                ! it begins from negative number not to satisfy first iteration
       var_location =  0
       
       do k = 1 , num_col                                               ! loop over columns of local matrix
          
          if( col_index( k ) .eq. 0 ) cycle                             ! k-th variable is not involved in the matrix
          
          if( abs( local_matrix( i , k ) ) < effective_epsilon ) cycle ! the value is physically zero
                
                
          if( col_index( k )  >  prev_var_id + 1  ) then
        
              call get_variable_location( row_index( i ) , col_index( k ) , var_location , .false. )
                      
          else if( col_index( k ) < prev_var_id + 1 ) then
                         
              call get_variable_location( row_index( i ) , col_index( k ) , var_location , .true.  )
            
          else
              var_location = var_location + 1  
              
          end if
          
          prev_var_id = col_index( k )
          
          destination_matrix( var_location ) = destination_matrix( var_location ) + local_matrix( i , k ) 
          
       end do ! over k

    end do ! over i
   
  
  end subroutine add_to_global_matrices


!!====================================================================!!
!!
!! add elemental force vector to resultant force
!! format of elemental force vector is :
!!
!! row1 : f_x_1 , f_y_1 , f_z_1 
!! row2 : f_x_2 , f_y_2 , f_z_2
!! row3 : f_x_3 , f_y_3 , f_z_3
!! row4 : f_x_4 , f_y_4 , f_z_4 , so on.
!! 
!! The set { row1 , row2 , row3 , ... } is collected in row_index.
!!
!! if the row_index is nonzero then its value shows variable index.
!! if it is zero then there is no need to add that entry to global force.
!! 
!!====================================================================!!


  subroutine add_to_resultant_force( local_force_vector, row_index , num_var )
    
    use  mod_geometry , only : num_dim
    implicit none
    real( kind = r_kind ) , dimension(:,:) , intent( in ) :: local_force_vector
    integer               , dimension(:)   , intent( in ) :: row_index
    integer               , intent( in ) :: num_var 

    integer :: i 

    
    do i = 1 , num_var
    
       if( row_index(i) .eq. 0 ) cycle                                  ! no force for this variable.
       
       force( row_index(i) , 1:num_dim ) = force( row_index(i) , 1:num_dim ) + local_force_vector( i , 1:num_dim )
       
    end do 
  
  
  end subroutine add_to_resultant_force
  
  
!====================================================================!!
!!  get_variable_location in dependency array
!!  get ( row , col ) element of a generic matrix 
!!  return location in dictionary 
!!
!!  The search algorithm is binary in first step and then linear as data set is small
!!  In 2D problem search space is then about 20 and in 3D it is about 150
!!====================================================================!!

  
 subroutine  get_variable_location( row , col , var_loc , left_side_search ) 
     
    implicit none
    integer , intent( in ) :: row , col 
    integer , intent( inout ) :: var_loc
    logical , intent( in ) :: left_side_search
     
    integer :: low , high
     
     
    !! suitable section of index interval
     
    if( left_side_search ) then        
                               
        low  = dictionary( row ) + 1                                   ! lower bound in searching interval
        high = min( dictionary( row + 1 ) , var_loc )                  ! upper bound in searching interval
         
    else
         
        low  = max( dictionary( row ) + 1 , var_loc )
        high = dictionary( row + 1 ) 
         
    end if
     
     !! use binary search in second step to divide the section into two subsection
     
     var_loc  = ( low + high ) / 2
     
     if( dependency( var_loc ) .eq. col ) return
     
     if( dependency( var_loc )   <  col ) then
           low  = var_loc + 1
     else 
           high = var_loc - 1
     end if
     
     !! use a linear search in third step to get the location
     do var_loc = low  , high
        if( dependency( var_loc ) .eq. col ) return
     end do
     
     write(*,*) 'row , col = ' , row , col
     stop 'Error: failed in finding variable location in physics module.' 
     
     
   end subroutine get_variable_location


!!====================================================================!!
!!
!!  find_DRM_internal_face
!!
!!====================================================================!!  

 subroutine find_DRM_internal_face()
    use  mod_geometry , only : num_dim , node, element_matrix
    implicit none
    integer :: i,j,k,n1, n2, d1, d2  , mm, node_base(3) , node_new(3)
    real( kind = r_kind) :: x1  , x2, xmin , xmax
    integer, dimension(4,3) :: faces = reshape( [ 1 , 4 , 2 , 1 , &
                                                & 2 , 3 , 3 , 4 , &
                                                & 5 , 7 , 6 , 8 ] , [4,3] ) ;
    
    
    !! find number of nodes in common : every element has at most two neighbors
    if( drm_layer_ndyn .eq. 0 ) return 
    
    
    do  i = 1 ,  drm_layer_info%num_element_in_layer
        mm = 0 ;
        n1 = drm_layer_info%element_index(i)
        do  d1 = 1 , 4
            node_base(1:3) = element_matrix( n1 , faces( d1 , 1:3) );
            
            do j = 1 , drm_layer_info%num_element_in_layer
                if( i .eq. j ) cycle 
                n2 = drm_layer_info%element_index(j)
                do  d2 = 1 , 4
                    node_new(1:3) = element_matrix( n2 , faces( d2 , 1:3) );
                    if( (node_base(1) .eq. node_new(1)) .and. (node_base(2) .eq. node_new(2)) .and. &
                      & (node_base(3) .eq. node_new(3)) ) then
                        
                        mm = mm + 1 ;                                   !! num neighbor
                        if( mm > 2 ) then
                            stop 'Error: DRM element has more than two nearby DRM elements in physics module.'
                        end if
                        
                        drm_layer_info%internal_face( i , mm ) = d1 ;
                    end if 
                end do
            end do ! over j
             
        end do ! over d1
        
        if( mm .eq. 0 ) then
            stop 'Error: DRM element has no nearby DRM elements in physics module.'
        end if
    end do ! over i
    
    
    
    !! get those elem_ids with pairs(4,0) in faces and store them in node_base(1:3) =0
    node_base(1:3) =0
    mm = 0
    
    do  i = 1 , drm_layer_info%num_element_in_layer
        if( (drm_layer_info%internal_face( i , 1 ) .eq. 4) .and. &
          & (drm_layer_info%internal_face( i , 2 ) .eq. 0) ) then
            
            mm = mm + 1
            if( mm > 2 ) then
                stop 'Error: There are more than 2 DRM elements with (4,0) common faces'
            end if
            node_base(mm) = i ;
        end if
    end do 
    
    if( mm .eq. 0 ) then !! There is no vertical set of elements: only horizontal row of elements
        do  i = 1 , drm_layer_info%num_element_in_layer
            drm_layer_info%internal_face( i , 1 ) = 3 ;          !! internal face
            drm_layer_info%internal_face( i , 2 ) = 0 ;          !! internal face
        end do
        return ;
    end if
    
    !! put smaller value: left most top element: there are only two elements
    n1 = element_matrix( drm_layer_info%element_index(node_base(1))  , 1)
    x1   = node(n1 , 1) 
    
    n2 = element_matrix( drm_layer_info%element_index(node_base(2))  , 1)
    x2   = node(n2 , 1)
    
    if( x1 < x2 ) then
        drm_layer_info%internal_face( node_base(1) , 1 ) = 1 ;          !! internal face
        drm_layer_info%internal_face( node_base(1) , 2 ) = 0 ;          !! internal face
        drm_layer_info%internal_face( node_base(2) , 1 ) = 2 ;          !! internal face 
        drm_layer_info%internal_face( node_base(2) , 2 ) = 0 ;          !! internal face
        xmin = x1 ;
        xmax = x2 ;
    else
        drm_layer_info%internal_face( node_base(1) , 1 ) = 2 ;          !! internal face
        drm_layer_info%internal_face( node_base(1) , 2 ) = 0 ;          !! internal face
        drm_layer_info%internal_face( node_base(2) , 1 ) = 1 ;          !! internal face 
        drm_layer_info%internal_face( node_base(2) , 2 ) = 0 ;          !! internal face
        
        xmin = x2 ;
        xmax = x1 ;  
    end if  
    
    
    do  i = 1 , drm_layer_info%num_element_in_layer
        if( ( i .eq. node_base(1) ) .or. ( i .eq. node_base(2) ) ) cycle ;
        
        !! vertical elements : internal face = 1 or 2: depending on x-x1 value
        if( (drm_layer_info%internal_face( i , 1 ) .eq. 3) .and. &
          & (drm_layer_info%internal_face( i , 2 ) .eq. 4) ) then
            
            !! get x location of first node
            x1 = node( element_matrix( drm_layer_info%element_index(i)  , 1) , 1 )
            if( abs(x1 - xmin) < abs(x1 - xmax) ) then
                drm_layer_info%internal_face( i , 1 ) = 1 ;          !! internal face
                drm_layer_info%internal_face( i , 2 ) = 0 ;          !! internal face
            else
                drm_layer_info%internal_face( i , 1 ) = 2 ;          !! internal face
                drm_layer_info%internal_face( i , 2 ) = 0 ;          !! internal face
            end if 
            
        !! left  corner element: face = (1,3) : no need to change
        
        else if ( (drm_layer_info%internal_face( i , 1 ) .eq. 1) .and. &
                & (drm_layer_info%internal_face( i , 2 ) .eq. 3) ) then
                
            drm_layer_info%internal_face( i , 1 ) = 1 ;          !! internal face
            drm_layer_info%internal_face( i , 2 ) = 3 ;          !! internal face
            
         
        !! right corner element: face = (2,3) : no need to change
        
        else if ( (drm_layer_info%internal_face( i , 1 ) .eq. 2) .and. &
                & (drm_layer_info%internal_face( i , 2 ) .eq. 3) ) then
                
            drm_layer_info%internal_face( i , 1 ) = 2 ;          !! internal face
            drm_layer_info%internal_face( i , 2 ) = 3 ;          !! internal face
            
            
        !! Horizontal elements : face = 3
        else if ( (drm_layer_info%internal_face( i , 1 ) .eq. 1) .and. &
                & (drm_layer_info%internal_face( i , 2 ) .eq. 2) ) then
                
            drm_layer_info%internal_face( i , 1 ) = 3 ;          !! internal face
            drm_layer_info%internal_face( i , 2 ) = 0 ;          !! internal face
        
        else
            stop 'Error: DRM layer face could not be specified'       
        end if
    end do  
     
    
 
 end subroutine find_DRM_internal_face
 
!!====================================================================!!


 subroutine  clean_up_physics()
    implicit none
    integer :: state = 0 , state1 = 0 ,i
     

     
    if( allocated(  material ) )  then
        deallocate( material , stat = state )
        state1 = state1 + state 
    end if


     
    if( allocated(  boundary_load ) ) then
        deallocate( boundary_load , stat = state )
         state1 = state1 + state 
    end if

     
    if( allocated(  element_property ) ) then
         deallocate( element_property , stat = state )
         state1 = state1 + state 
    end if
     
     
    if( allocated(  rigidity_matrix ) ) then
         deallocate( rigidity_matrix , stat = state )
         state1 = state1 + state 
    end if
     
    if( allocated(  dependency ) ) then
         deallocate( dependency , stat = state  )
         state1 = state1 + state 
    end if
     
    if( allocated(  dictionary ) ) then
         deallocate( dictionary , stat = state  )
         state1 = state1 + state 
    end if
     
    if( associated( destination_matrix ) ) then
         nullify( destination_matrix )
    end if
     
    if( allocated(  mass ) ) then
         deallocate( mass , stat = state  )
         state1 = state1 + state 
    end if
     
    if( allocated(  stiffness ) ) then
         deallocate( stiffness , stat = state  )
         state1 = state1 + state 
    end if
     
    if( allocated(  damping ) ) then
         deallocate( damping , stat = state  )
         state1 = state1 + state 
    end if
     
    if( allocated(  translational_variable ) ) then
         deallocate( translational_variable , stat = state  )
         state1 = state1 + state 
    end if
     
    if (state1 /= 0 ) stop 'Error: failed deallocating memory in physics class.'
     
    !if( analysis_kind /= 2 ) return ;                                   !! time domain analysis
    
    call  enforced_displacement%free_enforced_disp() ;
    call  boundary_springs%free_boundary_spring() ;
    call  nodal_loads%free_nodal_load() ;
    call  load_cases%free_load_case() ;
    
   
    
    if( allocated( dead_weight_force ) ) then
        deallocate( dead_weight_force , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating dead_weight_force vector in physics module.' 
    end if
    
    
    if( allocated( surface_pressure_force ) ) then
        deallocate( surface_pressure_force , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating surface_pressure_force vector in physics module.' 
    end if
    
  
    if( allocated(tension ) ) then
        deallocate( tension , stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating tension in physics module.'  
    end if
    
    
    if( allocated(  compression )) then
        deallocate( compression, stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating compression in physics module.'    
    end if
    
    
    if( drm_layer_ndyn > 0 ) then
        call drm_layer_info%free_drm_layer() ;
        
        if( allocated( drm_imported_data ) ) then
            deallocate( drm_imported_data , stat = state ) ;
            if( state /= 0 ) stop 'Failed in deallocating memory in clean_up_analysis in analysis module' ;
        end if
    end if
    
    if( allocated(  u_uddot_njtf_drm ) ) then
        deallocate( u_uddot_njtf_drm, stat = state )
        if( state /= 0 ) stop 'Error: failed in deallocating u_uddot_njtf_drm in physics module.' 
    end if 
    
  end subroutine clean_up_physics


    
      
      
  end  module  mod_physics
