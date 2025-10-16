!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!


!!
!!  This module claculates coefficient matrices in the equations.
!!  These matrices are mass, stiffness and damping matrix.
!!
!!  Also suitable set of variables are constructed and the data are stored
!!  efficiently. 
!!  !! first all variables are defined here.Depending on the model some of which are allocated.
!!  !! second every part is specialized to an abstract subroutine
!!     whose implementation is done in a specarate submodule/
!!     The aim: all variables inside this base module can be accessed through submodules.
!!     This module contains the fem_for_elastic_accoustic_and dam_res part.
!!     Only bboundary parts are implemented from elsewhere.
!!  
!!=====================================================================!!


  module  mod_fem
  
  use mod_utils    , only : r_kind 
  use mod_element
  
  implicit none
  private
  
  
  
    !! define all element types here.
    type(  t_Q8_element    ) , private , target  :: a_Q8_element          ! uninitilized instance of subclass
    type(  t_brick_element ) , private , target  :: a_brick_element       ! uninitilized instance of subclass
     
    class( t_element  )      , private , pointer :: the_element           ! pointer to an arbitrary element 
   
    
  
 
    
    integer , private :: current_num_Gauss_point = 0
    integer , private :: current_material_id = 0
    real( kind = r_kind ) , dimension(4) , public :: gauss_points          = 0.0_r_kind
    real( kind = r_kind ) , dimension(4) , public :: weight_of_gauss_point = 0.0_r_kind
    
     
  
    real( kind = r_kind ) , allocatable , private :: elemental_mass_matrix(:,:)
    real( kind = r_kind ) , allocatable , private :: elemental_stiffness_matrix(:,:)
    real( kind = r_kind ) , allocatable , private :: elemental_damping_matrix(:,:)
    real( kind = r_kind ) , allocatable , private :: elemental_force_vector(:,:)
           
    integer ,  allocatable , private :: row_variable_set(:)
    integer ,  allocatable , private :: col_variable_set(:)
    integer ,  allocatable , private :: indexed_set(:)    
    
    integer , private :: res_inf_found_node_index               ! the node common in res-inf and res-found boundaries
    integer , private :: found_inf_corner_node(4)               ! the boundary corner nodes 
     
     
    real( kind = r_kind ) , allocatable , private :: dN(:,:)           !! derivative of shape function     
    real( kind = r_kind ) , allocatable , private :: dN_rigidity(:,:)  !! product of dN and rigidity matrices  
    real( kind = r_kind ) , allocatable , private :: dN_D_dN(:,:)      !! dN * rigidity * transpose( dN )
     
     
    
    real( kind = r_kind ) , allocatable , private :: rotation(:,:)      ! Lysmer rotation matrix: only for Lysmer BC
    real( kind = r_kind ) , allocatable , private :: C_Lysmer(:,:)      ! only for Lysmer matrix
    
        
    real( kind = r_kind ) , allocatable, private :: boundary_shape_matrix(:,:)       !  only for HW model
    real( kind = r_kind ) , allocatable, private :: boundary_der_shape_matrix(:,:)   !  only for HW model
    real( kind = r_kind ) , allocatable, private :: alpha_coeffs(:)                  !  only for HW model
    
    
    !!-----------------------------------------------------------------!!
    !!  subroutines for fem analysis
    
    public   :: do_fem_analysis
    private  :: construct_fem_model
    public   :: destruct_fem_model 
    
    private  :: fem_for_acoustic_bulk 
    private  :: fem_for_elastic_bulk 
    private  :: fem_for_dam_reservior_boundary
    private  :: fem_for_reservior_foundation_boundary
    private  :: fem_for_foundation_truncation_boundary
    private  :: fem_for_reservior_truncation_boundary  
    
    private  :: find_res_found_inf_common_node                          !! used in higher order BC models. 
    private  :: find_found_inf_corner_node
    
    public  :: set_gauss_points   
    public  :: set_rigidity_matrix 
    public  :: find_strain_stress    
    private :: elemental_strain_stress
    
    public  :: force_in_drm_layer
    
    public  :: freq_domain_force_in_drm_layer
    
    interface
    
        module subroutine basic_model_elastic_bulk( elem_id , is_drm_calc )
          integer, intent( in ) :: elem_id
          logical , intent( in ) :: is_drm_calc
        end subroutine basic_model_elastic_bulk
        
        
        module subroutine elastic_surface_pressure_force( elem_id  , bc_id )
          integer, intent( in ) :: elem_id
          integer, intent( in ) :: bc_id
        end subroutine elastic_surface_pressure_force
        
       
    
       module subroutine basic_model_acoustic_bulk( elem_id )
          integer, intent( in ) :: elem_id
       end subroutine basic_model_acoustic_bulk
       
       
       module subroutine basic_model_dam_reservior_boundary( elem_id , bc_id )
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine basic_model_dam_reservior_boundary
       
       
       module subroutine basic_model_reservior_foundation_boundary( elem_id , bc_id )
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine  basic_model_reservior_foundation_boundary
       
       
       
       
       !!-------------------------------------------------------------!!
       !! Models of boundary condition in foundation truncation boundary:
       
       
       module subroutine basic_model_foundation_truncation_boundary( elem_id , bc_id )       !! Lysmer
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine basic_model_foundation_truncation_boundary
       
       
       
       module subroutine HW_model_foundation_truncation_boundary( elem_id , bc_id )           !! HW found-inf
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine HW_model_foundation_truncation_boundary
       
       
       
       !!-------------------------------------------------------------!!
       !! Models of boundary condition in reservior truncation boundary:
       
       
       module subroutine basic_model_reservior_truncation_boundary( elem_id , bc_id )        !! Somerfeld
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine basic_model_reservior_truncation_boundary
       
       module subroutine HW_model_reservior_truncation_boundary( elem_id , bc_id )           !! HW-res-inf
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine HW_model_reservior_truncation_boundary
       
       
       module  subroutine GN_model_reservior_truncation_boundary( elem_id , bc_id )           !! GN-res-inf
          integer, intent( in ) :: elem_id 
          integer, intent( in ) :: bc_id
       end subroutine GN_model_reservior_truncation_boundary
       
       
       
    end interface
    
    
  
  contains
  
  
   



!!=====================================================================!!
!!
!!  the fem analysis
!!
!!=====================================================================!!
    
    
    
  subroutine  do_fem_analysis()
    use mod_physics , only : the_model_description, analysis_kind
    implicit none

    call  construct_fem_model()                                         !! assign memory to needed arrays
    
    !!----------------------------------------------------------------!!
    !! Model dependent subroutines
    
    
    if( the_model_description%foundation_inf_boundary_model_id > 0 ) then  !! needs for HW models of found-inf boundary
        call  find_found_inf_corner_node()
    end if
    
    if( the_model_description%reservior_inf_boundary_model_id > 0 ) then
       call find_res_found_inf_common_node() ;
    end if
    
    
    
    !!----------------------------------------------------------------!!
    !! Model independent subroutines
    
    call  fem_for_elastic_bulk()
    
    call  fem_for_acoustic_bulk()
    
    call  fem_for_dam_reservior_boundary()
    
    call  fem_for_reservior_foundation_boundary()
    
    call  fem_for_foundation_truncation_boundary()
    
    call  fem_for_reservior_truncation_boundary()
    
     
  end subroutine  do_fem_analysis
 
 
!!====================================================================!!
!!
!! force in DRM layer in freq domain analysis
!!
!! 
!!====================================================================!!

    subroutine freq_domain_force_in_drm_layer( omega , a_g , beta, dir )
        use  mod_geometry , only : element_matrix , num_dim, degree_of_freedom
        use  mod_physics  , only : drm_layer_info 
        use  mod_solver   , only : cplx_response
        
        implicit none
        real( kind = r_kind) , intent(in) :: omega
        real( kind = r_kind) , intent(in) :: beta                       ! hysteresis_damping_coeff
        real( kind = r_kind) , intent(in) :: a_g
        integer , intent( in ) :: dir                                   ! earthquake_direction
        
        integer :: i,j,k, n, num_var , elem_id , is_boundary_var(16) 
        real( kind = r_kind ) :: tmp_re
        complex( kind = r_kind ) :: tmp_cplx
        
        tmp_re = 1.0_r_kind / ( omega ** 2) ;
        
        tmp_cplx = cmplx( tmp_re , 2.0_r_kind * tmp_re * beta )
        
        
        num_var = the_element%num_node_in_element * num_dim 
        
     
        do  elem_id = 1 ,  drm_layer_info%num_element_in_layer
        
            call  basic_model_elastic_bulk( drm_layer_info%element_index(elem_id) , .true. ) 
              
            
            !! find boundary nodes
            is_boundary_var(:) = 0 ;
            
            !! for nodes in internal face
            if(  ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 1 )   .and. &
              &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 0 ) ) then
                
                is_boundary_var(1:4)  = 1 ;   !! node 1 ,2
                is_boundary_var(9:10) = 1 ;   !! node 5
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 2 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 0 ) ) then
                
                is_boundary_var(5:8)   = 1 ;   !! node 3 ,4
                is_boundary_var(13:14) = 1 ;   !! node 7
                
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 3 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 0 ) ) then
                
                is_boundary_var(3:6)   = 1 ;   !! node 2 ,  3
                is_boundary_var(11:12) = 1 ;   !! node 6
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 1 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 3 ) ) then
                
                is_boundary_var(3:4)   = 1 ;   !! node 2 
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 2 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 3 ) ) then
                
                is_boundary_var(5:6)   = 1 ;   !! node 3 
            else
                stop 'Error in specifying DRM face in fem module '
            end if 
            
            
            !! calculate effective forces: Mbar J, put in effective force
            elemental_force_vector(:,:) = 0.0_r_kind                    !! set force zero
            
            do  i = 1 , num_var
                if( is_boundary_var(i) .eq. 1 ) then
                    do  j = 1 , num_var                                 !! calculate force over non-boundary nodes
                        if( is_boundary_var(j) .eq. 1 ) cycle 
                        if( mod(j , 2) .eq. 1 ) then                    !! calculate -Mbar * J for M_{b,e} * J
                            elemental_force_vector(i,1) = elemental_force_vector(i,1) - elemental_mass_matrix(i,j)
                        else
                            elemental_force_vector(i,2) = elemental_force_vector(i,2) - elemental_mass_matrix(i,j)
                        end if
                    end do 
                else
                    do  j = 1 , num_var                                 !! calculate force over boundary nodes
                        if( is_boundary_var(j) .eq. 0 ) cycle 
                        if( mod(j , 2) .eq. 1 ) then                    !! calculate + Mbar * J for M_{e,b} * J
                            elemental_force_vector(i,1) = elemental_force_vector(i,1) + elemental_mass_matrix(i,j)
                        else
                            elemental_force_vector(i,2) = elemental_force_vector(i,2) + elemental_mass_matrix(i,j)
                        end if
                    end do 
                end if 
                
            end do  ! over i
            
            
            
            !! add Mbar * J to response vector
            do  i = 1 , num_var
                n = row_variable_set(i)
                if( n .eq. 0 ) cycle
                cplx_response(n) = cplx_response(n) + cmplx( elemental_force_vector(i,dir) * a_g , 0.0_r_kind)
            end do 

            
            !! calculate effective forces: K * J, put in effective force again
            elemental_force_vector(:,:) = 0.0_r_kind                    !! set force zero
            
            do  i = 1 , num_var
                if( is_boundary_var(i) .eq. 1 ) then
                    do  j = 1 , num_var                                 !! calculate force over non-boundary nodes
                        if( is_boundary_var(j) .eq. 1 ) cycle 
                        if( mod(j , 2) .eq. 1 ) then                    !! calculate -K * J for K_{b,e} * J
                            elemental_force_vector(i,1) = elemental_force_vector(i,1) + &
                        &   elemental_stiffness_matrix(i,j)
                        else
                            elemental_force_vector(i,2) = elemental_force_vector(i,2) + &
                        &   elemental_stiffness_matrix(i,j)
                        end if
                    end do 
                else
                    do  j = 1 , num_var                                 !! calculate force over boundary nodes
                        if( is_boundary_var(j) .eq. 0 ) cycle 
                        if( mod(j , 2) .eq. 1 ) then                    !! calculate + K * J for K_{e,b} * J
                            elemental_force_vector(i,1) = elemental_force_vector(i,1) - &
                        &   elemental_stiffness_matrix(i,j)
                        else
                            elemental_force_vector(i,2) = elemental_force_vector(i,2) - &
                        &   elemental_stiffness_matrix(i,j)
                        end if
                    end do 
                end if 
                
            end do  ! over i
            
            
            !! add Kbar * J to response vector
            do  i = 1 , num_var
                n = row_variable_set(i)
                if( n .eq. 0 ) cycle
                cplx_response(n) = cplx_response(n) + tmp_cplx *  elemental_force_vector(i,dir) * a_g
            end do 
            
        end do ! over elem_id
        
    end subroutine freq_domain_force_in_drm_layer

!!====================================================================!!
!!
!! calculate forces of DRM Layer
!!
!!====================================================================!!

    subroutine force_in_drm_layer( istep , file_a_id , file_b_id , alpha_m , alpha_k )
        use  mod_geometry , only : element_matrix , num_dim, degree_of_freedom , node
        use  mod_physics  , only : drm_layer_info , drm_imported_data , drm_layer_ndyn , load_cases
        use  mod_solver   , only : response 
        
        implicit none
        integer , intent(in) :: istep 
        integer , intent(in) :: file_a_id , file_b_id    !! input files
        real( kind = r_kind) , intent(in) :: alpha_m
        real( kind = r_kind) , intent(in) :: alpha_k
        
        real( kind = r_kind) :: tmp_re , ue_vec(16) , ve_vec(16) , ae_vec(16) , a_g = 9.810_r_kind 
        integer :: i,j,k, n, num_var , elem_id, elem_nodes(8) , is_boundary_var(16) 
        
        
        num_var = the_element%num_node_in_element * num_dim
        
        if( drm_layer_ndyn .eq. 1 ) then
            !! allocate required matrix 
            if( .not. allocated( drm_imported_data ) ) then 
                allocate( drm_imported_data(drm_layer_info%num_node_in_layer , 6) , stat = i )
                if( i /=0 ) stop 'Error in allocating drm_imported_data array in fem module '
            end if 
        
            !! read a block of free field data
            drm_imported_data(:,:) = 0.0_r_kind ;
        
            do  i = 1 , drm_layer_info%num_node_in_layer
                read( file_a_id , *) j , drm_imported_data(i,1:6)
            end do 
            
        else if( drm_layer_ndyn .eq. 2 ) then
        
            !! allocate required matrix 
            if( .not. allocated( drm_imported_data ) ) then 
                allocate( drm_imported_data(drm_layer_info%num_cols , 2) , stat = i )
                if( i /=0 ) stop 'Error in allocating drm_imported_data array in fem module '
            end if 
            
            !! read a block of free field data
            drm_imported_data(:,:) = 0.0_r_kind ;
            
            if( file_a_id > 0 ) then
                read( file_a_id , *) drm_imported_data( 1:drm_layer_info%num_cols , 1)  !! read x data
            end if
            
            if( file_b_id > 0 ) then
                read( file_b_id , *) drm_imported_data( 1:drm_layer_info%num_cols , 2)  !! read y data
            end if
            
        end if
         
        
        
     
        do  elem_id = 1 ,  drm_layer_info%num_element_in_layer
        
            call  basic_model_elastic_bulk( drm_layer_info%element_index(elem_id) , .true. ) 
            
    elemental_damping_matrix(1:num_var,1:num_var) = alpha_m * elemental_mass_matrix( 1:num_var,1:num_var) &
                                            & + alpha_k * elemental_stiffness_matrix(1:num_var,1:num_var)
            
            !!  find nodal elements
            elem_nodes(:) = 0 ;
            
            if( drm_layer_ndyn .eq. 1 ) then
                do  i = 1 , the_element%num_node_in_element
                    j = element_matrix( drm_layer_info%element_index(elem_id) , i)
               
                    !! search in pair node index
                    do  k = 1 , drm_layer_info%num_node_in_layer
                        if( j .eq. drm_layer_info%paired_node_index(k,2) ) then
                            elem_nodes(i) = k ;
                            exit ;
                        end if
                    end do 
                end do 
                
            else if( drm_layer_ndyn .eq. 2 ) then
            
                do  i = 1 , the_element%num_node_in_element
                    j = element_matrix( drm_layer_info%element_index(elem_id) , i)
               
                    !! compare y-value with levels 
                    do  k = 1 , drm_layer_info%num_levels
                        tmp_re = abs( abs( node(j,2) ) - real( drm_layer_info%levels(k) , kind = r_kind) ) ;! diff of depth
                        if( tmp_re < 1.0_r_kind + 1.0e-9 ) then
                            elem_nodes(i) = k
                            exit 
                        end if  
                    end do 
                end do 
            
            end if 
            
            
            !! check if all nodes are found
            do  i = 1 , the_element%num_node_in_element
                if( elem_nodes(i) .eq. 0 ) stop 'Error in matching element array in fem module '
            end do 
            
            
            
            
            !! find boundary nodes
            is_boundary_var(:) = 0 ;
            
            !! for nodes in internal face
            if(  ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 1 )   .and. &
              &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 0 ) ) then
                
                is_boundary_var(1:4)  = 1 ;   !! node 1 ,2
                is_boundary_var(9:10) = 1 ;   !! node 5
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 2 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 0 ) ) then
                
                is_boundary_var(5:8)   = 1 ;   !! node 3 ,4
                is_boundary_var(13:14) = 1 ;   !! node 7
                
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 3 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 0 ) ) then
                
                is_boundary_var(3:6)   = 1 ;   !! node 2 ,  3
                is_boundary_var(11:12) = 1 ;   !! node 6
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 1 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 3 ) ) then
                
                is_boundary_var(3:4)   = 1 ;   !! node 2 
                
            else if( ( drm_layer_info%internal_face( elem_id , 1 ) .eq. 2 )   .and. &
                  &  ( drm_layer_info%internal_face( elem_id , 2 ) .eq. 3 ) ) then
                
                is_boundary_var(5:6)   = 1 ;   !! node 3 
            else
                stop 'Error in specifying DRM face in fem module '
            end if 
            
             
            
            !! put u , udot and uddot in elemental vectors
            ue_vec(:) = 0.0_r_kind 
            ve_vec(:) = 0.0_r_kind 
            ae_vec(:) = 0.0_r_kind 
            
            if( drm_layer_ndyn .eq. 1 ) then !! read from free field data
                do  i = 1 , the_element%num_node_in_element
                    j = 2 *(i-1)
                    ue_vec( j + 1:j + 2) = drm_imported_data( elem_nodes(i),1:2)       !! u_x , uy
                    ve_vec( j + 1:j + 2) = drm_imported_data( elem_nodes(i),3:4)       !! udot_x ,  udot_y
                    ae_vec( j + 1:j + 2) = drm_imported_data( elem_nodes(i),5:6) * a_g !! uddot_x , uddot_y
                end do 
            else    !! 
                
                do  i = 1 , the_element%num_node_in_element
                    j = 2 *(i-1)
                    if( elem_nodes(i) .eq. 1 ) then                     !! node at free surface
                        ae_vec( j + 1 ) = load_cases%load_coeffs(5) * load_cases%time_load(istep,6)  !! uddot_x
                        ae_vec( j + 2 ) = load_cases%load_coeffs(6) * load_cases%time_load(istep,7)  !! uddot_y
                        ue_vec( j + 1 ) = drm_imported_data( 4,1)        !!  u_x 
                        ue_vec( j + 2 ) = drm_imported_data( 4,2)        !!  u_y
                        ve_vec( j + 1 ) = drm_imported_data( 3,1)        !!  v_x 
                        ve_vec( j + 2 ) = drm_imported_data( 3,2)        !!  v_y
                    else
                        k = 3*( elem_nodes(i) - 2 ) + 1
                        ue_vec( j + 1 ) = drm_imported_data( k + 3 , 1 )  !! u_x
                        ue_vec( j + 2 ) = drm_imported_data( k + 3 , 2 )  !! u_y
                        ve_vec( j + 1 ) = drm_imported_data( k + 2 , 1 )  !! v_x
                        ve_vec( j + 2 ) = drm_imported_data( k + 2 , 2 )  !! v_y
                        ae_vec( j + 1 ) = drm_imported_data( k + 1 , 1 ) * a_g  !! a_x
                        ae_vec( j + 2 ) = drm_imported_data( k + 1 , 2 ) * a_g  !! a_y
                    end if
                end do 
                
            end if
            
            
            !! calculate effective forces
            elemental_force_vector(:,:) = 0.0_r_kind                    !! set force zero
            
            do  i = 1 , num_var
                tmp_re = 0.0_r_kind 
                
                if( is_boundary_var(i) .eq. 1 ) then
                    do  j = 1 , num_var                                 !! calculate force over non-boundary nodes
                        if( is_boundary_var(j) .eq. 1 ) cycle 
                        tmp_re = tmp_re + elemental_mass_matrix(i,j)      * ae_vec(j) ; !! M * a 
                        tmp_re = tmp_re + elemental_damping_matrix(i,j)   * ve_vec(j) ; !! C * v 
                        tmp_re = tmp_re + elemental_stiffness_matrix(i,j) * ue_vec(j) ; !! K * u 
                    end do 
                    elemental_force_vector(i,1) = - tmp_re
                    
                else
                    do  j = 1 , num_var                                 !! calculate force over boundary nodes
                        if( is_boundary_var(j) .eq. 0 ) cycle 
                        tmp_re = tmp_re + elemental_mass_matrix(i,j)      * ae_vec(j) ; !! M * a 
                        tmp_re = tmp_re + elemental_damping_matrix(i,j)   * ve_vec(j) ; !! C * v 
                        tmp_re = tmp_re + elemental_stiffness_matrix(i,j) * ue_vec(j) ; !! K * u 
                    end do 
                    elemental_force_vector(i,1) =  tmp_re
                end if 
                
            end do  ! over i
            
            
                

            !! add to response vector
            do  i = 1 , the_element%num_node_in_element
                j = element_matrix( drm_layer_info%element_index(elem_id) , i) ! node id
                do  k = 1 , num_dim
                    n = degree_of_freedom( j , 2*k)
                    if( n /= 0) then
                        response(n) = response(n) + elemental_force_vector( 2*(i-1) + k,1)
                    end if
                end do 
            end do 
            
        end do 
        
        
    end subroutine force_in_drm_layer

!!=====================================================================!!
!!
!! calculate stress and strain data
!!
!!=====================================================================!!

    subroutine  find_strain_stress( istep , time , is_data_printed , analysis_kind, num_region_element )
        use  mod_geometry , only : node , element_matrix, num_element, num_node, num_dim, degree_of_freedom
        use  mod_physics  , only : material , element_property, tension, compression
        use  mod_solver   , only : solution 
        
        implicit none
        integer             , intent(in) :: istep , analysis_kind, num_region_element(3)
        real(kind = r_kind) , intent(in) :: time  ;
        logical             , intent(in) :: is_data_printed
        real( kind = r_kind) :: tmp_re , small_u_vec(60) , position_vec(60)
        integer :: i,j,k, d, n , elem_id, n_elem , igauss, region_id
       
        !! 1. print some data
        if(is_data_printed) then
            write(21,'(A16,F7.2,A14,I5,A1)') 'STRESSES AT TIME',time,'   :   (ISTEP=',istep,')'
        end if
        
        igauss = 0 ;
        
        
        do  elem_id = 1 , size( element_property , 1 )
            
            region_id = element_property(elem_id , 1 ) ;
            
            !! 2. extract data from position and solution
            !! 2.1. get position
            position_vec(:) = 0.0_r_kind ;
            
            !! loop over all nodes in element
            n = 0 ;
            do  i = 1 , size( element_matrix , 2)  
                j = element_matrix( elem_id , i) ;                      !! node_id
                do  d = 1 , num_dim                                     !! loop over dimension
                    n = n + 1 ;
                    if( j /= 0 ) then
                        position_vec(n) = node(j,d) ;                   !! get node position 
                    end if
                end do
            end do
                
                
            !! 2.2. get solution
            small_u_vec(:)  = 0.0_r_kind ;
            
            if( region_id < 3 ) then  !! the element is elastic
                
                n = 0 ;
                do  i = 1 , size( element_matrix , 2)                   !! loop over all nodes in element 
                    j = element_matrix( elem_id , i) ;                  !! node_id
                    do  d = 1 , num_dim                                 !! loop over dimension
                        n = n + 1 ;
                        k = degree_of_freedom(j, 2*d) ;                 !! var_id of the dof at given node
                        if( k .eq. 0 ) cycle
                        small_u_vec(n) = solution(k) ;
                    end do
                end do
            end if
            
            
            
            call elemental_strain_stress( elem_id, region_id , small_u_vec, position_vec,  & 
                                & igauss, analysis_kind , num_region_element , is_data_printed );
            
            
      
        end do ! over el_id
        
        
   
    end subroutine find_strain_stress
  
!!=====================================================================!!
!!
!!  calculate elemental stress strain matrix
!!
!!
!!=====================================================================!!

    subroutine  elemental_strain_stress(  elem_id, region_id , small_u_vec, position_vec,  &
                                &  igauss , analysis_kind , num_region_element , is_data_printed )
                                
        use  mod_physics  , only : material , element_property, rigidity_matrix, element_type
        use  mod_physics  , only : tension, compression
        use  mod_geometry , only : num_element , num_dim , degree_of_freedom
        use  mod_solver   , only : solution 
        implicit none
        
        integer, intent(in)  :: elem_id, region_id, analysis_kind, num_region_element(3)
        real( kind = r_kind ), intent( inout ) ::  small_u_vec(60), position_vec(60)
        integer, intent( inout ) :: igauss
        logical, intent(in) :: is_data_printed
        integer :: n_gauss_point, r, s , t , k , n , d, local_elem_id , elem_gauss_id
        real( kind = r_kind ) :: d_volume , ep(6), sig(6), xyzT(3) ,s1,s2,s3, poi,x1,x2,x3,zzt
        
       
        
        !! get local element id
        if( element_property(elem_id , 1 ) .eq. 1 ) then
            local_elem_id = elem_id
        elseif( element_property(elem_id , 1 ) .eq. 2 ) then
            local_elem_id = elem_id - num_region_element(1)
        else
            local_elem_id = elem_id - num_region_element(1) - num_region_element(2)
        end if
        
        n_gauss_point = element_property(elem_id , 3 ) ;
        call  set_gauss_points( n_gauss_point ) ;
        call  set_rigidity_matrix( elem_id , region_id ) ;
        poi = material( element_property( elem_id , 2) )%poisson_ratio
       
        
        elem_gauss_id = 0
        
        if( element_type(region_id) .eq. 12 ) then
        
            if( .not. is_data_printed ) return ;
            
            zzt = 0.0_r_kind ;
            
            !! loop over gauss points
            do  r = 1 , n_gauss_point
                do  s = 1 , n_gauss_point
                    elem_gauss_id = elem_gauss_id + 1 ;
                    
                    if( elem_gauss_id .eq. 1 ) then
                        write(21,'(I12,I10)') n_gauss_point**2 , 1
                    end if 
                        
    
                    write(21,'(3I4,I2,1X,20E11.4)') local_elem_id , region_id , elem_gauss_id ,1, &
                    &    zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,  &
                    &    zzt,zzt,gauss_points(r),gauss_points(s),zzt,zzt,zzt,zzt,zzt,zzt,zzt
                              
                end do ! over s
            end do ! over r
            
            return 
        end if
        
        
        if( element_type(region_id) .eq. 14 ) then
        
            if( .not. is_data_printed ) return ;
            
            zzt = 0.0_r_kind ;
            
            !! loop over gauss points
            do  r = 1 , n_gauss_point
                do  s = 1 , n_gauss_point
                    do  t = 1 , n_gauss_point
                    elem_gauss_id = elem_gauss_id + 1 ;
                    
                    if ( elem_gauss_id .eq. 1 ) then
                            write(21,'(I12,I10)') n_gauss_point**3 , 1
                    end if 
                    
    
                    write(21,'(3I4,I2,1X,20E11.4)') local_elem_id , region_id , elem_gauss_id ,1, &
                    &    zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,  &
                    &    zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt,zzt
                    
                    end do ! over t          
                end do ! over s
            end do ! over r
            
            return 
        end if
        
        
        
        
        if( num_dim .eq. 3 ) then
            
            !! loop over gauss points
            do  r = 1 , n_gauss_point
                do  s = 1 , n_gauss_point
                    do  t = 1 , n_gauss_point
                    
                        call the_element%get_shape_function( elem_id , gauss_points(r) , &
                                                         &   gauss_points(s) , gauss_points(t) )
                   
                        call the_element%get_derivative_shape_function( elem_id , gauss_points(r) ,  &
                                                                    & gauss_points(s) , gauss_points(t) )
                   
                        !! do sum over all nodes
                        
                        igauss = igauss + 1 ;
                        elem_gauss_id = elem_gauss_id + 1
                        ep(:)  = 0.0_r_kind ;
                        sig(:) = 0.0_r_kind ;
                        xyzT(:) = 0.0_r_kind ;
                        
                        do  k = 1 , size( the_element%shape_vector )
                       
                            dN = 0.0_r_kind
                            dN( 1 , 1 ) = the_element%derivative_shape_matrix( 1 , k )  ! D_x N_k  
                            dN( 2 , 2 ) = the_element%derivative_shape_matrix( 2 , k )  ! D_y N_k
                            dN( 3 , 3 ) = the_element%derivative_shape_matrix( 3 , k )  ! D_z N_k
                       
                            dN( 1 , 4 ) = dN( 2 , 2 )                             ! D_y N_k 
                            dN( 2 , 4 ) = dN( 1 , 1 )                             ! D_x N_k 
                            dN( 2 , 5 ) = dN( 3 , 3 )                             ! D_z N_k 
                            dN( 3 , 5 ) = dN( 2 , 2 )                             ! D_y N_k 
                            dN( 1 , 6 ) = dN( 3 , 3 )                             ! D_z N_k 
                            dN( 3 , 6 ) = dN( 1 , 1 )                             ! D_x N_k 
                            
                            do  n = 1 , 6
                                do  d = 1 , num_dim
                                    ep(n) = ep(n) + dN( d , n ) * small_u_vec( (k-1) * num_dim + d ) ;
                                end do
                            end do
                            
                            do  d = 1 , num_dim
                    xyzT(d) = xyzT(d) + the_element%shape_vector(k) * position_vec( (k-1) * num_dim + d );
                            end do
                           

                        end do ! over k
                        
                        !! calculate sigma = rigidity_matrix * strain
                        do  n = 1 , 6
                            do  k = 1 , 6
                                sig(n) = sig(n) + rigidity_matrix(n,k) * ep(k) ;
                            end do
                        end do
                        
                        
                        zzt = 0.0_r_kind ;
                        
                        if( region_id /= 1 ) then
                            tension( igauss , 1:9) = 0.0_r_kind ;
                            compression( igauss , 1:9) = 0.0_r_kind ;
                        else
                            
                            call  eigsss(sig , s1 ,s2,s3 , region_id , poi) ;
                        
                            if( s1 .ge. tension( igauss,7) ) then
                                tension( igauss , 1:2)= sig(1:2)
                                tension( igauss , 4 ) = sig(3)
                                tension( igauss , 7 ) = s1
                                tension( igauss , 9 ) = s3
                            end if
            
                            if( s3 <= compression(igauss,9) ) then
                                compression( igauss , 1:6)= sig(1:6)
                                compression( igauss , 7 ) = s1
                                compression( igauss , 9 ) = s3 
                            end if
                        end if
                        
                        if ( elem_gauss_id .eq. 1 .and. is_data_printed ) then
                            write(21,'(I12,I10)') n_gauss_point**3 , 1
                        end if 
                        
                        if( is_data_printed ) then
                            write(21,'(3I4,I2,1X,20E11.4)') local_elem_id , region_id , elem_gauss_id ,1, &
                                 & sig(1:6),s1,s2,s3,zzt,zzt, gauss_points(r), gauss_points(s) ,  &
                                & gauss_points(t), zzt,zzt,zzt,zzt,zzt,zzt
                        end if
                              
                    end do ! over t
                end do ! over s
            end do ! over r
       
            
        else !! num_dim == 2
        
            !! loop over gauss points
            do  r = 1 , n_gauss_point
                do  s = 1 , n_gauss_point
                   
                        call the_element%get_shape_function( elem_id , gauss_points(s) ,  gauss_points(r)  )
                   
                        call the_element%get_derivative_shape_function( elem_id, gauss_points(s), gauss_points(r))
                   
                        !! do sum over all nodes
                        
                        igauss = igauss + 1 ;
                        elem_gauss_id = elem_gauss_id + 1
                        
                        ep(:)  = 0.0_r_kind ;
                        sig(:) = 0.0_r_kind ;
                        xyzT(:) = 0.0_r_kind ;
                        
                        do  k = 1 , size( the_element%shape_vector )
                       
                            dN = 0.0_r_kind
                            dN( 1 , 1 ) = the_element%derivative_shape_matrix( 1 , k )    ! D_x N_k  
                            dN( 2 , 2 ) = the_element%derivative_shape_matrix( 2 , k )    ! D_y N_k 
                    
                    
                            dN( 1 , 3 ) = dN( 2 , 2 )                           ! D_y N_k
                            dN( 2 , 3 ) = dN( 1 , 1 )                           ! D_x N_k
                    
                            
                            do  n = 1 , 3
                                do  d = 1 , num_dim
                                    ep(n) = ep(n) + dN( d , n ) * small_u_vec( (k-1) * num_dim + d ) ;
                                end do
                            end do
                            
                            do  d = 1 , num_dim
                    xyzT(d) = xyzT(d) + the_element%shape_vector(k) * position_vec( (k-1) * num_dim + d );
                            end do
                           

                        end do ! over k
                        
                        !! calculate sigma = rigidity_matrix * strain
                        do  n = 1 , 6
                            do  k = 1 , 6
                                sig(n) = sig(n) + rigidity_matrix(n,k) * ep(k) ;
                            end do
                        end do
                        
                        
                        
                        zzt = 0.0_r_kind ;
                        
                         if( region_id /= 1 ) then
                            tension( igauss , 1:9) = 0.0_r_kind ;
                            compression( igauss , 1:9) = 0.0_r_kind ;
                        else
                            call  eigsss(sig , s1 ,s2,s3 , region_id , poi) ;
                            if( (s1 > tension( igauss,7)) .and. (region_id .eq. 1) ) then
                                tension( igauss , 1:2)= sig(1:2)
                                tension( igauss , 4 ) = sig(3)
                                tension( igauss , 7 ) = s1
                                tension( igauss , 8 ) = s2
                                tension( igauss , 9 ) = s3
                            end if
            
                            if( (s3 <= compression(igauss,9)) .and. (region_id .eq. 1) ) then
                                compression( igauss , 1:2)= sig(1:2)
                                compression( igauss , 4 ) = sig(3)
                                compression( igauss , 7 ) = s1
                                compression( igauss , 8 ) = s2
                                compression( igauss , 9 ) = s3 
                            end if
                        end if
                        
                        
                        if ( (elem_gauss_id .eq. 1) .and. is_data_printed ) then
                            write(21,'(I12,I10)') n_gauss_point**2 , 1
                        end if 
                        
                        if( is_data_printed ) then
                            write(21,'(3I4,I2,1X,20E11.4)') local_elem_id , region_id , elem_gauss_id ,1, &
                            &    sig(1:2),zzt,sig(3),zzt , zzt, s1,s2,s3,zzt,zzt, gauss_points(s), &
                            &    gauss_points(r) , zzt,zzt,zzt,zzt,zzt,zzt,zzt
                        end if      
                      
                    
                end do ! over s
            end do ! over r
        
        
        end if !! num_dim
        
    
    end subroutine elemental_strain_stress  


!!=====================================================================!!
!!
!! EIGSSS subroutine
!!
!!=====================================================================!!    
    
    subroutine eigsss( sig , s1 ,s2,s3 , region_id , poi )
    use  mod_physics  , only : element_type
    use  mod_geometry , only : num_element , num_dim , degree_of_freedom
    use  mod_solver   , only : solution 
    implicit none
    
    real( kind = r_kind ), intent(inout) :: sig(6),s1,s2,s3
    integer, intent(in) :: region_id 
    real( kind = r_kind ), intent(in) :: poi
    real( kind = r_kind ) :: error , A(3,3) , AMS(3,3), B2, BAR,BETA,COEFF,S,C,CS,SC
    integer :: i, j, k , N , IELTYPE , I1,I3
    
    error = 1.0e-6
    N = 3
    s1 = -1.0e5
    s2 = 0.0_r_kind
    s3 = 1.0e5 ;
    A(:,:) = 0.0_r_kind
    IELTYPE = element_type(region_id)
    
    if ( IELTYPE .eq. 6 ) then
        A(1,1)=SIG(1)
        A(2,2)=SIG(2)
        A(1,2)=SIG(3)
        A(2,1)=SIG(3)
    elseif( IELTYPE .eq. 5) then
        A(1,1)=SIG(1)
        A(2,2)=SIG(2)
        A(3,3)= poi *(SIG(1)+SIG(2))
        A(1,2)=SIG(3)
        A(2,1)=SIG(3)
    elseif ( (IELTYPE /= 5 ) .and. (IELTYPE /= 6 ) ) then
        A(1,1)=SIG(1)
        A(2,2)=SIG(2)
        A(3,3)=SIG(3)
        A(1,2)=SIG(4)
        A(2,1)=SIG(4)
        A(1,3)=SIG(6)
        A(3,1)=SIG(6)
        A(2,3)=SIG(5)
        A(3,2)=SIG(5)
    end if
    
    
    AMS(:,:) = 0.0_r_kind
    do i =1,N
       AMS(i,i) = 1.0_r_kind
    end do
    
    B2 = 0.0_r_kind
    do i=1,N
        do j=1,N
          if (i.NE.j) B2 = B2 + A(i,j)**2
        end do 
    end do
    
    if (IELTYPE.EQ.6 .OR. IELTYPE.EQ.5) N=2                  !!!
    
    if( B2 > ERROR) then
        BAR = 0.50_r_kind  * B2 /real(N*N , kind = r_kind) ;
        DO WHILE (B2 .GT. ERROR)
            DO  I=1,N-1
                DO J=I+1,N
                    IF (A(J,I)**2 <= BAR) CYCLE  ! do not touch small elements
                    B2 = B2 - 2.0_r_kind *A(J,I)**2
                    BAR = 0.50_r_kind * B2 /real(N*N , kind = r_kind)
                    BETA = (A(J,J)-A(I,I))/(2.0_r_kind *A(J,I))
                    COEFF = 0.50_r_kind * BETA / SQRT(1.0_r_kind + BETA**2)
                    S = SQRT(max(0.50_r_kind + COEFF , 0.0_r_kind))
                    C = SQRT(max(0.50_r_kind - COEFF , 0.0_r_kind))
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    DO K=1,N
                        CS= C*A(I,K)+S*A(J,K)
                        SC=-S*A(I,K)+C*A(J,K)
                        A(I,K) = CS
                        A(J,K) = SC
                    END DO
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    DO  K=1,N
                        CS =  C*A(K,I)+S*A(K,J)
                        SC = -S*A(K,I)+C*A(K,J)
                        A(K,I) = CS
                        A(K,J) = SC
                        CS =  C*AMS(K,I)+S*AMS(K,J)
                        SC = -S*AMS(K,I)+C*AMS(K,J)
                        AMS(K,I) = CS
                        AMS(K,J) = SC
                    END DO
                END DO
            END DO
        END DO
    end if
    

    N=3
    I1 = 0
    I3 = 0
    do I=1,N
        if (S1.LE.A(I,I)) then
            S1=A(I,I) ;
            I1 = I ;
        end if
        
        if( S3.GT.A(I,I) ) then
            I3=I
            S3=A(I,I)
        end if
    end do
    
    
    do  I=1,N
        IF (I.NE.I1 .AND. I.NE.I3) S2=A(I,I)
    end do


    end subroutine eigsss
  
  
!!===========================================================================!!
!!
!!  submodule subroutine: fem_for_elastic_bulk
!!
!!===========================================================================!!

    subroutine fem_for_elastic_bulk()
       
        use mod_physics , only  : element_property, is_time_domian , boundary_load , is_reservior_present
        use mod_geometry , only : num_element
        implicit none
        integer :: i , k , r , n_bc_max
                   
                       
        !! basic model: every translation dof in bulk is described by one state function.
        do i = 1  ,  num_element 
            if( element_property( i , 1 ) <  3 ) then   !region_id=1,2 : dam and foundation                                                   
                call  basic_model_elastic_bulk( i , .false. )
            end if
        end do

        
        !! surface pressure force appears only in dynamic analysis, i.e. analysis_kind = 2
        !! unless drm freq domain
        if( .not. is_time_domian ) return ;
   
        if( .not. is_reservior_present ) return ;                       !! only when reservoir exists
        
        n_bc_max = size( element_property , 2 ) - 3
     
        if( n_bc_max .eq. 0 ) return                                    ! if there is no boundary condition.
     
        do  i = 1  ,   num_element 
            if( element_property( i , 1 ) /= 1 ) cycle                  ! region_id = 1 : dam  
          
            do  k = 1 , n_bc_max
                r = element_property( i , 3 + k )
                if( r .eq. 0 ) cycle
               
                if( boundary_load( r )%k_type /= 3 ) then               ! if(k_type not.eq. 3 ) surface pressure force
                    call  elastic_surface_pressure_force( i , r )
                end if
                                                                        
            end do ! over k
        end do ! over i
        
        
    end subroutine fem_for_elastic_bulk
  
  
  
   

!!===========================================================================!!
!!
!!  submodule subroutine: fem_for_acoustic_bulk
!!
!!===========================================================================!!

  subroutine fem_for_acoustic_bulk()
       
        use mod_physics , only  : element_property
        use mod_geometry , only : num_element
        implicit none
        integer :: i 
            
       !! basic model: every pressure dof in bulk is described by one state function.
        do i = 1  , num_element 
            if( element_property( i , 1 ) .eq. 3 ) then       ! region_id = 3 : reservior                                                     
                call  basic_model_acoustic_bulk( i )
            end if
          
        end do

       
  end subroutine fem_for_acoustic_bulk
  
  
  


!!===========================================================================!!
!!
!!  submodule subroutine: fem_for_dam_reservior_boundary
!!  basic model: every pressure and translational dof in dam-reservior boundary 
!!               are described by one state function.
!!===========================================================================!!

    subroutine fem_for_dam_reservior_boundary()
       
        use mod_physics , only  : element_property , boundary_load
        use mod_geometry , only : num_element
        implicit none
        integer :: i , k , r , n_bc_max
       
        n_bc_max = size( element_property , 2 ) - 3
     
        if( n_bc_max .eq. 0 ) return                                     ! if there is no boundary condition.
     
        do  i = 1  , num_element 
            if( element_property( i , 1 ) /= 3 ) cycle                    ! region_id = 3 : reservior  
          
            do  k = 1 , n_bc_max
                r = element_property( i , 3 + k )
                if( r .eq. 0 ) cycle
               
                if( boundary_load( r )%k_type .eq. 3 ) then              ! dam-reservior boundary
                    call  basic_model_dam_reservior_boundary( i , r )
                end if
                                                                        
            end do ! over k
        end do ! over i

       
    end subroutine fem_for_dam_reservior_boundary
  
  
  



!!===========================================================================!!
!!
!!  submodule subroutine: fem_for_reservior_foundation_boundary
!!  basic model: every pressure dof in this boundary 
!!               is described by one state function.
!!
!!===========================================================================!!

  subroutine fem_for_reservior_foundation_boundary()
       
       use mod_physics , only  : element_property , boundary_load
       use mod_geometry , only : num_element
       implicit none
       integer :: i , k , r , n_bc_max
       
       n_bc_max = size( element_property , 2 ) - 3
     
       if( n_bc_max .eq. 0 ) return                                     ! if there is no boundary condition.
     
       do i = 1 , num_element 
          if( element_property( i , 1 ) /= 3 ) cycle                    ! region_id = 3 : reservior  
          
          do k = 1 , n_bc_max
               r = element_property( i , 3 + k )
               if( r .eq. 0 ) cycle
               
               if( boundary_load( r )%k_type .eq. 2 ) then              ! reservior_foundation_boundary
                   call  basic_model_reservior_foundation_boundary( i , r )
               end if
                                                                        
          end do ! over k

       end do ! over i

       
  end subroutine fem_for_reservior_foundation_boundary
  
  
  
  

!!===========================================================================!!
!!
!!  submodule subroutine: fem_for_foundation_truncation_boundary:
!!
!!  basic model: every translational dof in foundation_truncation_boundary
!!               is described by one state function.
!!
!!  
!!  HW model: every translational dof in foundation_truncation_boundary
!!               is described O(N,M) HW model.
!!
!!===========================================================================!!

  subroutine fem_for_foundation_truncation_boundary()
       
       use mod_physics , only  : element_property , boundary_load, the_model_description
       use mod_geometry , only : num_element
       implicit none
       integer :: i , k , r , n_bc_max
       
       n_bc_max = size( element_property , 2 ) - 3
     
       if( n_bc_max .eq. 0 ) return                                     ! if there is no boundary condition.
     
        
       do i =  1  ,  num_element 
          if( element_property( i , 1 ) /= 2 ) cycle                    ! region_id =  2 : foundation
          
          do k = 1 , n_bc_max
               r = element_property( i , 3 + k )
               if( r .eq. 0 ) cycle
               
               
               !! condition of HW model comes here.
                if( boundary_load( r )%k_type .eq. 3 ) then              !  foundation-infinity boundary
                   
                    if( the_model_description%foundation_inf_boundary_model_id .eq. 0 ) then   ! Lysmer BC
                   
                       call  basic_model_foundation_truncation_boundary( i , r )
                       
                    elseif( the_model_description%foundation_inf_boundary_model_id .eq. 1 ) then ! HW BC
                    
                       call  HW_model_foundation_truncation_boundary( i , r )
                     

                       
                    end if
                   
                end if
                                                                        
          end do ! over k

       end do ! over i

       
  end subroutine fem_for_foundation_truncation_boundary
  
   
  
  

!!===========================================================================!!
!!
!!  submodule subroutine: fem_for_reservior_truncation_boundary:
!!
!!  basic model: every pressure dof in reservior_truncation_boundary
!!               is described by one state function.
!!
!!  
!!  HW model: every pressure dof in reservior_truncation_boundary
!!               is described O(N,M) HW model.
!!
!!  
!!  
!!  GN model: every pressure dof in reservior_truncation_boundary
!!               is described O(N,M) GN model.
!!===========================================================================!!

  subroutine fem_for_reservior_truncation_boundary()
       
       use mod_physics , only  : element_property , boundary_load , the_model_description
       use mod_geometry , only : num_element
       implicit none
       integer :: i , k , r , n_bc_max
       
       n_bc_max = size( element_property , 2 ) - 3
     
       if( n_bc_max .eq. 0 ) return                                     ! if there is no boundary condition.
       
       
     
       do i = 1 , num_element 
          if( element_property( i , 1 ) /= 3 ) cycle                    ! region_id =  3 : reservior
          
          do k = 1 , n_bc_max
               r = element_property( i , 3 + k )
               if( r .eq. 0 ) cycle
               
               if( boundary_load( r )%k_type .eq. 1 ) then               ! reservior-infinity boundary
                   
                   !! check BC model at this point.
                   if( the_model_description%reservior_inf_boundary_model_id .eq. 0 ) then   ! Somerfel BC
                   
                       call  basic_model_reservior_truncation_boundary( i , r )
                       
                   elseif( the_model_description%reservior_inf_boundary_model_id .eq. 1 ) then ! HW BC
                   
                       call  HW_model_reservior_truncation_boundary( i , r )
                       
                   elseif( the_model_description%reservior_inf_boundary_model_id .eq. 2 ) then ! GN BC
                    
                       call  GN_model_reservior_truncation_boundary( i , r )
                       
                   end if
                   
               end if
                                                                        
          end do ! over k

       end do ! over i

       
  end subroutine fem_for_reservior_truncation_boundary
  



!!====================================================================!!
!!
!!  find_res_found_inf_common_node and set alpha_reflection_coefficient
!!  Such information is required in higher order model in reservior.
!!====================================================================!!
  
  subroutine  find_res_found_inf_common_node()
     use mod_physics  , only : element_property , boundary_load , alpha_reflection_coefficient
     use mod_geometry , only : num_element, degree_of_freedom, element_matrix
     use mod_utils    , only : get_nodes_on_element_face 
     implicit none
     
     integer :: i , k , r , n_bc_max , n_face, n_bc_nodes , j , n
     
     res_inf_found_node_index = -1
     alpha_reflection_coefficient = 0.0_r_kind
     
     
     n_bc_max = size( element_property , 2 ) - 3
     
     if( n_bc_max .eq. 0 ) return
        
     do i =  1 , num_element 
        
        if( element_property( i , 1 ) /= 3 ) cycle                     ! element not in region_id =  3 : reservior
                                                                                                                                   
        do k = 1 , n_bc_max
           r = element_property( i , 3 + k )
               
           if( r .eq. 0 ) cycle
           if( boundary_load( r )%k_type /= 2 ) cycle                   ! otherwise it is not res-found boundary
           
           alpha_reflection_coefficient = boundary_load( r )%strength  ! used as alpha in res-found boundary
           n_face = boundary_load( r )%face_id
           indexed_set = 0
           n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
           do j = 1 , n_bc_nodes
              n = element_matrix( i , indexed_set( j ) ) 
              if( n .eq. 0 ) cycle                                      ! if node is absent.
              if( degree_of_freedom(n, 1 ) .eq. 3 ) then                ! class_id = 3 : node in res-inf boundary
                  res_inf_found_node_index = n 
                  exit
              end if
           end do ! over j
    
           if( res_inf_found_node_index > 0 ) exit                      !! exit k-loop                                                    
        end do ! over k
        
        if( res_inf_found_node_index > 0 ) exit                        !! exit i-loop  
        
     end do ! over i  
     
     if( res_inf_found_node_index < 0 ) then
         stop 'The common node of res-found-inf boundaries cannot be found in higher order model' 
     end if
     
  
  end subroutine find_res_found_inf_common_node    
        
        
        
  
  
!!====================================================================!!
!!
!!  find corner nodes in foundation truncation boundary
!!  The method assumes that foundation geometry is a rectangle.
!!
!!====================================================================!!
  
  subroutine  find_found_inf_corner_node()
     use mod_geometry , only : num_node , node , degree_of_freedom
     implicit none
     
     integer :: i , j  
     real( kind = r_kind ) :: y_min , y_max , x_min , x_max             ! enclosing rectangle
     real( kind = r_kind ) :: eps_y_min , eps_y_max , eps_x_min , eps_x_max , dy , dx
     
     !! initial value
     found_inf_corner_node(1:4) = -2 ;
     
     !! to avoid comparing floating point numbers, choose x_max , ..., y_min as mean values 
     !! of all nodes at foundation truncation boundary
     
     j = 0 ;  !! a counter
     
     !! find enclosing box
     do  i = 1 , num_node
         if( degree_of_freedom(i, 1) /=  2 ) cycle                      !! node is not in found-inf class
         if( j .eq. 0 ) then
             x_min = node(i,1) ;
             x_max = x_min ;
             y_min = node(i,2) ;
             x_max = y_min ;
             j =  1 ;
             cycle ;
         end if
         
         x_max = max( x_max , node(i,1) ) ;
         x_min = min( x_min , node(i,1) ) ;
         y_max = max( y_max , node(i,2) ) ;
         y_min = min( y_min , node(i,2) ) ;
     end do ! over i
     

     eps_x_min = epsilon( x_min ) ;
     eps_x_max = epsilon( x_max ) ;
     eps_y_min = epsilon( y_min ) ;
     eps_y_max = epsilon( y_max ) ;
     
     
     !! find corner nodes
     do  i = 1 , num_node
         if( degree_of_freedom(i, 1) /=  2 ) cycle                      !! node is not in found-inf class
         
         !! first node
         dx = node(i,1) - x_max ;
         dy = node(i,2) - y_min ;
         
         !! first corner node
         if( abs( dx ) + abs( dy) < 10.0_r_kind * ( eps_x_max  +  eps_y_min ) ) then
             found_inf_corner_node(1) = i ;
         end if 
         
         !! second node
         dy = node(i,2) - y_max ;
         if( abs( dx ) + abs( dy) < 10.0_r_kind * ( eps_x_max  +  eps_y_max ) ) then
             found_inf_corner_node(2) = i ;
         end if 
         
         !! third node
         dx = node(i,1) - x_min ;
         if( abs( dx ) + abs( dy) < 10.0_r_kind * ( eps_x_min  +  eps_y_max ) ) then
             found_inf_corner_node(3) = i ;
         end if  
         
         !! fourth node
         dy = node(i,2) - y_min ;
         if( abs( dx ) + abs( dy) < 10.0_r_kind * ( eps_x_min  +  eps_y_min ) ) then
             found_inf_corner_node(4) = i ;
         end if 
          
     end do ! over i
     
     
     do  i = 1 , 4
         if( found_inf_corner_node(i) < 0 ) then
             stop 'Failed in finding corner nodes in foundation truncation boundary.' 
         end if
     end do
     
  
  end subroutine find_found_inf_corner_node   
           
           
               
      

!!=====================================================================!!
!!
!!  the fem model constructor
!!
!!=====================================================================!!
     
  
   subroutine  construct_fem_model()
   
     use mod_geometry , only : num_dim
     use mod_physics  , only : r => num_max_variables_in_eqns 
     use mod_physics  , only : the_model_description
     implicit none   
      
     integer :: state 
      
     call destruct_fem_model()                                         ! check if the element arrays are null
                       
     
     !! set element type  
     if( num_dim .eq. 2 ) then
        the_element  => a_Q8_element                                    ! pointer points to 2D subclass                          
     else
        the_element  => a_brick_element                                 ! pointer points to 3D subclass 
     end if
     
     
     call the_element%construct_element()                               ! element constructor   
     
     !! following variables are common between different fem models

     allocate( elemental_mass_matrix( r  , r ) , stat = state )              
     if( state /= 0 ) stop 'Failed in allocating elemental_mass_matrix in fem module ' 
     
     allocate( elemental_stiffness_matrix( r  , r ) , stat = state )
     if( state /= 0 ) stop 'Failed in allocating elemental_stiffness_matrix in fem module ' 
     
     allocate( elemental_damping_matrix( r  , r ) , stat = state )
     if( state /= 0 ) stop 'Failed in allocating elemental_damping_matrix in fem module ' 
     
     allocate( elemental_force_vector( r  , num_dim ) , stat = state )
     if( state /= 0 ) stop 'Failed in allocating elemental_force_vector in fem module '
     
     allocate( row_variable_set( r ) , stat = state )
     if( state /= 0 ) stop 'Failed in allocating row_variable_set in fem module ' 
     
     allocate( col_variable_set( r ) , stat = state )
     if( state /= 0 ) stop 'Failed in allocating col_variable_set in fem module ' 
     
     
     !! the variables used in stiffness of elastic media: common in all models
     
     if( num_dim .eq. 2 ) then
         allocate( dN( 2 , 3 ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating dN matrix in fem module ' 
         
         allocate( dN_rigidity( 2 , 3 ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating dN_rigidity in fem  module ' 
     else
         allocate( dN( 3 , 6 ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating dN matrix in fem module ' 
         
         allocate( dN_rigidity( 3 , 6 ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating dN_rigidity in fem module ' 
     end if
     
     allocate( dN_D_dN( num_dim , num_dim ) , stat = state )                     
     if( state /= 0 ) stop 'Failed in allocating dN_D_dN matrix in fem module ' 
     
     
      ! the indexed set is max node in element
     allocate( indexed_set( the_element%num_node_in_element ) , stat = state )                     
     if( state /= 0 ) stop 'Failed in allocating indexed_set vector in fem module ' 
          
     
     !! used when foundation truncation boundary is Lysmer
     
     if( the_model_description%foundation_inf_boundary_model_id  .eq. 0 ) then    !! Lysmer BC
         allocate( rotation( num_dim , num_dim ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating rotation in fem module ' 
     
         allocate( C_Lysmer( num_dim , num_dim ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating C_Lysmer in fem module ' 
     end if
     
     
              
     !! used if reservior truncation boundary is GN or HW model
     
     if( the_model_description%reservior_inf_boundary_model_id  >  0 ) then    !! GN or HW models
     
         allocate( boundary_shape_matrix( r , r ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating boundary_shape_matrix in fem module ' 
     
         allocate( boundary_der_shape_matrix( r , r ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating boundary_der_shape_matrix in fem module '
         
         !! used for both HW and GN res-inf
         allocate( alpha_coeffs( r ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating alpha_coeffs in fem module ' 
         
     end if
     
     !! used if foundation truncation boundary model is HW 
     
     if( the_model_description%foundation_inf_boundary_model_id  >  0 ) then    !!  HW model
     
         allocate( alpha_coeffs( r ) , stat = state )                     
         if( state /= 0 ) stop 'Failed in allocating alpha_coeffs in fem module ' 
        
     end if
     



   end subroutine construct_fem_model
   

!!=====================================================================!!
!!
!!  the fem model destructor
!!
!!=====================================================================!!
   
   
   
   subroutine destruct_fem_model ()
   
     implicit none
     integer :: state  , state1 = 0
     
     
     if( associated( the_element ) ) then
         call  the_element%destruct_element()
         nullify( the_element )
     end if
     
     
     if( allocated(  elemental_mass_matrix  ) ) then
         deallocate( elemental_mass_matrix , stat = state )
         state1 = state1 + state
     end if
     
     if( allocated(  elemental_stiffness_matrix)) then
         deallocate( elemental_stiffness_matrix , stat = state )
         state1 = state1 + state
     end if

     if( allocated(  elemental_damping_matrix ) ) then
         deallocate( elemental_damping_matrix , stat = state )
         state1 = state1 + state
     end if

     if( allocated(  elemental_force_vector ) ) then
         deallocate( elemental_force_vector , stat = state )
         state1 = state1 + state
     end if

     
     if( allocated(  row_variable_set ) )  then
         deallocate( row_variable_set , stat = state )
         state1 = state1 + state
     end if

     if( allocated(  col_variable_set ) )  then
         deallocate( col_variable_set , stat = state )
         state1 = state1 + state
     end if

     if( allocated(  indexed_set ) ) then
         deallocate( indexed_set , stat = state )
         state1 = state1 + state
     end if

          
      !! destrcut matrices used in calculating stiffness of elastic media
                        
      if( allocated(  dN ) )  then
          deallocate(  dN , stat = state )
          state1 = state1 + state
      end if
     
      if( allocated(  dN_rigidity ) ) then
          deallocate( dN_rigidity , stat = state )
          state1 = state1 + state
      end if

      if( allocated(  dN_D_dN ) ) then
          deallocate( dN_D_dN , stat = state )
          state1 = state1 + state
      end if
   
      if( allocated( indexed_set ) ) then
          deallocate( indexed_set , stat = state)
          state1 = state1 + state
      end if 

      
      
      !! free matrices used for Lysmer BC
      
      if( allocated( rotation ) )  then
          deallocate( rotation , stat = state)
          state1 = state1 + state 
      end if

      
      if( allocated( C_Lysmer ) )  then
          deallocate( C_Lysmer , stat = state)
          state1 = state1 + state
      end if 

      
      !! free matrices used for reservior HW BC
      
      if( allocated( boundary_shape_matrix ) )  then
          deallocate( boundary_shape_matrix , stat = state)
          state1 = state1 + state 
      end if

      
      if( allocated( boundary_der_shape_matrix ) ) then
          deallocate( boundary_der_shape_matrix , stat = state)
          state1 = state1 + state 
      end if

      
      if( allocated( alpha_coeffs ) ) then
          deallocate( alpha_coeffs , stat = state)
          state1 = state1 + state 
      end if
      
      
      
      if( state1 /= 0 ) stop ' Error in deallocating arrays in fem module'
     
   
   end subroutine destruct_fem_model 
   
   
   
     
    
 
     
     
              
   

   
   
!!=====================================================================!!
!!
!!  set Gauss points and their weight
!!
!!====================================================================!!
  
  subroutine  set_gauss_points( n_point )
     
     implicit  none
     integer , intent( in ) :: n_point
     
     if( n_point .eq. current_num_Gauss_point ) return
     
     current_num_Gauss_point = n_point
     
     if( n_point .eq. 2 ) then
     
         gauss_points(1) = -1.0_r_kind / sqrt( 3.0_r_kind) 
         gauss_points(2) = -gauss_points(1)
         weight_of_gauss_point(1:2) = 1.0_r_kind
         
     else if( n_point .eq. 3 ) then
        
        gauss_points(1) =  -sqrt( 3.0_r_kind / 5.0_r_kind )
        gauss_points(2) =   0.0_r_kind
        gauss_points(3) = - gauss_points(1)
        weight_of_gauss_point(1)   = 5.0_r_kind / 9.0_r_kind
        weight_of_gauss_point(2)   = 8.0_r_kind / 9.0_r_kind
        weight_of_gauss_point(3)   = 5.0_r_kind / 9.0_r_kind
     
     else if( n_point .eq. 4 ) then
        
        gauss_points(1) =  -sqrt( 3.0_r_kind / 7.0_r_kind + ( 2.0_r_kind / 7.0_r_kind ) * &
                         &  sqrt( 6.0_r_kind / 5.0_r_kind ) )
                            
        gauss_points(2) =  -sqrt( 3.0_r_kind / 7.0_r_kind - ( 2.0_r_kind / 7.0_r_kind ) * &
                         &  sqrt( 6.0_r_kind / 5.0_r_kind ) )
         
        gauss_points(3) =  -gauss_points(2)
                            
        gauss_points(4) =  -gauss_points(1)
         
        weight_of_gauss_point(1) = ( 18.0_r_kind - sqrt( 30.0_r_kind) ) / 36.0_r_kind
        weight_of_gauss_point(2) = ( 18.0_r_kind + sqrt( 30.0_r_kind) ) / 36.0_r_kind
        weight_of_gauss_point(3) = weight_of_gauss_point(2)
        weight_of_gauss_point(4) = weight_of_gauss_point(1)
        
     else 
         stop 'Number of Gauss points is more than 4 in fem module. ' 
     end if
     
  end subroutine set_gauss_points
  
  
!!======================================================================
!!
!!  set rigidity matrix
!!
!!======================================================================

  subroutine  set_rigidity_matrix( elm_id , reg_id)
  
     use mod_physics  , only : element_property , material , rigidity_matrix , element_type
     use mod_geometry , only : num_dim
     
     implicit  none
     integer , intent( in ) :: elm_id , reg_id
     logical  :: is_plane_stress 
     real( kind = r_kind ) :: E , nu  , a
     
     if( reg_id .eq. 3 ) return                                             ! if region_id = 3 is reservior, do nothing.
     if( element_property( elm_id , 2 ) .eq. current_material_id  ) return  ! if mat_id are the same, do nothing
    
         
     current_material_id = element_property( elm_id , 2 )
     
     rigidity_matrix = 0.0_r_kind
     
     E  = material( element_property( elm_id , 2 ) )%elastic_module
     nu = material( element_property( elm_id , 2 ) )%poisson_ratio
     
     if( num_dim .eq. 2 ) then
         if( element_type( reg_id ) .eq. 5 ) then                       ! plane strain 
             is_plane_stress = .false.
         elseif( element_type( reg_id ) .eq. 6 ) then                   ! plane stress
             is_plane_stress = .true.
         else
             stop 'The element type has not defined properly, plane stress or plane strain'
         end if
     end if
     
     if( num_dim .eq. 2 ) then
         if( is_plane_stress )  then
             rigidity_matrix( 1 , 1 ) = E / ( 1.0_r_kind - nu ** 2 )
             rigidity_matrix( 2 , 2 ) = rigidity_matrix( 1 , 1 )
         
             rigidity_matrix( 1 , 2 ) = nu * rigidity_matrix( 1 , 1 )
             rigidity_matrix( 2 , 1 ) = rigidity_matrix( 1 , 2 )
         
             rigidity_matrix( 3 , 3 ) = 0.50_r_kind * E / ( 1.0_r_kind + nu )
         
         else
             a = E / ( ( 1.0_r_kind + nu  ) * ( 1.0_r_kind -  2.0_r_kind * nu  ) )
         
             rigidity_matrix( 1 , 1 ) = a * ( 1.0_r_kind - nu  ) 
             rigidity_matrix( 2 , 2 ) = rigidity_matrix( 1 , 1 )
         
             rigidity_matrix( 1 , 2 ) = a * nu
             rigidity_matrix( 2 , 1 ) = rigidity_matrix( 1 , 2 )
         
             rigidity_matrix( 3 , 3 ) = 0.50_r_kind * E / ( 1.0_r_kind + nu  )
         
         end if
         
         return
     end if
     
     !! num_dim = 3  in the following.

     a = E / ( ( 1.0_r_kind + nu  ) * ( 1.0_r_kind -  2.0_r_kind * nu  ) )
         
     rigidity_matrix( 1 , 1 ) = ( 1.0_r_kind - nu ) * a
     rigidity_matrix( 2 , 2 ) = rigidity_matrix( 1 , 1 )
     rigidity_matrix( 3 , 3 ) = rigidity_matrix( 1 , 1 )
         
     rigidity_matrix( 1 , 2 ) = nu * a
     rigidity_matrix( 1 , 3 ) = rigidity_matrix( 1 , 2 )
     rigidity_matrix( 2 , 1 ) = rigidity_matrix( 1 , 2 )
     rigidity_matrix( 2 , 3 ) = rigidity_matrix( 1 , 2 )
     rigidity_matrix( 3 , 1 ) = rigidity_matrix( 1 , 2 )
     rigidity_matrix( 3 , 2 ) = rigidity_matrix( 1 , 2 )
         
     rigidity_matrix( 4 , 4 ) = E / ( 2.0_r_kind * ( 1.0_r_kind + nu  ) )
     rigidity_matrix( 5 , 5 ) = rigidity_matrix( 4 , 4 )
     rigidity_matrix( 6 , 6 ) = rigidity_matrix( 4 , 4 )

     
  
  end subroutine set_rigidity_matrix
 
  
  
  
  
  
  end module mod_fem
