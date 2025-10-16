!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!

!!==============================================================================!!
!!
!!
!!  This submodule evaluates FEM calculation in dam-reservior boundary.
!! 
!!
!!===============================================================================!!

  
  submodule ( mod_fem ) smod_dam_res_boundary
  
    use  mod_utils , only : r_kind
    implicit none
    
  
    contains
    
    
  
  
!!====================================================================!!
!!
!!  The fem analysis in dam-reservior boundary in basic model,
!!  i.e. the pressure and every component of translational dofs have only one state function.
!!  The equation is 
!!   G \ddot(p) + H p = R_III
!!   where 
!!       R_III = - int( N n^T Q^T )( \ddot(r) + J a ) = - B_{III} \ddot(r) - B_{III} J a_g
!!   
!!  Here the matrix B_{III} plays the role of mass matrix for translation dof of nodes on this boundary.
!!  It has the form
!!  
!!  This matrix has the following shape in 2D problems
!!
!!   B_{III} = [ N_1 n_1 N_1 , N_1 n_2 N_1  ,    N_1 n_1 N_2 , N_1 n_2 N_2 , ..., N_1 n_1 N_k , N_1 n_2 N_k ... ;
!!               N_2 n_1 N_1 , N_2 n_2 N_1  ,    N_2 n_1 N_2 , N_2 n_2 N_2 , ..., N_2 n_1 N_k , N_2 n_2 N_k ... ;
!!               N_3 n_1 N_1 , N_3 n_2 N_1  ,    N_3 n_1 N_2 , N_3 n_2 N_2 , ..., N_3 n_1 N_k , N_3 n_2 N_k ... ;
!!            ... ]
!!   
!!  In 3D it is of the form:
!
!!  B_{III} = [ N_1 n_1 N_1 , N_1 n_2 N_1  ,  N_1 n_3 N_1   ,      N_1 n_1 N_2 , N_1 n_2 N_2 , N_1 n_3 N_2 , ... ;
!!              N_2 n_1 N_1 , N_2 n_2 N_1  ,  N_2 n_3 N_1    ,     N_2 n_1 N_2 , N_2 n_2 N_2 , N_2 n_3 N_2 , ... ;
!!              N_3 n_1 N_1 , N_3 n_2 N_1  ,  N_3 n_3 N_1    ,     N_3 n_1 N_2 , N_3 n_2 N_2 , N_3 n_3 N_2 , ... ;
!!            ... ] 
!!
!! 
!! Then 
!!     B_{III} * J = [ N_1 n_1 (N_1 + N_2 + ... + N_k)  ,  N_1 n_2 (N_1 + N_2 + ... + N_k)  ,  N_1 n_3 (N_1 + N_2 + ... + N_k) ;
!!                     N_2 n_1 (N_1 + N_2 + ... + N_k)  ,  N_2 n_2 (N_1 + N_2 + ... + N_k)  ,  N_2 n_3 (N_1 + N_2 + ... + N_k) ;
!!                     N_3 n_1 (N_1 + N_2 + ... + N_k)  ,  N_3 n_2 (N_1 + N_2 + ... + N_k)  ,  N_3 n_3 (N_1 + N_2 + ... + N_k) ;
!!                    ...]
!!
!! That can also be written as
!! 
!!   S = sum( N_I ) = N_1 + N_2 + ...+ N_k : I runs over facial nodes with pressure dof.
!!   
!!   B_{III} * J = [ N_1 n_1 S  ,  N_1 n_2 S  ,  N_1 n_3 S ;
!!                   N_2 n_1 S  ,  N_2 n_2 S  ,  N_2 n_3 S ;
!!                   N_3 n_1 S  ,  N_3 n_2 S  ,  N_3 n_3 S ;
!!                    ...]
!! 
!!
!!  Once this matrix is found, one should also add F = B^T * p in the elastic equation
!!
!!      M \ddot( u ) + K u = B^T * p
!!
!!  The two sets of equations become
!!
!!  G \ddot(p) + B_{III} \ddot(r) + H p = F_III    : here B is added to global mass matrix
!!  M \ddot( u ) + K u - B_{III}^T * p = 0         : here B is added to global stiffness matrix
!!
!!  where F_{III} = - B_{III} * J a_g
!!
!!  In Q8 element, B_III is  8 * 16 matrix where its domain is 8 * 2 = 16 translational dof of nodes and
!!  its image is 8 pressure dof of nodes.
!!  
!!  In contrast, B^T is 16*8 matrix where its domain is 8 pressure dof and output ( image ) is 8 * 2 = 16
!!  where it is multiplied by pressure and gives translatonal dof.
!!
!!====================================================================!!
   
  module subroutine basic_model_dam_reservior_boundary( elem_id , bc_id )
   
    use mod_physics  , only : element_property , boundary_load , add_to_global_matrices , add_to_resultant_force
    use mod_physics  , only : drm_layer_ndyn
    use mod_geometry , only : num_dim , degree_of_freedom , element_matrix
    use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face 
    use  mod_fem     , only : set_gauss_points 
    implicit none
    
    integer , intent( in ) :: elem_id , bc_id
    integer :: n_bc_nodes , n_face , direction_value(2) , n_gauss_point, num_col , num_row
    integer :: k , n , r, s , d, q1 
    real( kind = r_kind ) :: fixed_value , d_volume , temp
    
     !!----------------------------------------------------------------!
     !! set initial values for damping matrix, force vector and variabl_set
     
     elemental_mass_matrix      = 0.0_r_kind
     
     elemental_stiffness_matrix = 0.0_r_kind
     
     elemental_force_vector     = 0.0_r_kind
     
     indexed_set = 0
     
   
     !!----------------------------------------------------------------!
     !! set gauss points
     
     n_gauss_point = element_property( elem_id , 3 )
      
     call  set_gauss_points( n_gauss_point )
     
     n_face = boundary_load( bc_id )%face_id
     n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
     
     call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
     fixed_value = real( direction_value(2) , kind = r_kind )           !! either xi or eta or zeta is fixed at +1 or -1
     
   
     
     !!----------------------------------------------------------------!
     !! set suitable variable set for mass matrix in fluid equation.
     
     row_variable_set = 0
     col_variable_set = 0
         
     do k = 1 , n_bc_nodes
            
        n = element_matrix( elem_id , indexed_set(k) ) 
            
        if( n .eq. 0 ) cycle                                            ! if node is absent or no dof
        if( degree_of_freedom( n , 1 ) .eq. 0 ) cycle                   ! node has no dof
        
        q1 = ( indexed_set(k) - 1 ) * num_dim
        
        !! row var set is pressure dof of nodes
        row_variable_set( indexed_set(k) ) =  degree_of_freedom( n , 2 * num_dim + 2 )                     
        
        !! col var set is translational dof of nodes
        do  d =  1 , num_dim
            col_variable_set( q1 + d ) = degree_of_freedom( n , 2 * d )
        end do          
     end do ! over k
     
     num_row = the_element%num_node_in_element
     num_col =  num_dim * the_element%num_node_in_element
    
   
     
     !!----------------------------------------------------------------!!
     !! form damping matrix and force vector
     
     if( num_dim .eq. 2 ) then                                          ! there is only one direction to do integration
                 
         do r = 1 ,  n_gauss_point 
              
            if( direction_value(1) .eq. 1 ) then                        !! find shape function on face
               
                call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r)  )      ! xi = fixed
                d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) ) 
                d_volume = d_volume  * weight_of_gauss_point(r)  
            else
               
                call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value  )      ! eta  = fixed
                d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value ) 
                d_volume =  d_volume * weight_of_gauss_point(r)   
            end if 
            
            
            call the_element%get_normal_at_face( elem_id , direction_value(1)  ,  &                      ! get outward normal
                                              &  indexed_set(1:n_bc_nodes) , is_outward = .true. )     
            
            !! calculate B_{III} and store in mass matrix
            
            do  k = 1 , n_bc_nodes   
                
                do n = 1 , n_bc_nodes 
                   
                   q1 = ( indexed_set( n ) - 1 ) * num_dim
                   
                   do  d = 1 , num_dim
                           elemental_mass_matrix( indexed_set( k ) , q1 + d ) =  &
                        &  elemental_mass_matrix( indexed_set( k ) , q1 + d ) +  &
                        &  the_element%shape_vector( indexed_set(k) ) *  d_volume * &
                        &  the_element%shape_vector( indexed_set(n) ) * the_element%normal_at_face( d ) 
                   end do ! over d
                   
                          
                end do ! over n
            end do ! over k
            
            
            
            !! calculate force
            !! get sum of shape functions of facial nodes
            temp = 0.0_r_kind
            do  k = 1 , n_bc_nodes   
                temp = temp + the_element%shape_vector( indexed_set(k) )
            end do ! over k
                 
            
            !! force vector
            do k = 1 , n_bc_nodes
               n = indexed_set( k )
               
               do  d = 1 , num_dim
                   elemental_force_vector( n , d ) = elemental_force_vector( n , d ) - &
                &  the_element%shape_vector( n ) * temp * the_element%normal_at_face( d ) * d_volume
               end do
            end do ! over k         
            
                             
         end do ! over r
          
          
     else  !! 3D problem
      
          do r = 1 , n_gauss_point                                      
             do s = 1 ,  n_gauss_point
             
                if( direction_value(1) .eq. 1 ) then                    !! xi = fixed_value
               
      call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r) , gauss_points(s) ) 
      d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) , gauss_points(s) ) 
      d_volume = abs( d_volume ) * weight_of_gauss_point(r) * weight_of_gauss_point(s)
                    
                            
                elseif( direction_value(1) .eq. 2 ) then                !! eta = fixed_value
               
      call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value , gauss_points(s) )   
      d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value , gauss_points(s) ) 
      d_volume = abs( d_volume ) * weight_of_gauss_point(r) * weight_of_gauss_point(s) 
                             
                else                                                    !! zeta = fixed_value
                       
      call the_element%get_shape_function( elem_id , gauss_points(r) , gauss_points(s) , fixed_value  )
      d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , gauss_points(s) , fixed_value ) 
      d_volume = abs( d_volume ) * weight_of_gauss_point(r) * weight_of_gauss_point(s)
                                              
                end if 
                
                
      call the_element%get_normal_at_face( elem_id , direction_value(1)  ,  &                      ! get outward normal
                                              &  indexed_set(1:n_bc_nodes) , is_outward = .true. )     
            
                      
                !! calculate B_{III} and storing in mass matrix
            
                do  k = 1 , n_bc_nodes   
                
                    do n = 1 , n_bc_nodes 
                   
                       q1 = ( indexed_set( n ) - 1 ) * num_dim
                   
                       do  d = 1 , num_dim
                           elemental_mass_matrix( indexed_set( k ) , q1 + d ) =  &
                        &  elemental_mass_matrix( indexed_set( k ) , q1 + d ) +  &
                        &  the_element%shape_vector( indexed_set(k) ) *  d_volume * &
                        &  the_element%shape_vector( indexed_set(n) ) * the_element%normal_at_face( d ) 
                       end do ! over d
                   
                          
                    end do ! over n
                end do ! over k
            
                
            
                !! calculate force
                !! get sum of shape functions of facial nodes
                temp = 0.0_r_kind
                do  k = 1 , n_bc_nodes   
                    if( indexed_set( k ) .eq. 0 ) cycle                     ! node is absent or no dof 
                    temp = temp + the_element%shape_vector( indexed_set(k) )
                end do ! over k
            
            
            
                !! force vector
                do k = 1 , n_bc_nodes
                   if( indexed_set( k ) .eq. 0 ) cycle                      ! node is absent or no dof 
                   n = indexed_set( k )
               
                   do  d = 1 , num_dim
                       elemental_force_vector( n , d ) = elemental_force_vector( n , d ) - &
                    &  the_element%shape_vector( n ) * temp * the_element%normal_at_face( d ) * d_volume
                   end do
                end do ! over k
                
             end do ! over s
          end do ! over r
          
     end if
     

     
       
     !! assemble to global matrix
     call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                 & num_row , col_variable_set , num_col  )
       
    !! assemble to global force vector
    if( drm_layer_ndyn .eq. 0 ) then
        call  add_to_resultant_force( elemental_force_vector , row_variable_set  , num_row )
    end if  

     
     !! now form -B^T * p and add to stiffness matrix
     
     elemental_stiffness_matrix = 0.0_r_kind
     do  k = 1 , num_row
         do  n = 1 , num_col
             elemental_stiffness_matrix( n , k ) = - elemental_mass_matrix( k , n )  ! transpose      
         end do ! over n
     end do ! over k
     
     !! for this input (col_var) is prossure dof and output( row_var ) is translational dof  
     
     call  add_to_global_matrices( elemental_stiffness_matrix , 2 , col_variable_set , & 
                                 & num_col , row_variable_set , num_row  )
     
 
     
     
  end subroutine basic_model_dam_reservior_boundary
   
   
   
   
   
   
  
  end submodule smod_dam_res_boundary


