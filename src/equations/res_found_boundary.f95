!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!
!
!!==============================================================================!!
!!
!!
!!  This submodule evaluates FEM calculation in reservior_foundation boundary.
!! 
!!
!!===============================================================================!!

  
  submodule ( mod_fem ) smod_res_found_boundary
  
    use  mod_utils , only : r_kind
    implicit none
    
  
    contains
    
    
  
  !!====================================================================!!
!!
!!  the fem analysis in reservior-foundation boundary in basic model
!!   G \ddot(p) + H p = R_II = (-q/rho ) * \int( N N^T) \dot(p)  - \int( N * n^T * Q^T ) J * a_g 
!!   or
!!   G \ddot(p) + C_II \dot(p) + H p = - \int( N * n^T * Q^T ) J * a_g
!! 
!!   where the damping matrix is C_II = + (q/rho ) * int( N N^T).
!!   
!!   and the force term is F_II = - B_II * a_g
!!   where 
!!        B_II = \int( N * n^T * Q^T ) = \int( N * n^T * N^T)
!!
!!   where N is shape function matrix of the boundary face of underlying element in fluid and
!!   and Q is the shape function matrix of the boundary face of nearby element in foundation.
!!   It has been shown that on boundary face Q = N
!!   that is, for example, N1 = Q4 where both referes to the same node in global nodal index set.
!!   Thus, here, the force term is written as
!!
!!     B_II = \int( N * n^T * N^T ) * J
!!   This matrix has the following shape in 2D problems
!!
!!   B_{II} = [ N_1 n_1 N_1 , N_1 n_2 N_1 ,    N_1 n_1 N_2 , N_1 n_2 N_2 , ..., N_1 n_1 N_k , N_1 n_2 N_k , ... ;
!!              N_2 n_1 N_1 , N_2 n_2 N_1 ,    N_2 n_1 N_2 , N_2 n_2 N_2 , ..., N_2 n_1 N_k , N_2 n_2 N_k , ... ;
!!              N_3 n_1 N_1 , N_3 n_2 N_1 ,    N_3 n_1 N_2 , N_3 n_2 N_2 , ..., N_3 n_1 N_k , N_3 n_2 N_k , ... ;
!!            ... ]
!!   
!!  In 3D it is of the form:
!
!!  B_{II} = [ N_1 n_1 N_1 , N_1 n_2 N_1  ,  N_1 n_3 N_1   ,      N_1 n_1 N_2 , N_1 n_2 N_2 , N_1 n_3 N_2 , ... ;
!!             N_2 n_1 N_1 , N_2 n_2 N_1  ,  N_2 n_3 N_1   ,      N_2 n_1 N_2 , N_2 n_2 N_2 , N_2 n_3 N_2 , ... ;
!!             N_3 n_1 N_1 , N_3 n_2 N_1  ,  N_3 n_3 N_1   ,      N_3 n_1 N_2 , N_3 n_2 N_2 , N_3 n_3 N_2 , ... ;
!!           ... ] 
!!
!! On the other hand J is, in 3D,
!!             J^T = [ 1 , 0 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , ... ;
!!                     0 , 1 , 0 , 0 , 1 , 0 , 0 , 1 , 0 , ....;
!!                     0 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , 1 , ....]
!! So
!!     B_{II} * J = [ N_1 n_1 (N_1 + N_2 + ... + N_k)  ,  N_1 n_2 (N_1 + N_2 + ... + N_k) , N_1 n_3 (N_1 + N_2 + ... + N_k) ;
!!                    N_2 n_1 (N_1 + N_2 + ... + N_k)  ,  N_2 n_2 (N_1 + N_2 + ... + N_k) , N_2 n_3 (N_1 + N_2 + ... + N_k) ;
!!                    N_3 n_1 (N_1 + N_2 + ... + N_k)  ,  N_3 n_2 (N_1 + N_2 + ... + N_k) , N_3 n_3 (N_1 + N_2 + ... + N_k) ;
!!                    ...]
!!
!! That can also be written as
!! 
!!   S = sum( N_I ) = N_1 + N_2 + ...+ N_k : I runs over facial nodes with pressure dof.
!!   
!!   B_{II} * J = [ N_1 n_1 S  , N_1 n_2 S , N_1 n_3 S ;
!!                  N_2 n_1 S  , N_2 n_2 S , N_2 n_3 S ;
!!                  N_3 n_1 S  , N_3 n_2 S , N_3 n_3 S ;
!!                 ...]
!!
!! The sum is automatically over facial nodes as, on the boundary face, the shape function is zero for other nodes.
!!
!!====================================================================!!
   
  subroutine basic_model_reservior_foundation_boundary( elem_id , bc_id )
   
    use mod_physics  , only : element_property , material , boundary_load , the_model_description 
    use mod_physics  , only : add_to_global_matrices , add_to_resultant_force
    use mod_geometry , only : num_dim , degree_of_freedom , element_matrix
    use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face 
    use  mod_fem     , only : set_gauss_points 
    implicit none
    
  
    integer , intent( in ) :: elem_id , bc_id
    integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point
    integer :: k , n , r, s , d 
    real( kind = r_kind ) :: rho , c, fixed_value , d_volume ,  q , alpha , temp
    
     !!----------------------------------------------------------------!
     !! set initial values for damping matrix, force vector and variabl_set
     
     elemental_damping_matrix = 0.0_r_kind
     
     elemental_force_vector   = 0.0_r_kind
     
     indexed_set = 0
     
     
     !!----------------------------------------------------------------!
     !! get physical property of the element
     
     alpha = boundary_load( bc_id )%strength                             ! used as alpha in res-found boundary
     c   = material( element_property( elem_id , 2) )%elastic_module     ! used as c (speed of sound)
     rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
   
     q = ( ( 1.0_r_kind - alpha ) / ( 1.0_r_kind + alpha ) ) / ( c * rho )
    
     
     !!----------------------------------------------------------------!
     !! set gauss points
     
     n_gauss_point = element_property( elem_id , 3 )
      
     call  set_gauss_points( n_gauss_point )
     
     n_face = boundary_load( bc_id )%face_id
     n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
     
     call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
     fixed_value = real( direction_value(2) , kind = r_kind )           !! either xi or eta or zeta is fixed at +1 or -1
     
   
     
     
     !!----------------------------------------------------------------!
     !! set suitable variable set for damping matrix.
     
     
     row_variable_set = 0
         
     do k = 1 , n_bc_nodes
            
        n = element_matrix( elem_id , indexed_set( k ) ) 
            
        if( n .eq. 0 ) cycle                                            ! if node is absent
        if( degree_of_freedom( n , 1 ) .eq. 0 ) cycle                   ! node has no dof
            
        row_variable_set( indexed_set( k ) ) = degree_of_freedom( n , 2 * num_dim + 2 ) 
        
     end do ! over k
     
     num_var =  the_element%num_node_in_element
     
     col_variable_set = row_variable_set 
     
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
            
                      
            !! calculate damping matrix
            
            do  k = 1 , n_bc_nodes   
                if( indexed_set( k ) .eq. 0 ) cycle                     ! node is absent or no dof 
                
                do n = 1 , n_bc_nodes 
                   if( indexed_set( n ) .eq. 0 ) cycle                  ! node is absent or no dof
                   
                   elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) =  &
                &  elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) +  &
                &  the_element%shape_vector( indexed_set(k) ) *  d_volume * &   
                &  the_element%shape_vector( indexed_set(n) ) *  q
                          
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
            
                      
                !! calculate damping matrix
            
                do  k = 1 , n_bc_nodes   
                    if( indexed_set( k ) .eq. 0 ) cycle                     ! node is absent or no dof 
                
                    do n = 1 , n_bc_nodes 
                       if( indexed_set( n ) .eq. 0 ) cycle                  ! node is absent or no dof
                   
                       elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) =  &
                    &  elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) +  &
                    &  the_element%shape_vector( indexed_set(k) ) *  d_volume * &   
                    &  the_element%shape_vector( indexed_set(n) ) *  q
                          
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
     call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                 & num_var , col_variable_set , num_var  ) 
           
     
     !! assemble to global force vector
     call  add_to_resultant_force( elemental_force_vector , row_variable_set , num_var )
     
    
   
   end subroutine basic_model_reservior_foundation_boundary
   
   
   
   
   
  
  end submodule smod_res_found_boundary



