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
!!  This submodule evaluates FEM calculation in bulk of acoustic media, i.e. reservior.
!!  No boundary condition is taken into account here.
!!
!!===============================================================================!!

  
  submodule ( mod_fem ) smod_acoustic_bulk
  
    use  mod_utils , only : r_kind
    implicit none
    
  
    contains
    
    

!!====================================================================!!
!!
!!  the fem analysis in bulk of acoustic media in basic model, i.e. the pressure has only one state function.
!!  The equation is, here,
!!
!!  G \ddot{p} + H p = 0
!!
!!  Every pressure dof is described by 1 state function in this model.
!!
!!====================================================================!!

   
  module subroutine basic_model_acoustic_bulk( elem_id )
   
     use mod_physics  , only : element_property , material, add_to_global_matrices
     use mod_geometry , only : element_matrix, degree_of_freedom , num_dim 
     use  mod_fem     , only : set_gauss_points 
     implicit none
     integer , intent( in ) :: elem_id 
     
     integer :: k , n, r , s , t , n_gauss_point , d , num_var
     real( kind = r_kind ) ::  d_volume , c2 , rho
     
      
     !-----------------------------------------------------------------!
     !! set initial values for matrices
     
     elemental_mass_matrix      = 0.0_r_kind
      
     elemental_stiffness_matrix = 0.0_r_kind
      
     !----------------------------------------------------------------!
     !! get physical property of the element
      
     c2  = material( element_property( elem_id , 2) )%elastic_module ** 2           ! used as c (speed of sound)
     rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
     n_gauss_point = element_property( elem_id , 3 )
      
     call  set_gauss_points( n_gauss_point )
     
    
     
     
     !---------------------------------------------------------------!
     !! set suitable variables set: every dof has 1 state function in this model
     
    
     row_variable_set = 0
     num_var = the_element%num_node_in_element 
   
     do k = 1 , the_element%num_node_in_element
      
        n = element_matrix( elem_id , k )
         
        if( n .eq. 0 ) cycle                                            ! if node is absent.
        if( degree_of_freedom( n , 1 ) .eq. 0 ) cycle                   ! node has no dof
        
        row_variable_set( k ) = degree_of_freedom( n , 2 * num_dim + 2 )
         
     end do ! over k
     
     col_variable_set = row_variable_set
      
     !!---------------------------------------------------------------!!
     !!  2D element : calculate elemental mass (G) and stiffness matrices ( H )
       
     if( num_dim .eq. 2 ) then
                 
         do r = 1 , n_gauss_point
            do s = 1 , n_gauss_point
             
               call the_element%get_shape_function( elem_id , gauss_points(r) , gauss_points(s)  )
                
               call the_element%get_derivative_shape_function( elem_id , gauss_points(r) , gauss_points(s) )
                
               d_volume = weight_of_gauss_point(r) * weight_of_gauss_point(s) *  &
                       &  abs( the_element%determinant_of_jacobian )
                 
                
               do  k = 1 , size( the_element%shape_vector )             !! do sum over all nodes                
                   do n = k , size( the_element%shape_vector )          !! elemental matrices are symmetric
                                                     
                        
                       !! mass matrix (G in fluid)
                       elemental_mass_matrix(k,n)  = elemental_mass_matrix(k,n)  +  &
                     & the_element%shape_vector(k) * the_element%shape_vector(n) * d_volume
                       
           
                       !! stiffness matrix ( H in fluid )
                        
                       do d = 1 , num_dim
                          elemental_stiffness_matrix( k  , n ) = elemental_stiffness_matrix( k  , n ) + &
                       &  the_element%derivative_shape_matrix( d , k ) *  &
                       &  the_element%derivative_shape_matrix( d , n ) * d_volume
                       end do
                        
                     
                   end do ! over n
               end do ! over k
                  
            end do ! over s
         end do ! over r
          
          
     else  !! 3D problem
      
         do r = 1 , n_gauss_point
            do s = 1 ,  n_gauss_point
               do t = 1 , n_gauss_point
                
                  call the_element%get_shape_function( elem_id , gauss_points(r) , &
                                       &  gauss_points(s) , gauss_points(t) )
                   
                  call the_element%get_derivative_shape_function(  elem_id , gauss_points(r) ,  &
                                                   &  gauss_points(s) , gauss_points(t) )
                   
                  d_volume = weight_of_gauss_point(r) * weight_of_gauss_point(s) *  &
                           & weight_of_gauss_point(t) * abs( the_element%determinant_of_jacobian )
                        
                    
                  !! do sum over all nodes
                
                  do  k = 1 , size( the_element%shape_vector )
                       
                      do n = k , size( the_element%shape_vector )       ! elemental matrices are symmetric                    
                       
                           !! mass matrix
                           elemental_mass_matrix(k,n) = elemental_mass_matrix(k,n)  +  &
                        &  the_element%shape_vector(k) * the_element%shape_vector(n) * d_volume
                       
                       
                           !! stiffness matrix
                           
                           do d = 1 , num_dim
                              elemental_stiffness_matrix( k  , n ) = elemental_stiffness_matrix( k  , n ) + &
                           &  the_element%derivative_shape_matrix( d , k ) * &
                           &  the_element%derivative_shape_matrix( d , n ) * d_volume
                           end do
                        
                      end do ! over n
                  end do ! over k
                
                
               end do ! over t
            end do ! over s
         end do !  over r

     end if
     
     
     !! symmetrize matrices
     
    do  k = 1 , size( the_element%shape_vector )
     
        elemental_mass_matrix( k , k )      = elemental_mass_matrix( k , k ) / ( rho * c2 )
        elemental_stiffness_matrix( k , k ) = elemental_stiffness_matrix( k , k ) / rho
         
        do n = k + 1 , size( the_element%shape_vector )
         
            elemental_mass_matrix( k , n )      = elemental_mass_matrix( k , n ) / ( rho * c2 )
            elemental_stiffness_matrix( k , n ) = elemental_stiffness_matrix( k , n ) / rho
            
            elemental_mass_matrix( n , k )      = elemental_mass_matrix( k , n )
            elemental_stiffness_matrix( n , k ) = elemental_stiffness_matrix( k , n )
            
        end do ! over n
     
    end do ! over k
     
     
      
     !! put G ( elemental mass matrix ) in global mass matrix
     
      
     call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                   & num_var , col_variable_set , num_var  )
       
    
     !! put H ( elemental stiffness matrix ) in global stiffness matrix                    
     call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                 & num_var , col_variable_set , num_var  ) 
     
     
   
   end subroutine basic_model_acoustic_bulk
   
   
   
   
   
   
   
   
  
  end submodule smod_acoustic_bulk

