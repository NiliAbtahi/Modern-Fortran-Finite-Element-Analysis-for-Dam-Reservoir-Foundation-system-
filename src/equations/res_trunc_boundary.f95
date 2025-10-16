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
!!  This submodule evaluates FEM calculation in reservior truncation boundary.
!! 
!!
!!===============================================================================!!

  
  submodule ( mod_fem ) smod_res_trunc_boundary
  
    use  mod_utils , only : r_kind
    implicit none
    
  
    contains
    
    

  
!!====================================================================!!
!!
!!  The fem analysis in reservior_truncation_boundary in basic model:                   Somerfeld model.
!!  i.e. every component of pressure dof has only one state function :                  Somerfeld model.
!!  The equation is 
!!   G \ddot(p) + H p = R_I = (-1/rho c ) * int( N N^T) \dot(p)
!!   or
!!   G \ddot(p) + C_I \dot(p) + H p = 0,
!! 
!!   where the damping matrix is C_I = + (1/rho c ) * int( N N^T).
!!
!!   The pressure is specified by 1 state function in linear model.
!! 
!!====================================================================!!
   
  module subroutine basic_model_reservior_truncation_boundary( elem_id , bc_id )  
   
    use mod_physics  , only : element_property , material , boundary_load , add_to_global_matrices
    use mod_geometry , only : num_dim , degree_of_freedom , element_matrix
    use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face
    use mod_fem      , only : set_gauss_points 
    implicit none


    integer , intent( in ) :: elem_id , bc_id
    integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point
    integer :: k , n , r, s
    real( kind = r_kind ) :: rho , c, fixed_value , coeff , d_volume
    
     !!----------------------------------------------------------------!
     !! set initial values for damping matrix
     
     elemental_damping_matrix = 0.0_r_kind
     indexed_set = 0
     
     
     !!----------------------------------------------------------------!
     !! get physical property of the element
     c   = material( element_property( elem_id , 2) )%elastic_module     ! used as c (speed of sound)
     rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
    
     coeff = 1.0_r_kind / ( rho * c )
     
     
     !!----------------------------------------------------------------!
     !! set gauss points
     
     n_gauss_point = element_property( elem_id , 3 )
      
     call  set_gauss_points( n_gauss_point )
     
     n_face = boundary_load( bc_id )%face_id
     n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
     
     call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
     fixed_value = real( direction_value(2) , kind = r_kind )           !! either xi or eta or zeta is fixed at +1 or -1
     
     
     !!----------------------------------------------------------------!
     !! set suitable variable set.
     
     row_variable_set = 0
     
     do k = 1 , n_bc_nodes
        n = element_matrix( elem_id , indexed_set(k) )
        if( n .eq. 0 )  cycle
        if( degree_of_freedom( n , 1 ) .eq. 0 ) cycle 
        row_variable_set( indexed_set(k) ) = degree_of_freedom( n , 2 * num_dim + 2 ) 
     end do
     

     num_var =  the_element%num_node_in_element
     col_variable_set = row_variable_set 
     
    
     
     !!----------------------------------------------------------------!!
     !! form damping matrix
     
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
            
            
           
            do  k = 1 , n_bc_nodes   
                
                do n = 1 , n_bc_nodes 
                   
                   elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) =  &
                &  elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) +  &
                &  the_element%shape_vector( indexed_set(k) ) *  d_volume * &   
                &  the_element%shape_vector( indexed_set(n) ) *  coeff
                          
                end do ! over n
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
                
                
               
                
                
                !! form damping matrix from Lysmer boundary condition.
                
                do  k = 1 , n_bc_nodes   
                
                    do n = 1 , n_bc_nodes 
                   
                       elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) =  &
                    &  elemental_damping_matrix( indexed_set( k ) , indexed_set( n ) ) +  &
                    &  the_element%shape_vector( indexed_set(k) ) *  d_volume * &   
                    &  the_element%shape_vector( indexed_set(n) ) *  coeff
                          
                    end do ! over n
                end do ! over k
                  
                
             end do ! over s
          end do ! over r
          
     end if
       
     !! assemble to global matrix
     call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                 & num_var , col_variable_set , num_var  ) 
      

  end subroutine basic_model_reservior_truncation_boundary
   
   
   
   


!!====================================================================!!
!!
!!  The fem analysis in reservior_truncation_boundary in GN model.
!!
!! 
!!====================================================================!!
   
 module subroutine GN_model_reservior_truncation_boundary( elem_id , bc_id )  
   
    use mod_physics  , only : element_property , material , boundary_load , alpha_reflection_coefficient
    use mod_geometry , only : num_dim , degree_of_freedom , element_matrix
    use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face
    use mod_physics  , only : the_model_description , add_to_resultant_force , add_to_global_matrices
    use mod_fem      , only : set_gauss_points
    implicit none
 
    integer , intent( in ) :: elem_id , bc_id
    integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point
    integer :: k , n , r, s , phi_id ,  special_node_index
    real( kind = r_kind ) :: rho , c, fixed_value , d_volume , q 
    
    
     
    !!----------------------------------------------------------------!
    !! step 1. get physical property of the element.
    
    c   = material( element_property( elem_id , 2) )%elastic_module     ! used as c (speed of sound)
    rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
     
    q = ( 1.0_r_kind - alpha_reflection_coefficient ) / ( 1.0_r_kind + alpha_reflection_coefficient ) / c 
     
     
!!======================================================================================================!!
!!
!!  Step 2. Form boundary_shape_matrix and boundary_der_shape_matrix.
!!
!!======================================================================================================!!    
     
     !!---------------------------------------------------------------!!
     !! 2.1. get local node index of boundary nodes
     
     indexed_set = 0
     
     n_face = boundary_load( bc_id )%face_id
     n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
   
     call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
     fixed_value = real( direction_value(2) , kind = r_kind )           !! either xi or eta or zeta is fixed at +1 or -1
     
     
      
     !!---------------------------------------------------------------!!
     !! 2.2.  calcuate local matrices for the three facial nodes: 3 , 4 , 7 in local index for Q8 element
     
     boundary_shape_matrix     = 0.0_r_kind
     boundary_der_shape_matrix = 0.0_r_kind
     
     
     n_gauss_point = element_property( elem_id , 3 ) 
     call  set_gauss_points( n_gauss_point )
     
     
     
     do r = 1 ,  n_gauss_point 
              
        if( direction_value(1) .eq. 1 ) then                        !! find shape function on face
               
            call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r)  )      ! xi = fixed
            call the_element%get_derivative_shape_function( elem_id , fixed_value , gauss_points(r) )
            d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) ) 
            d_volume = d_volume  * weight_of_gauss_point(r)  
        else
               
            call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value  )      ! eta  = fixed
            call the_element%get_derivative_shape_function( elem_id , gauss_points(r) , fixed_value )
            d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value ) 
            d_volume =  d_volume * weight_of_gauss_point(r)   
        end if 
            
       
        do  k = 1 , n_bc_nodes   
            do n = k , n_bc_nodes 
              
               boundary_shape_matrix(k,n) = boundary_shape_matrix(k,n) +   &
            &  the_element%shape_vector(indexed_set(k)) * the_element%shape_vector(indexed_set(n)) * d_volume 
              
              
               boundary_der_shape_matrix(k , n) = boundary_der_shape_matrix(k,n) + &
            &  the_element%derivative_shape_matrix(2,indexed_set(k)) *  &
            &  the_element%derivative_shape_matrix(2,indexed_set(n)) * d_volume
          
                          
            end do ! over n
        end do ! over k
                  
     end do ! over r
          
     
     !!---------------------------------------------------------------!!
     !! 2.3. symmetrize local matrices
     do r  = 1 , n_bc_nodes 
        do s = 1 , r-1
           boundary_shape_matrix(r,s)     = boundary_shape_matrix(s,r)
           boundary_der_shape_matrix(r,s) = boundary_der_shape_matrix(s,r)
        end do
     end do ! over r
     
     
!!============================================================================================================   
!!  
!! Step 3. Equation involving total pressure : eqn (3.5.a)
!!   G \ddot{P} + H P + ( a_1 / rho c ) N_{loc} * [ \dot(P) - \dot(P^i) ] + ( 1/rho c ) N_{loc} * \phi_1 = 0
!! The first two terms have already been inculded during reservior-bulk calculation.
!! The next two terms terms are now included in damping matrix.
!! The last term is added to stiffness matrix.
!! The phi_1 term appears only when N>1, for otherwise phi_1 = 0.
!! a_1 exists if N > 0 
!!
!!============================================================================================================
    
  
     col_variable_set  = 0
     row_variable_set  = 0
     
     !! 3.1. form row_var_set : total pressure variable of boundary nodes
     
     do  k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index                    
         n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
         if( n .eq. 0 ) cycle                                           ! node is absent or no dof
         row_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )  ! pressure dof of boundary node
     end do 
     
     !! 3.2. form col_var_set: total pressure and induced pressure variables of boundary nodes
     !! 3.2.1 put pressure var index in col_var_set
     
     col_variable_set(1:n_bc_nodes)  = row_variable_set(1:n_bc_nodes)
     
     num_var = n_bc_nodes
     
     !! 3.2.2 put P^i var index in col_var_set: last item in dof matrix of the node
     do  k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
         num_var = num_var + 1
         n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
         if( n .eq. 0 ) cycle                                           ! node is absent or no dof
         col_variable_set(num_var) = degree_of_freedom(n , 2 * num_dim + 3 )  ! P^i var of boundary nodes  
     end do 
     
     !! 2.3. prepare damping and stiffness matrix matrix
     elemental_damping_matrix    = 0.0_r_kind
 
     
     !! 2.3.1. set  coefficients
     alpha_coeffs = 0.0_r_kind
     if( the_model_description%order_of_reservior_higher_order_bc(1) > 0 ) then     ! if N > 0  then a_1 exists
         alpha_coeffs(1) = the_model_description%reservior_incident_coefficient(1) / ( rho * c ) ! a_1/(rho *c)
     end if
     alpha_coeffs(2) = 1.0_r_kind / ( rho * c )                                            ! 1 / (rho * c )
     
     r = n_bc_nodes
     
     !! 2.3.1.1. put local shape matrix in elemental damping matrix
     do  k = 1 , n_bc_nodes                                             ! size of row_var_set: each row a node: 3 , 4 , 7
         elemental_damping_matrix( k ,1:n_bc_nodes     ) =  alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes )   
         elemental_damping_matrix( k ,r+1:r+n_bc_nodes ) = -alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes )   
     end do 
     
     
     !! 2.3.1.2 assemble this matrix to global damping matrix
     
     if( the_model_description%order_of_reservior_higher_order_bc(1) > 0 ) then  
         call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , num_var  ) 
     end if
     
     
     !! 2.4. put term phi_1 to stiffness matrix: this needs N + M > 1
     
     !! 2.4.1. put phi_1 var index in col_var_set: for simplicity put N+M in phi_id
     
     phi_id = the_model_description%order_of_reservior_higher_order_bc(1)  + &
            & the_model_description%order_of_reservior_higher_order_bc(2)
       
     if( phi_id > 1 ) then                                              !( N + M > 1 )
         col_variable_set = 0
         
         do  k = 1 , n_bc_nodes                                         ! loop over boundary nodes with local index
             n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
             if( n .eq. 0 ) cycle                                           ! node is absent or no dof
             col_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )  ! put pressure instead of phi_1
             if( col_variable_set(k) /= 0 ) then                            ! node has pressure dof
                 col_variable_set(k) = col_variable_set(k) + 1              ! change P to phi_1 var
             end if
         end do 
         
         !! 2.4.2. put 1 / (rho * c ) * N_{loc} phi_1 in stiffness matrix
         
         elemental_stiffness_matrix = 0.0_r_kind
         
         do  k = 1 , n_bc_nodes                                        
             elemental_stiffness_matrix( k ,1:n_bc_nodes) = alpha_coeffs(2) * boundary_shape_matrix(k, 1:n_bc_nodes )     
         end do 
         
         !! 2.4.3. put  stiffness matrix in global matrices
         
         call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , n_bc_nodes  ) 
         
     end if !( if phi_id = N + M > 0 )
     
     
     
     
   
     
!!--------------------------------------------------------------------!!
!!
!!  Step 4. check if res_inf_found_node is in the current element.
!!  
!!
!!--------------------------------------------------------------------!!   
    
     special_node_index = -2
     
     do k = 1 , n_bc_nodes
          n = element_matrix( elem_id , indexed_set( k ) )              ! get global index of boundary nodes
          if( n .eq. res_inf_found_node_index ) then
              special_node_index = k ;
              exit
          end if
     end do  ! over k
     


!!---------------------------------------------------------------------!!
!!
!!  Step 5. Add equation governing induced pressure P^i
!!  
!!          The induced pressure exists even if N =0 and M = 0 
!!---------------------------------------------------------------------!!  
       
    !! The induced pressure is latest variable in var-set of boundary nodes
    !! 5.1. form row_var_set
    
    row_variable_set = 0
    col_variable_set = 0 
    
    do  k = 1 , n_bc_nodes                                              ! loop over boundary nodes with local index
                            
        n = element_matrix( elem_id , indexed_set( k ) )                ! get global index of boundary nodes
        if( n .eq. 0 ) cycle                                            ! node is absent or no dof
        if( degree_of_freedom(n,1) .eq. 0 ) cycle                       ! node has no dof
        row_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 3 )   ! induced pressure dof of boundary node
    end do 
    
    col_variable_set( 1:n_bc_nodes ) = row_variable_set( 1:n_bc_nodes )
    
    elemental_mass_matrix       = 0.0_r_kind
    elemental_stiffness_matrix  = 0.0_r_kind
    elemental_damping_matrix    = 0.0_r_kind
    elemental_force_vector      = 0.0_r_kind
    
    !! 5.2. put N_{loc} \ddot{ \PI }^i in mass matrix
    
    do  k = 1 , n_bc_nodes
        elemental_mass_matrix(k , 1 : n_bc_nodes) = boundary_shape_matrix(k, 1:n_bc_nodes )
    end do
    
    !! 5.3. put c^2 D_{loc} PI^i in stiffness matrix
    do  k = 1 , n_bc_nodes
        elemental_stiffness_matrix(k , 1 : n_bc_nodes) = (c ** 2 ) * boundary_der_shape_matrix(k, 1:n_bc_nodes )
    end do
    

    
    !! 5.4. put qc^2 \dot{ PI }^i in elemental damping matrix also force term
    
    if( special_node_index > 0 ) then
        elemental_damping_matrix = 0.0_r_kind
        elemental_force_vector   = 0.0_r_kind
        
        elemental_damping_matrix(1 , 1 ) = q * (c ** 2 )          
        elemental_force_vector(  1 , 2 ) = rho * ( c ** 2 )             !! index 2 means it is in y-direction  
    end if
    

    
    !! 5.5. send data to global matrices
    
    !! 5.5.1. send elemental mass matrix to global mass matrix
       
     call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                 & n_bc_nodes , col_variable_set , n_bc_nodes  ) 
     
     !! 5.5.2. send elemental stiffness matrix to global stiffness matrix
       
     call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                 & n_bc_nodes , col_variable_set , n_bc_nodes  ) 
      
     !! 5.5.3. send elemental damping matrix to global damping matrix : if special node is in element
     !! 5.5.4. send elemental force vector to global force vector     : if special node is in element
             
     if( special_node_index > 0 ) then
         n = element_matrix( elem_id , indexed_set( special_node_index ) ) 
         row_variable_set(1) = row_variable_set(special_node_index) 
         
         if( row_variable_set(1) /= 0  ) then                                  ! node has dof
             call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                 & 1 , row_variable_set , 1  ) 
        
             call  add_to_resultant_force( elemental_force_vector , row_variable_set , 1 )   ! only first index
         end if !( n /=0 )

     end if  ! ( special_node_index > 0 )
      
     
        

!!----------------------------------------------------------------------!!
!!
!!  Step 5. Add equation governing auxilary functions \{ \phi_{j-1} , phi_j , phi_{j+1} \} for 1 <= j <= N-1
!!          
!!  
!! if N >  1 then the equation is 
!!        (1- a_k^2) S \ddot{ theta }_{k-1} - ( a_{k+1} + a_k ) * S \dot{ theta }_k
!!      - S theta_{k+1} + c^2 D theta_{k-1} + qc^2( theta_{k-1} )_{y=0} = 0 , 
!!  
!! if N = 1 then it becomes
!!
!!
!!---------------------------------------------------------------------!! 
   
   
   do  phi_id = 1 , the_model_description%order_of_reservior_higher_order_bc(1) - 1   ! 1 <=  j <= N-1
   
       elemental_mass_matrix      = 0.0_r_kind
       elemental_damping_matrix   = 0.0_r_kind
       elemental_stiffness_matrix = 0.0_r_kind
       row_variable_set = 0
       col_variable_set = 0
       alpha_coeffs = 0.0_r_kind
       
       !! 5.1. put the terms (1- a_k^2) S \ddot{ theta }_{k-1} + c^2 D theta_{k-1} + qc^2 \dot( theta_{k-1} )_{y=0}
       !! 5.1.1. set coefficieint
       
       alpha_coeffs(1) = 1.0_r_kind - the_model_description%reservior_incident_coefficient(phi_id) ** 2  !! 1 - a_k^2
       
       !! 5.1.2. set row var set: phi_{ k }
       
       do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
          n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
          if( n .eq. 0 ) cycle                                           ! node is absent or no dof
          row_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )  ! P var of boundary nodes  
          if( row_variable_set(k) /= 0 ) then
              row_variable_set(k) = row_variable_set(k) + phi_id         ! change P to phi_k
          end if
       end do 
       
       !! 5.1.2. set col var set: if k = 1 then P and P^i appears. Otherwise phi_{k-1}
       !! 5.1.2.1. put pressure in col
       do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
          n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
          if( n .eq. 0 ) cycle                                           ! node is absent or no dof
          col_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )  ! P var of boundary nodes  
       end do
       
       !! 5.1.2.2. if k = 1 add P^i to col_var_set otherwise change P var to phi_{k-1}
       num_var = n_bc_nodes
       if( phi_id .eq. 1 ) then
           do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
              num_var = num_var + 1
              n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
              if( n .eq. 0 ) cycle                                           ! node is absent or no dof
              col_variable_set(num_var) = degree_of_freedom(n , 2 * num_dim + 3 )  ! P^i var of boundary nodes  
           end do
       else 
           do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
              if( col_variable_set(k) /= 0 ) then
                  col_variable_set(k) = col_variable_set(k) + phi_id - 1   ! change pressure var to phi_{k-1}
              end if  
           end do
       end if
       
       !! 5.1.3. form mass matrix and stiffness matrix for vector theta_{k-1}
       !! 5.1.3.1. put one copy of mass matrix 
       do  k = 1 , n_bc_nodes                                      
           elemental_mass_matrix( k ,1:n_bc_nodes ) =  alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes ) 
           elemental_stiffness_matrix( k ,1:n_bc_nodes ) =  (c ** 2) * boundary_der_shape_matrix(k, 1:n_bc_nodes )     
       end do 
       
       !! 5.1.3.2. if k = 1 then add share of P^i pressure too.
       if( phi_id .eq. 1 ) then
           r = n_bc_nodes
           do  k = 1 , n_bc_nodes                                      
               elemental_mass_matrix( k ,r+1:r+n_bc_nodes ) = -alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes )
               elemental_stiffness_matrix( k ,r+1:r+n_bc_nodes) = -(c ** 2) * boundary_der_shape_matrix(k, 1:n_bc_nodes )     
           end do 
       end if
       
       !! 5.1.4. add to global mass and stiffness matrices
       
       call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , num_var  ) 
       
       call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , num_var  )
                                     
                                     
       !! 5.2. put damping term
       !! 5.2.1. set coefficieint: if k = N then alpha_coeffs(1) = -a_k otherwise alpha_coeffs(1) = -( a_k + a_{k+1} )
       alpha_coeffs(1) = 0.0_r_kind
       alpha_coeffs(1) = - the_model_description%reservior_incident_coefficient( phi_id )    &  !! - ( a_k + a_{k+1} )
                    &    - the_model_description%reservior_incident_coefficient( phi_id + 1 )   
       
       
       !! 5.2.2. set col_var_set: theta_k
       col_variable_set = 0
       col_variable_set(1:n_bc_nodes ) = row_variable_set(1:n_bc_nodes )
       
       
       !! 5.2.3. form damping matrix
       do  k = 1 , n_bc_nodes                                      
           elemental_damping_matrix( k ,1:n_bc_nodes ) = alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes )     
       end do 
       
       !! 5.2.4. add to global matrix
       
       call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , n_bc_nodes  ) 
       
       !! 5.3. form stiffness matrix. If k < N-1 then the term phi_{k+1} appears unconditionally.
       !!      form stiffness matrix. If k = N-1 then the term phi_{k+1} appears when M /= 0.
       
        
       
       elemental_stiffness_matrix = 0.0_r_kind
       col_variable_set = 0
       num_var = 0
         
       !! put phi_{k+1} in col_var_set
       do k = 1 , n_bc_nodes                                            ! loop over boundary nodes with local index
          if( row_variable_set(k) /= 0 ) then                           ! so node has pressure dof
              col_variable_set(k) = row_variable_set(k) + 1             ! change pressure P to phi_{k+1}
          end if 
       end do
           
       !! calculate stiffness matrix - S phi_{k+1}
       do  k = 1 , n_bc_nodes                                      
           elemental_stiffness_matrix( k ,1:n_bc_nodes) = - boundary_shape_matrix(k, 1:n_bc_nodes )     
       end do 
           
       if( phi_id < the_model_description%order_of_reservior_higher_order_bc(1) - 1 ) then  !! enters unconditionally
           
           call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , n_bc_nodes  ) 
       else  !! phi_id = N-1
               
       !! if k = N - 1 : add term -S theta_{k+1} if M /= 0 as phi_{k+1} = phi_{N} is zero if M = 0.
          if( the_model_description%order_of_reservior_higher_order_bc(2) > 0 ) then
               
              call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                           & n_bc_nodes , col_variable_set , n_bc_nodes  )                    
          end if !! ( M > 0 )
           
       end if !( phi_id < N - 1 )                              
 
       
       if( special_node_index < 1 ) cycle 
       
       !! the term qc^2 (theta_{k-1} )_{y=0} : it exists for all values of k as terms like phi_0 and phi_{N-1} exist.
       
       row_variable_set(1) = row_variable_set( special_node_index )
       col_variable_set = 0
       elemental_damping_matrix = 0.0_r_kind
       n = element_matrix( elem_id , indexed_set( special_node_index ) )      ! get global index of boundary nodes

       
       if( phi_id .eq. 1 ) then
           col_variable_set(1) = degree_of_freedom(n , 2 * num_dim + 2 ) ! put P in col var set
           col_variable_set(2) = degree_of_freedom(n , 2 * num_dim + 3 ) ! put P^i in col var set
           elemental_damping_matrix(1,1) = q * ( c**2)
           elemental_damping_matrix(1,2) = - elemental_damping_matrix(1,1)
           
           call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & 1 , col_variable_set , 2  )
       else 
           col_variable_set(1) = degree_of_freedom(n , 2 * num_dim + 2 ) ! put P instead of phi_{k-1}
           if( col_variable_set(1) /= 0 ) then
               col_variable_set(1) = col_variable_set(1) + phi_id - 1
               elemental_damping_matrix(1,1) = q * ( c**2)
               call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & 1 , col_variable_set , 1  )
           end if
       end if !!( phi_id = 1 )
                                  
       
   end do ! over phi_id
  
     
!!====================================================================!!
!!
!!  6. case phi_id = N : it combines incident and evanscent waves.
!!                    If M =0 then the remaining equations will not appears. 
!!  In the special case of O(1,1) the following equation uses phi_0 = P - P^i                  
!!  This case is not covered in the previous loop.
!!====================================================================!!
  
  !! if M = 0 none of subsequent equations contribute.
  
   if( the_model_description%order_of_reservior_higher_order_bc(2) .eq. 0 ) return
   

   phi_id = the_model_description%order_of_reservior_higher_order_bc(1)     !! put N in phi_id
       
   elemental_mass_matrix      = 0.0_r_kind
   elemental_damping_matrix   = 0.0_r_kind
   elemental_stiffness_matrix = 0.0_r_kind
   row_variable_set = 0
   col_variable_set = 0
   alpha_coeffs = 0.0_r_kind
       
   !!  put the terms (1- a_k^2) S \ddot{ theta }_{k-1} + c^2 D theta_{k-1} + qc^2 \dot( theta_{k-1} )_{y=0}
   !!  set coefficieint
       
   alpha_coeffs(1) = 1.0_r_kind - the_model_description%reservior_incident_coefficient(phi_id) ** 2  !! 1 - a_k^2
       
   !!  set row var set: phi_{ k } = phi_N
       
   do k = 1 , n_bc_nodes                                                ! loop over boundary nodes with local index
      n = element_matrix( elem_id , indexed_set( k ) )                  ! get global index of boundary nodes
      if( n .eq. 0 ) cycle                                              ! node is absent or no dof
      row_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )     ! P var of boundary nodes  
      if( row_variable_set(k) /= 0 ) then
          row_variable_set(k) = row_variable_set(k) + phi_id            ! change P to phi_k
      end if
    end do 
       
    !!  set col var set: if k = 1 then P and P^i appears. Otherwise phi_{k-1}
    !!  put pressure in col
    do k = 1 , n_bc_nodes                                               ! loop over boundary nodes with local index
       n = element_matrix( elem_id , indexed_set( k ) )                 ! get global index of boundary nodes
       if( n .eq. 0 ) cycle                                             ! node is absent or no dof
       col_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )    ! P var of boundary nodes  
    end do
       
   !! if phi_id = N = 1 add P^i to col_var_set otherwise change P var to phi_{k-1}
   num_var = n_bc_nodes
   if( phi_id .eq. 1 ) then
       do k = 1 , n_bc_nodes                                            ! loop over boundary nodes with local index
          num_var = num_var + 1
          n = element_matrix( elem_id , indexed_set( k ) )              ! get global index of boundary nodes
          if( n .eq. 0 ) cycle                                          ! node is absent or no dof
          col_variable_set(num_var) = degree_of_freedom(n , 2 * num_dim + 3 )  ! P^i var of boundary nodes  
       end do
   else 
       do k = 1 , n_bc_nodes                                            ! loop over boundary nodes with local index
          if( col_variable_set(k) /= 0 ) then
              col_variable_set(k) = col_variable_set(k) + phi_id - 1    ! change pressure var to phi_{k-1}
          end if  
       end do
   end if
       
  !!  form mass matrix and stiffness matrix for vector theta_{k-1}
  !!  put one copy of mass matrix 
  
  do  k = 1 , n_bc_nodes                                      
      elemental_mass_matrix( k ,1:n_bc_nodes ) =  alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes ) 
      elemental_stiffness_matrix( k ,1:n_bc_nodes ) =  (c ** 2) * boundary_der_shape_matrix(k, 1:n_bc_nodes )     
  end do 
       
  !! if phi_id = N = 1 then add share of P^i pressure too.
  if( phi_id .eq. 1 ) then
      r = n_bc_nodes
      do  k = 1 , n_bc_nodes                                      
          elemental_mass_matrix( k ,r+1:r+n_bc_nodes ) = -alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes )
          elemental_stiffness_matrix( k ,r+1:r+n_bc_nodes) = -(c ** 2) * boundary_der_shape_matrix(k, 1:n_bc_nodes )     
      end do 
  end if
       
  !! add to global mass and stiffness matrices
       
  call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                              & n_bc_nodes , col_variable_set , num_var  ) 
       
  call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                              & n_bc_nodes , col_variable_set , num_var  )
                                     
                                     
  !!  put damping term -a_k S dot{ phi}_N 
  !!  set coefficieint: alpha_coeffs(1) = -a_N 
  alpha_coeffs(1) = 0.0_r_kind
  alpha_coeffs(1) = -the_model_description%reservior_incident_coefficient( phi_id ) !! -a_N
  
       
  !!  set col_var_set: theta_k
  col_variable_set = 0
  col_variable_set(1:n_bc_nodes ) = row_variable_set(1:n_bc_nodes )
       
  elemental_stiffness_matrix = 0.0_r_kind  
    
  !! 5.2.3. form damping and stiffness matrices for these terms.
  do  k = 1 , n_bc_nodes                                      
      elemental_damping_matrix( k ,1:n_bc_nodes )   = alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes )       
  end do 
       
  !! add to global matrix
       
  call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , n_bc_nodes  ) 
                                     
       
  !! 5.3. form stiffness matrix. 
  !!      If k = N then two terms appear: phi_{N} with coeff b_1 enters when M> 0
  !!                                      and phi_{N+1} contributes if M > 1  as phi_{N+M} = 0
        
       
   elemental_stiffness_matrix = 0.0_r_kind
   col_variable_set = 0
   num_var = 0
       
  
   !!  If k = N then two terms appear: phi_{N} with coeff -b_1 enters when M> 0 ( this condition is satisfied here)
   !!                                      and phi_{N+1} contributes if M > 1  as phi_{N+M} = 0
           
   num_var = 0
   col_variable_set = 0
           
   
   !! set coefficieint -b_1
   alpha_coeffs(1) = -the_model_description%reservior_evanescent_coefficient(1)
               
   !! put -b_1 phi_{N} in stiffness matrix
   col_variable_set(1:n_bc_nodes ) = row_variable_set(1:n_bc_nodes)
               
   num_var = n_bc_nodes
               
   !! term -S theta_{N+1} appears if  M > 1
   if( the_model_description%order_of_reservior_higher_order_bc(2) > 1 ) then
       do k = 1 , n_bc_nodes                                            ! loop over boundary nodes with local index
          num_var = num_var + 1
          n = element_matrix( elem_id , indexed_set( k ) )              ! get global index of boundary nodes
          if( n .eq. 0 ) cycle                                          ! node is absent or no dof
          col_variable_set(num_var) = degree_of_freedom(n , 2 * num_dim + 2 ) ! put P instead of phi_{k+1}
          if( col_variable_set(num_var) /= 0 ) then                           ! so node has pressure dof
              col_variable_set(num_var) = col_variable_set(num_var) + phi_id + 1  ! change pressure P to phi_{N+1}
          end if 
       end do
   end if !! ( M > 1 )
               
               
   !! form stiffness matrix for term -b_1 S phi_N
   do  k = 1 , n_bc_nodes                                      
       elemental_stiffness_matrix( k ,1:n_bc_nodes) = alpha_coeffs(1)  * &
                                   &  boundary_shape_matrix(k, 1:n_bc_nodes )     
   end do 
               
   !! form stiffness matrix for term - S phi_{N+1}
   r = n_bc_nodes
   if( the_model_description%order_of_reservior_higher_order_bc(2) > 1 ) then
       do  k = 1 , n_bc_nodes                                      
           elemental_stiffness_matrix( k ,r + 1 :r + n_bc_nodes) =  &
                                       &  - boundary_shape_matrix(k, 1:n_bc_nodes )     
       end do 
   end if !! ( M > 1 )
               
   call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                               & n_bc_nodes , col_variable_set , num_var  )
           
  
       
   if( special_node_index > 0  ) then
       
       !! the term qc^2 (theta_{k-1} )_{y=0} : it exists for all values of k as terms like phi_0 and phi_{N-1} exist.
       
       row_variable_set(1) = row_variable_set( special_node_index )
       col_variable_set = 0
       elemental_damping_matrix = 0.0_r_kind
       n = element_matrix( elem_id , indexed_set( special_node_index ) )      ! get global index of boundary nodes

       
       if( phi_id .eq. 1 ) then
           col_variable_set(1) = degree_of_freedom(n , 2 * num_dim + 2 ) ! put P in col var set
           col_variable_set(2) = degree_of_freedom(n , 2 * num_dim + 3 ) ! put P^i in col var set
           elemental_damping_matrix(1,1) = q * ( c**2)
           elemental_damping_matrix(1,2) = - elemental_damping_matrix(1,1)
           
           call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & 1 , col_variable_set , 2  )
       else 
           col_variable_set(1) = degree_of_freedom(n , 2 * num_dim + 2 ) ! put P instead of phi_{k-1}
           if( col_variable_set(1) /= 0 ) then
               col_variable_set(1) = col_variable_set(1) + phi_id - 1
               elemental_damping_matrix(1,1) = q * ( c**2)
               call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & 1 , col_variable_set , 1  )
           end if
       end if !!( phi_id = 1 )
                                  
     
   end if !! special node index > 0
  
     
     
          
!!----------------------------------------------------------------------!!
!!
!!  Step 6. Add equation governing auxilary functions \{ \phi_{N+ j-1} , phi_{N+j} , phi_{N + j+1} \} 
!!          for 1 <= j <= M-1.
!!
!!          If j = M-1 then the term phi_{N + j+1} = phi_{N+M} = 0 and will not enters.
!!
!!---------------------------------------------------------------------!!   
       
   do  phi_id = 1 , the_model_description%order_of_reservior_higher_order_bc(2) - 1
       
       elemental_mass_matrix      = 0.0_r_kind
       elemental_damping_matrix   = 0.0_r_kind
       elemental_stiffness_matrix = 0.0_r_kind
       row_variable_set = 0
       col_variable_set = 0
       alpha_coeffs     = 0.0_r_kind
       
       !! 6.1. form row_var_set: term phi_{N+k}
       do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
          n = element_matrix( elem_id , indexed_set( k ) )               ! get global index of boundary nodes
          if( n .eq. 0 ) cycle                                           ! node is absent or no dof
          row_variable_set(k) = degree_of_freedom(n , 2 * num_dim + 2 )  ! P var of boundary nodes  
          if( row_variable_set(k) /= 0 ) then
              row_variable_set(k) = row_variable_set(k) + phi_id  &       ! change P to phi_{N+k}
                               &  + the_model_description%order_of_reservior_higher_order_bc(1)  
          end if
       end do 
       
       !! 6.2. put term S \ddot{ theta }_{N+k-1} in mass matrix
       !! 6.2.1. form col_var_set
       
       do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index
          if( row_variable_set(k) /= 0 ) then
              col_variable_set(k) = row_variable_set(k) - 1             !! change phi_{N+k} to phi_{N+k-1}
          end if
       end do 
       
       !! 6.2.2. form mass matrix
       do  k = 1 , n_bc_nodes                                      
           elemental_mass_matrix( k , 1:n_bc_nodes) = boundary_shape_matrix(k, 1:n_bc_nodes )     
       end do 
       
       !! 6.2.3. add to global matrices
       call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , n_bc_nodes  )
       
       
       !! 6.3. add following terms to stiffness matrix
       !!      -b_k^2 S theta_{N+k-1} + c^2 D theta_{N+k-1} - ( b_k + b_{k+1} ) S theta_{N+k}
       !!      - S theta_{N+k+1}
       !! The last term appears when k < M-1
       
       !! 6.3.1. form col_var_set
       col_variable_set = 0
       
       !! 6.3.1.1. put term theta_{N+k-1} in col_var_set
       do k = 1 , n_bc_nodes                                             ! loop over boundary nodes with local index 
          if( row_variable_set(k) /= 0 ) then
              col_variable_set(k) = row_variable_set(k) -1 
          end if
       end do 
       
       num_var = n_bc_nodes
       
       !! 6.3.1.2. put term theta_{N+k} in col_var_set
       do k = 1 , n_bc_nodes 
          num_var = num_var + 1
          if( row_variable_set(k) /= 0 ) then
              col_variable_set(num_var) = row_variable_set(k)           !! phi_{N+k}
          end if
       end do 
       
       
       !! 6.3.1.2. put term theta_{N+k+1} in col_var_set if k < M-1
       if( phi_id < the_model_description%order_of_reservior_higher_order_bc(2) - 1 ) then
           do k = 1 , n_bc_nodes 
              num_var = num_var + 1
              if( row_variable_set(k) /= 0 ) then
                  col_variable_set(num_var) = row_variable_set(k) + 1   !! phi_{N+k+1}
              end if
           end do
       end if
       
       
       !! 6.4. form stiffness matrix
       !! 6.4.1. form coefficieints
       alpha_coeffs = 0.0_r_kind
       alpha_coeffs(1) = - the_model_description%reservior_evanescent_coefficient(phi_id) ** 2  !! -b_k^2
       alpha_coeffs(2) = - the_model_description%reservior_evanescent_coefficient(phi_id)  &   !! -b_k 
                       & - the_model_description%reservior_evanescent_coefficient(phi_id+1)
       
       
       !! 6.4.1. put terms -b_k^2 S theta_{N+k-1} + c^2 D theta_{N+k-1} in stiffness matrix
       do  k = 1 , n_bc_nodes                                      
           elemental_stiffness_matrix( k , 1:n_bc_nodes) = alpha_coeffs(1) * boundary_shape_matrix(k, 1:n_bc_nodes ) 
           elemental_stiffness_matrix( k , 1:n_bc_nodes) = elemental_stiffness_matrix( k , 1:n_bc_nodes)  &
                     &  + ( c ** 2 ) *   boundary_der_shape_matrix(k, 1:n_bc_nodes )  
       end do 
       
       r = n_bc_nodes
       
       !! 6.4.1. put terms - ( b_k + b_{k+1} ) S theta_{N+k} in stiffness matrix
       do  k = 1 , n_bc_nodes                                      
           elemental_stiffness_matrix( k , r+1: r + n_bc_nodes) = alpha_coeffs(2) *   &
                                       &   boundary_shape_matrix(k, 1:n_bc_nodes ) 
       end do 
       
       
       !! 6.4.2. put terms - S theta_{N+k+1} in stiffness matrix
       if( phi_id < the_model_description%order_of_reservior_higher_order_bc(2) - 1 ) then
           r = r + n_bc_nodes
           do  k = 1 , n_bc_nodes                                      
               elemental_stiffness_matrix( k , r+1: r + n_bc_nodes) = - boundary_shape_matrix(k, 1:n_bc_nodes ) 
           end do 
       end if !! ( phi_id < M-1)
       
       
       !! 6.5. add to global matrix
       call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                     & n_bc_nodes , col_variable_set , num_var  )
       
       
       if( special_node_index < 1 ) cycle
       
       row_variable_set(1) = row_variable_set( special_node_index )
       col_variable_set = 0
       elemental_damping_matrix = 0.0_r_kind
       n = element_matrix( elem_id , indexed_set( special_node_index ) )      ! get global index of boundary nodes
       col_variable_set(1) = degree_of_freedom(n , 2 * num_dim + 2 )          ! P var of boundary nodes 
       
       if( col_variable_set(1) /= 0 ) then
           col_variable_set(1) = col_variable_set(1) + phi_id - 1   &
                            &  + the_model_description%order_of_reservior_higher_order_bc(1)
                            
           elemental_damping_matrix(1,1) = q * ( c ** 2 )
           
           call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                     & 1 , col_variable_set , 1  )
       end if 
           
   end do !! over phi_id: 1 : M-1 
   
   
     

  end subroutine GN_model_reservior_truncation_boundary
   
   

  
 
 
 
 
!!====================================================================!!
!!
!!  The fem analysis in reservior_truncation_boundary in HW model 
!!
!! 
!!====================================================================!!
   
 module subroutine HW_model_reservior_truncation_boundary( elem_id , bc_id )  
   
    use mod_physics  , only : element_property , material , boundary_load , alpha_reflection_coefficient
    use mod_geometry , only : num_dim , degree_of_freedom , element_matrix
    use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face
    use mod_physics  , only : the_model_description , add_to_resultant_force , add_to_global_matrices
    use mod_fem      , only : set_gauss_points
    implicit none
 
    integer , intent( in ) :: elem_id , bc_id
    integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point
    integer :: i , k , n , r, s , phi_id , special_node_index, var_id, row , col, num_row , num_col
    real( kind = r_kind ) :: rho , c, fixed_value , d_volume , q 
    real( kind = r_kind ) :: local_mass_matrix(3,3) , local_stiffness_matrix(3,3) , local_damping_matrix(3,3) 
   
    
    !!----------------------------------------------------------------!
    !! step 1. get physical property of the element
     c   = material( element_property( elem_id , 2) )%elastic_module     ! used as c (speed of sound)
     rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
     
     q = ( 1.0_r_kind - alpha_reflection_coefficient ) / ( 1.0_r_kind + alpha_reflection_coefficient ) /c 

  
!!---------------------------------------------------------------------!!
!!
!!  Step 2. Form boundary_shape_matrix and boundary_der_shape_matrix.
!!
!!---------------------------------------------------------------------!!    
     
     !!---------------------------------------------------------------!!
     !! 2.1. get local node index of boundary nodes
     
     indexed_set = 0
     
     n_face = boundary_load( bc_id )%face_id
     n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
   
     call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
     fixed_value = real( direction_value(2) , kind = r_kind )           !! either xi or eta or zeta is fixed at +1 or -1
     
     
     
      
     !!---------------------------------------------------------------!!
     !! 2.2.  calcuate local matrices for the three facial nodes: 3 , 4 , 7 in local index for Q8 element
     
     local_mass_matrix      = 0.0_r_kind
     local_stiffness_matrix = 0.0_r_kind 
     
     
     
     n_gauss_point = element_property( elem_id , 3 ) 
     call  set_gauss_points( n_gauss_point )
     
     
     do r = 1 ,  n_gauss_point 
              
        if( direction_value(1) .eq. 1 ) then                        !! find shape function on face
               
            call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r)  )      ! xi = fixed
            call the_element%get_derivative_shape_function( elem_id , fixed_value , gauss_points(r) )
            d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) ) 
            d_volume = d_volume  * weight_of_gauss_point(r)  
        else
               
            call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value  )      ! eta  = fixed
            call the_element%get_derivative_shape_function( elem_id , gauss_points(r) , fixed_value )
            d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value ) 
            d_volume =  d_volume * weight_of_gauss_point(r)   
        end if 
            
      
        do  k = 1 , n_bc_nodes   
            do n = k , n_bc_nodes 
              
               local_mass_matrix(k,n) = local_mass_matrix(k,n) +   &
            &  the_element%shape_vector(indexed_set(k)) * the_element%shape_vector(indexed_set(n)) * d_volume 
              
              
               local_stiffness_matrix(k , n) = local_stiffness_matrix(k,n) + &
            &  the_element%derivative_shape_matrix(2,indexed_set(k)) *  &
            &  the_element%derivative_shape_matrix(2,indexed_set(n)) * d_volume
          
                          
            end do ! over n
        end do ! over k
                  
     end do ! over r
          
    
     
     !!---------------------------------------------------------------!!
     !! 2.3. symmetrize local matrices
     do r  = 1 , n_bc_nodes
        do s = 1 , r-1
           local_mass_matrix(r,s)     = local_mass_matrix(s,r)
           local_stiffness_matrix(r,s) = local_stiffness_matrix(s,r)
        end do
     end do ! over r
     
     
     !! 2.4. divide by rho
     local_mass_matrix       = ( 1.0_r_kind / rho ) * local_mass_matrix ;
     local_stiffness_matrix  = ( 1.0_r_kind / rho ) * local_stiffness_matrix ;
     local_damping_matrix   = 0.0_r_kind
     local_damping_matrix(2,2) = q / rho ;     ! for special node 
                               
     
!!---------------------------------------------------------------------!!
!!
!!  Step 3. Form row and col variable set index
!!
!!---------------------------------------------------------------------!!
    
     col_variable_set  = 0
     row_variable_set  = 0
     
     !! 3.1. get number of variables: N+M+3: get max num of variables for three boundary nodes
     num_var = 0 ;
     
     n = element_matrix( elem_id , indexed_set( 1 ) ) ;   !! first node
     
     
     do  i = 1 , n_bc_nodes
         n = element_matrix( elem_id , indexed_set( i ) ) ;   !! each boundary node
         k = degree_of_freedom (n, 2 * num_dim + 3)  - degree_of_freedom (n, 2 * num_dim + 2) + 1 ;
         num_var  = max( num_var , k ) ;
     end do
     
     ! ! 3.2. form index set for bigmatrices
     
     var_id = 0 ;
     do  i  = 1 , num_var
         do  k = 1 , n_bc_nodes
             var_id  = var_id + 1 ;
             n = element_matrix( elem_id , indexed_set( k ) ) ;   !! each boundary node
             
             if( degree_of_freedom (n, 2 * num_dim + 2) .eq. 0 ) cycle  ;  !! node at y =H
             
             row_variable_set( var_id ) = degree_of_freedom (n, 2 * num_dim + 2) + i -1 ;
         end do
     end do 
     
    num_row = var_id  ;
    num_col = num_row ;
 
     
     col_variable_set(1:num_col ) = row_variable_set(1:num_row ) ;
     
    

    !! 3.3. prepare elemental matrices
    
     elemental_damping_matrix    = 0.0_r_kind
     elemental_mass_matrix       = 0.0_r_kind
     elemental_stiffness_matrix  = 0.0_r_kind
     elemental_force_vector      = 0.0_r_kind
     
    
     
  
   
     
      
     
!!---------------------------------------------------------------------!! 
!!  
!! Step 4. Equation involving total pressure
!!   G \ddot{P} + H P + ( a_0 / c ) M_{loc} * [ \dot(P) - \dot(P^i) ] - ( 1/ c ) M_{loc} * \dot(\phi)_1 = 0
!! The first two terms have already been inculded during reservior-bulk calculation.
!! The  last terms are now included in damping matrix.
!!
!! 1. The term phi_1 comes from first incident wave, so if enters when N > 0
!! 2. This modification only appears in damping matrix
!!
!!---------------------------------------------------------------------!!
  
    
    !! step 4.1. from coefficient.
    alpha_coeffs = 0.0_r_kind
    alpha_coeffs(1) = the_model_description%reservior_incident_coefficient(1) / c  ;  ! a_0_bar =  a_0/c
    alpha_coeffs(2) = -1.0_r_kind / c  ;
    
    !! 4.2. P row ,  P col ,  term: ( a_0 / c ) M_{loc}  
    row = 0 ;
    col = 0 ;

    elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1)  * local_mass_matrix ;
    
    !! 4.3. P row , P_i col ,  term: -( a_0 / c ) M_{loc} 
    
    col  = ( num_var - 1 ) * n_bc_nodes  ;
    elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = - alpha_coeffs(1) * local_mass_matrix ;
    
    !! 4.4. P row , phi_1 col ,  term: -( 1 / c ) M_{loc}  if N > 0 
    
    if( the_model_description%order_of_reservior_higher_order_bc(1)  > 0 )  then
        col = n_bc_nodes 
        elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(2) * local_mass_matrix ;
    end if
    
   
!!--------------------------------------------------------------------!!
!!
!! 5. check if res_inf_found_node is in the current element.
!!
!!--------------------------------------------------------------------!! 

     special_node_index = -2
     
     do k = 1 , n_bc_nodes
          n = element_matrix( elem_id , indexed_set( k ) )              ! get global index of boundary nodes
          if( n .eq. res_inf_found_node_index ) then
              special_node_index = k ;
              exit
          end if
     end do  ! over k
  
     
    
   
  
!!---------------------------------------------------------------------!!
!!
!!  Step 6. Add equation governing induced pressure P^i
!!  
!!          The induced pressure exists even if N =0 and M = 0 
!!  
!!  Eqn: (1/c^2) * M * \ddot{P}_i +  K * P^i +  (q/rho) * C * \dot{P}^i  = a_y^g
!!
!!  Right hand side: number 1.0 is added to P^i_4, i.e. special node index
!!  Where is it?   row_var_set( ( num_var - 1 ) * n_bc_nodes + special_node_index ) 
!!
!!  The force is  F = [ 0  , 1.0 ] , i.e. in the y-direction
!!  At the end put:
!!  
!!                if( special_node_index > 0 ) then
!!                    col_var_set  = 0 ;
!!                    n  = row_var_set( ( num_var - 1 ) * n_bc_nodes + special_node_index ) 
!!                    col_var_set(1)  = n ;
!!                    F = 0.0 ;
!!                    F(1,2) = 1.0 ;
!!                    add_to_global_force( F , col_var_set , 1 ) ;
!!                end if
!!
!!---------------------------------------------------------------------!!  
       
    
    !! 6.1. set coeff
    alpha_coeffs = 0.0_r_kind ;
    alpha_coeffs(1) = 1.0_r_kind / ( c ** 2 ) ;
    
    
    !! 6.2. put (1/c^2) * M * \ddot{P}_i in mass matrix
    
    row = ( num_var - 1 ) * n_bc_nodes ;
    col = ( num_var - 1 ) * n_bc_nodes ;
    elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ;
    
    !! 6.3. put K * P^i in stiffness matrix
    
    elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) =  local_stiffness_matrix ;
    
    
    !! 6.4. put (q/rho) * C * \dot{P}^i  in damping matrix
    if( special_node_index > 0 )  then
        elemental_damping_matrix( row +1:row +3 , col +1:col +3 )=  local_damping_matrix
        elemental_force_vector(row + special_node_index , special_node_index )   = 1.0_r_kind ;
    end if
    
    
    
    
!!---------------------------------------------------------------------!!
!!
!!  Step 7. Add equation governing auxilary functions \{ \phi_ 0, \phi_1 , \phi_2 \}
!!          Output is phi_1: it is assigned to phi_1
!!  This equation is added if N > 0. 
!!  The term phi_2 is included if N>1
!!---------------------------------------------------------------------!!  
     
    if( the_model_description%order_of_reservior_higher_order_bc(1) > 0 ) then   !! if N > 0 
        

        !! 7.1. set coefficients used in this equation
        !!------------------------------------------------------------!!
        !! set alpha(1) = -2( a_bar_1 / c_bar ) * alpha_0 = 2 a_1 ( 1 - a_0^2 ) / c^2 , 
        !!     alpha(2) =  2 a_bar_0 * a_bar_1 + beta_1 =  ( 1 + a_1 ( a_1 + 2 a_0 ) )/ c^2 ,
        !!     alpha(3) = -alpha_1 = ( 1 - a_1^2 )/ c^2
        !!     alpha(4) =  2 a_bar_1 / c_bar =  2 * a_1 
        !!------------------------------------------------------------!!
        
        
        alpha_coeffs = 0.0_r_kind
        alpha_coeffs(1) = 1.0_r_kind - the_model_description%reservior_incident_coefficient(1) ** 2  !! ( 1 - a_0^2)
        alpha_coeffs(1) = 2.0_r_kind * alpha_coeffs(1) * the_model_description%reservior_incident_coefficient(2) !! *  2a_1
        alpha_coeffs(1) = alpha_coeffs(1) / ( c ** 2 ) ;                                                         !! / c^2
        
        
        alpha_coeffs(2) = the_model_description%reservior_incident_coefficient(2) +  &
                    &  2.0_r_kind * the_model_description%reservior_incident_coefficient(1)  ! a_1 + 2 a_0
        
        alpha_coeffs(2) = alpha_coeffs(2) * the_model_description%reservior_incident_coefficient(2) !! * a_1
        alpha_coeffs(2) = alpha_coeffs(2) + 1.0_r_kind
        alpha_coeffs(2) = alpha_coeffs(2) / ( c ** 2 ) ;                                            !! / c^2
        
        
        alpha_coeffs(3) = 1.0_r_kind - the_model_description%reservior_incident_coefficient(2) ** 2  !! ( 1 - a_1^2)
        alpha_coeffs(3) = alpha_coeffs(3) / ( c ** 2 ) ;                                             !! / c^2
        
        alpha_coeffs(4) = 2.0_r_kind * the_model_description%reservior_incident_coefficient(2)       !! 2 a_1 
        
        row = n_bc_nodes ;  !! added to phi_1 equation
    
        !! 7.1.2.  row: phi_1 , col: P
        col = 0 ;
        elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ;
        
        elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(4) * local_stiffness_matrix ;
        
        if( special_node_index > 0 ) then
            elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(4) * local_damping_matrix ;
        end if
        
        !! 7.1.3.  row: phi_1 , col: phi_1
        col = n_bc_nodes ;
        elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(2) * local_mass_matrix ;
        
        elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) =  local_stiffness_matrix ;
        
        if( special_node_index > 0 ) then
            elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) =  local_damping_matrix ;
        end if
        
        !! 7.1.4.  row: phi_1 , col: phi_2
        !! added if N+M+1 > 2 or N+M > 1
        
        if( the_model_description%order_of_reservior_higher_order_bc(1) + &     !! if ( N + M > 1 )
          & the_model_description%order_of_reservior_higher_order_bc(2) > 1 )  then
          
          col = 2 * n_bc_nodes ; 
          elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(3) * local_mass_matrix ;   
          
          elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) =  local_stiffness_matrix ;
          
          if( special_node_index > 0 ) then
              elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) =  local_damping_matrix ;
          end if
        
        end if !! N+M > 1
        
        
        !! 7.1.5.  row: phi_1 , col: P^i
        col = ( num_var - 1 ) * n_bc_nodes
        elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) =-alpha_coeffs(1) * local_mass_matrix ; 
        elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = -alpha_coeffs(4) * local_stiffness_matrix ;
        
        if( special_node_index > 0 ) then
            elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = -alpha_coeffs(4) * local_damping_matrix ;
        end if
 
        
    end if  ! ( N > 0 )  
     
    
    

     
!!---------------------------------------------------------------------!!
!!
!!  Step 8. Add equations governing auxilary functions \{ \phi_{j-1} , \phi_j , \phi_{j+1} \} 
!!          for  2 <= j <= N 
!!          This equation appears if N > 1 so writing loop over j=2:N 
!!          automatically does not include the case of N = 1.
!!
!!          The term phi_{N+1} does not appear if M = 0,
!! 
!!---------------------------------------------------------------------!!  

    do  phi_id = 2 , the_model_description%order_of_reservior_higher_order_bc(1)
        
        
        !! 8.1. set coefficients used in this equation
        !!------------------------------------------------------------!!
        !! set alpha(1) = - a_bar_j * alpha_{j-1} = a_j * ( 1 - a_{j-1}^2 ) / c^3 ;
        !!     alpha(2) = gamma_j =  ( a_{j-1} + a_j ) * ( 1 + a_j* a_{j-1} ) / c^3 ;
        !!     alpha(3) = -a_bar_{j-1} * alpha_j = a_{j-1} * ( 1 - a_{j} ^2 ) / c^3 ;
        !!     alpha(4) =  a_bar_j =  a_j / c ;
        !!     alpha(5) = xi_{j-1} = ( a_{j-1} +a_j ) / c ;
        !!     alpha(6) = a_bar_{j-1} = a_{j-1} / c ;
        !!------------------------------------------------------------!!
        
        alpha_coeffs = 0.0_r_kind 
        
        !! 8.1.1. put coeff of ddot{ theta }_{j-1} in alpha_coeffs(1): a_j * ( 1 - a_{j-1} ^2 ) / c^3
        alpha_coeffs(1) = 1.0_r_kind - the_model_description%reservior_incident_coefficient( phi_id ) ** 2
        alpha_coeffs(1) = alpha_coeffs(1) * the_model_description%reservior_incident_coefficient( phi_id + 1 )
        alpha_coeffs(1) = alpha_coeffs(1)  / ( c ** 3 ) ;
        
        !! 8.1.2. put coeff of ddot{ theta }_{j} in alpha_coeffs(2): ( a_{j-1} + a_j ) * ( 1 + a_j* a_{j-1} ) / c^3
        alpha_coeffs(2) = 1.0_r_kind + the_model_description%reservior_incident_coefficient( phi_id )  * &
                                     & the_model_description%reservior_incident_coefficient( phi_id + 1 )
                                    
        alpha_coeffs(2) = alpha_coeffs(2) * ( the_model_description%reservior_incident_coefficient( phi_id )  &
                                          +   the_model_description%reservior_incident_coefficient( phi_id + 1) ) ;
                                          
        alpha_coeffs(2) = alpha_coeffs(2) / ( c ** 3 ) ;
        
        
        !! 8.1.3. put coeff of ddot{ theta }_{j+1} in alpha_coeffs(3): a_{j-1} * ( 1 - a_{j} ^2 ) / c^3
        alpha_coeffs(3) = 1.0_r_kind - the_model_description%reservior_incident_coefficient( phi_id + 1 ) ** 2 ;
        alpha_coeffs(3) = alpha_coeffs(3) * the_model_description%reservior_incident_coefficient( phi_id ) ;
        
        alpha_coeffs(3) = alpha_coeffs(3) / ( c ** 3 ) ;
        
        !! 8.1.4. put coeff of  theta_{j-1} in alpha_coeffs(4):  a_j / c
        alpha_coeffs(4) = the_model_description%reservior_incident_coefficient( phi_id + 1 ) / c ;
        
        !! 8.1.5. put coeff of  theta_{j} in alpha_coeffs(5):  ( a_{j-1} +a_j ) / c
        alpha_coeffs(5) = the_model_description%reservior_incident_coefficient( phi_id )   &
                      & + the_model_description%reservior_incident_coefficient( phi_id + 1 )
                        
        alpha_coeffs(5) = alpha_coeffs(5) / c ;
        
        !! 8.1.6. put coeff of  theta_{j+1} in alpha_coeffs(6):  a_{j-1} / c
        alpha_coeffs(6) = the_model_description%reservior_incident_coefficient( phi_id ) / c ;
        
        !! This is added to phi_j term for j= 2 ,..., N
        
        row = phi_id * n_bc_nodes ;
        
        !! 8.2. term  phi_{j-1}
        col = ( phi_id - 1 ) * n_bc_nodes ;
        
        elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ; 
        elemental_stiffness_matrix(row +1:row +3 , col +1:col +3) = alpha_coeffs(4) * local_stiffness_matrix ;
        
        if( special_node_index > 0 ) then
            elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(4) * local_damping_matrix ;
        end if
        
        
        !! 8.3. term  phi_{j}
        col = phi_id  * n_bc_nodes ;
        
        elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(2) * local_mass_matrix ; 
        elemental_stiffness_matrix(row +1:row +3 , col +1:col +3) = alpha_coeffs(5) * local_stiffness_matrix ;
        
        if( special_node_index > 0 ) then
            elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(5) * local_damping_matrix ;
        end if
        
        
        !! 8.4. term phi_{j+1}: for the last term phi_{N+1} appears only if M > 0
        if( ( phi_id  .eq. the_model_description%order_of_reservior_higher_order_bc(1) )  .and. &
          & ( the_model_description%order_of_reservior_higher_order_bc(2) .eq. 0 ) ) cycle ;
          
        
        !! 8.4.1. if j < N this term is added to  col = ( phi_id  + 1 ) * n_bc_nodes ;
        !!        if J= N  this term is added to  col = ( phi_id  + 2 ) * n_bc_nodes  as psi is intermediate term
        
        col = ( phi_id  + 1 ) * n_bc_nodes ;
        if( phi_id  .eq. the_model_description%order_of_reservior_higher_order_bc(1) ) then
            col = col + n_bc_nodes ;                                                          !! to include psi
        end if
        
        elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(3) * local_mass_matrix ; 
        elemental_stiffness_matrix(row+1:row +3,col+1:col+3)= alpha_coeffs(6)* local_stiffness_matrix ;
                
        if( special_node_index > 0 ) then
            elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(6) * local_damping_matrix ;
        end if
        
             
                                 
    end do ! over phi_id
    
    

!!---------------------------------------------------------------------!!
!!
!!  Step 9. Add equations containing D_t psi
!!          All the remaining equations appears only if M > 0
!!  
!!---------------------------------------------------------------------!! 

    !! 9.1. check if M > 0 and N >= 1 as it used it used one of the equations phi_j where j \ne 0
    
    if ( ( the_model_description%order_of_reservior_higher_order_bc(2) > 0 ) .and. &
       & ( the_model_description%order_of_reservior_higher_order_bc(1) > 0 ) ) then
   
       !! 9.1. set coefficients used in this equation
       !!------------------------------------------------------------!!
       !! set alpha(1) = - alpha_N = ( 1 - a_N^2 ) / c^2 ;
       !!     alpha(2) =   beta_N  =  ( 1 + a_N^2 ) / c^2 ;
       !!     alpha(3) = 2 * a_bar_N = 2 * a_N / c  ;
       !!------------------------------------------------------------!!
       
       phi_id = the_model_description%order_of_reservior_higher_order_bc(1)      ! take as N 
   

       !! 9.2.1. put coeff of \ddot{ theta }_N in alpha_coeffs(1) : ( 1- a_N^2 ) / c^2
   
       alpha_coeffs(1) = 1.0_r_kind - the_model_description%reservior_incident_coefficient(phi_id + 1) ** 2  ! ( 1- a_N^2 )
       alpha_coeffs(1) = alpha_coeffs(1) / ( c ** 2 ) ;
       
       !! 9.2.2. put coeff of \ddot{ theta }_{N+1} in alpha_coeffs(3) :  1 + a_N^2 )
       alpha_coeffs(2) = 1.0_r_kind + the_model_description%reservior_incident_coefficient(phi_id + 1) ** 2  ! ( 1+ a_N^2 )
       alpha_coeffs(2) = alpha_coeffs(2) / ( c ** 2 ) ;
   
       !! 9.2.3. put coeff of \dot{ psi } in alpha_coeffs(2) : 2 * a_N /c
       alpha_coeffs(3) = 2.0_r_kind * the_model_description%reservior_incident_coefficient(phi_id + 1) / c ;
   
       
       !! This is added to psi term
        
       row = ( phi_id + 1 ) * n_bc_nodes ;
        
       !! 9.3. term  phi_N
       col =  phi_id  * n_bc_nodes ;
        
       elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ; 
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) =  local_stiffness_matrix ;
        
       if( special_node_index > 0 ) then
           elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = local_damping_matrix ;
       end if
       
       !! 9.4. term  psi
       col =  ( phi_id + 1 )  * n_bc_nodes ;
   
       elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(3) * local_mass_matrix ;
       
       !! 9.5. term  phi_{N+1}
       col =  ( phi_id + 2 )  * n_bc_nodes ;
        
       elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(2) * local_mass_matrix ; 
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) =  local_stiffness_matrix ;
        
       if( special_node_index > 0 ) then
           elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = local_damping_matrix ;
       end if
          
   end if
   
   
   
!!---------------------------------------------------------------------!!
!!
!!  Step 10. Add equations containing  psi.
!!          This equation appears only if M > 0, as it uses first equation of \phi_{N+j},
!!          The term phi_{N+2} appears only if M > 1, i.e. if N+2 < N + M+ 1 --> M > 1
!!---------------------------------------------------------------------!!
   
   if( the_model_description%order_of_reservior_higher_order_bc(2) > 0 ) then
     
       !! 10.1. set coefficients used in this equation
       !!------------------------------------------------------------!!
       !! set alpha(1) = c_bar^2 = 1 / c^2 ;
       !!     alpha(2) =  - 2 * b_bar_1  =  -2 * b_1 / c ;
       !!     alpha(3) = b_bar_1^2   = b_1^2 / c^2 ;
       !!------------------------------------------------------------!!
       
       !! 10.1.1. put  1/c^2
       
       alpha_coeffs(1) = 1.0_r_kind / ( c ** 2 ) ;  
       
       !! 10.1.2. put coeff of  psi in alpha_coeffs(1) : - 2 * b_1 / c
       alpha_coeffs(2) = -2.0_r_kind * the_model_description%reservior_evanescent_coefficient(1) / c ;  

       !! 10.1.3. put coeff of theta_{N+1} in alpha_coeffs(2) : b_1^2 /c^2
       alpha_coeffs(3) = ( the_model_description%reservior_evanescent_coefficient(1) ** 2 ) / ( c ** 2 ) ;
   
       
       phi_id = the_model_description%order_of_reservior_higher_order_bc(1) ;     ! take as N
       
       !! This is added phi_{N+1} term
        
       row = ( phi_id + 2 ) * n_bc_nodes ;
        
       !! 10.2. term  psi
       col =  ( phi_id + 1 ) * n_bc_nodes ;
        
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(2) * local_mass_matrix ; 
        
       
       !! 10.3. term  phi_{N+1}
       col = col + n_bc_nodes ;
        
       elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ; 
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) =  alpha_coeffs(3) * local_mass_matrix  &
                                                                 & + local_stiffness_matrix ;

       if( special_node_index > 0 ) then
           elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = local_damping_matrix ;
       end if
       
       !! 10.4. term  phi_{N+2}
       
       if( the_model_description%order_of_reservior_higher_order_bc(2) > 1 ) then
           
           col = col + n_bc_nodes ;
           elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ; 
           elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = -alpha_coeffs(3) * local_mass_matrix  &
                                                                     & + local_stiffness_matrix ;

           if( special_node_index > 0 ) then
               elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = local_damping_matrix ;
           end if
       
       end if 
       
       
   end if !! M > 0 
    
   
   
!!---------------------------------------------------------------------!!
!!
!!  Step 11. Add equations governing: phi_{N+j , N+j + 1 , N+j+2} : for k = 2 , M
!!  
!!  Also N+j+2 <=  N+M+1 --> j <= M-1, and j>= 1 , so  1 <= j <= M-1
!!  For the last term \phi_{N+j+2} does not appear.
!!
!!---------------------------------------------------------------------!!  

   do  phi_id = 1 , the_model_description%order_of_reservior_higher_order_bc(2)-1
       
       !! 11.1. set coefficients used in this equation
       !!------------------------------------------------------------!!
       !! set alpha(1) = b_bar_{j+1} * c_bar^2 = b_{j+1} / c^3 ;
       !!     alpha(2) = eta_{j} * c_bar^2  = ( b_j + b_{j+1} ) / c^3  ;
       !!     alpha(3) = b_bar_j * c_bar^2 = b_j / c^3 ;
       !!     alpha(4) = b_bar_{j+1} = b_{j+1} / c ;
       !!     alpha(5) = eta_{j}   = ( b_j + b_{j+1} ) / c  ;
       !!     alpha(6) = b_bar_j  = b_j / c ;
       !!     alpha(7) = -b_bar_{j+1} * b_bar_j^2 = - b_{j+1} * b_j^2 / c^3 ;
       !!     alpha(8) = b_bar_j * b_bar_{j+1} * eta_j = b_j * b_{j+1} *( b_j + b_{j+1} ) / c^3 ; 
       !!     alpha(9) = -b_bar_j * b_bar_{j+1}^2 = - b_{j+1}^2 * b_j / c^3 ; 
       !!
       !!---------------------------------------------------------------!!
   
       
       !! 11.1. form coefficients in equation
       
       !! 11.1.1. put  coeff of b_bar_{j+1} * c_bar^2 = b_{j+1} / c^3
       alpha_coeffs(1) = the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) / ( c ** 3 ) ;
        
       !! 11.1.2. put  eta_{j} * c_bar^2  = ( b_j + b_{j+1} ) / c^3  in alpha_coeffs(2)
       alpha_coeffs(2) = the_model_description%reservior_evanescent_coefficient( phi_id )  &
                    &  + the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) ; ! b_{j-1} + b_j
       
       alpha_coeffs(2) = alpha_coeffs(2) / ( c ** 3 ) ;
       
       !! 11.1.3. put  b_bar_j * c_bar^2 = b_j / c^3  in alpha_coeffs(3)                      
       alpha_coeffs(3) = the_model_description%reservior_evanescent_coefficient( phi_id ) / ( c ** 3 ) ;
       
       !! 11.1.4. put  b_bar_{j+1} = b_{j+1} / c  in alpha_coeffs(4)                    
       alpha_coeffs(4) = the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) / c ;
       
       !! 11.1.5. put  eta_{j}   = ( b_j + b_{j+1} ) / c  in alpha_coeffs(5)                    
       alpha_coeffs(5) = ( the_model_description%reservior_evanescent_coefficient( phi_id ) +  &
                         & the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) ) / c ;
    
       
       !! 11.1.6. put  b_bar_j  = b_j / c  in alpha_coeffs(6)   
       alpha_coeffs(6) = the_model_description%reservior_evanescent_coefficient( phi_id ) / c ;
       
       !! 11.1.7. put  -b_bar_{j+1} * b_bar_j^2 = - b_{j+1} * b_j^2 / c^3   in alpha_coeffs(7)   
       alpha_coeffs(7) = -the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) * &
                       &  the_model_description%reservior_evanescent_coefficient( phi_id ) ** 2 ;
                       
       alpha_coeffs(7) =  alpha_coeffs(7) / ( c ** 3 ) ;
       
       !! 11.1.8. put  b_bar_j * b_bar_{j+1} * eta_j = b_j * b_{j+1} *( b_j + b_{j+1} ) / c^3    in alpha_coeffs(8)   
       alpha_coeffs(8) = the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) + &
                       & the_model_description%reservior_evanescent_coefficient( phi_id )  ;
                       
       alpha_coeffs(8) = alpha_coeffs(8) * the_model_description%reservior_evanescent_coefficient( phi_id ) ;
       alpha_coeffs(8) = alpha_coeffs(8) * the_model_description%reservior_evanescent_coefficient( phi_id + 1 ) ;          
       alpha_coeffs(8) = alpha_coeffs(8) / ( c ** 3 ) ;
                     
       !! 11.1.9. put  -b_bar_j * b_bar_{j+1}^2 = - b_{j+1}^2 * b_j / c^3   in alpha_coeffs(8)   
       alpha_coeffs(9) = -the_model_description%reservior_evanescent_coefficient( phi_id ) * &
                       &  the_model_description%reservior_evanescent_coefficient( phi_id + 1) ** 2 ;
                       
       alpha_coeffs(9) =  alpha_coeffs(9) / ( c ** 3 ) ;
           
                    
       !! This is added phi_{N+j+1} term: 1(P) +N(phi_j) + 1( psi) + 1( phi_{N+1} ) + j
       row = ( the_model_description%order_of_reservior_higher_order_bc(1) + 2 + phi_id) * n_bc_nodes ;
    
        
       !! 11.2. term  phi_{N+j}: 
       col =  row - n_bc_nodes ;
        
       elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(1) * local_mass_matrix ; 
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(7) * local_mass_matrix  &
                                                                 & + alpha_coeffs(4) * local_stiffness_matrix ;

       if( special_node_index > 0 ) then
           elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(4) * local_damping_matrix ;
       end if
       
       !! 11.3. term  phi_{N+j+1}: 
       col =  col + n_bc_nodes ;
        
       elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(2) * local_mass_matrix ; 
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(8) * local_mass_matrix  &
                                                                 & + alpha_coeffs(5) * local_stiffness_matrix ;

       if( special_node_index > 0 ) then
           elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(5) * local_damping_matrix ;
       end if
       
       
       if( phi_id .eq. the_model_description%order_of_reservior_higher_order_bc(2)-1 ) cycle ;
       
       
       !! 11.3. term  phi_{N+j+2}: 
       col =  col + n_bc_nodes ;
        
       elemental_mass_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(3) * local_mass_matrix ; 
       elemental_stiffness_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(9) * local_mass_matrix  &
                                                                 & + alpha_coeffs(6) * local_stiffness_matrix ;

       if( special_node_index > 0 ) then
           elemental_damping_matrix( row +1:row +3 , col +1:col +3 ) = alpha_coeffs(6) * local_damping_matrix ;
       end if
       
         
       
   end do ! ( phi_id : 1 , M-1 )
   
   !! send to global matrices
   call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , num_var * n_bc_nodes , & 
                               & col_variable_set , num_var * n_bc_nodes  ) ;
                               
   call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , num_var * n_bc_nodes , & 
                               & col_variable_set , num_var * n_bc_nodes  ) ;
   
   call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , num_var * n_bc_nodes , & 
                               & col_variable_set , num_var * n_bc_nodes  ) ;
   
   
   !!  add force term in P^i equation
   if( special_node_index > 0 ) then
       call  add_to_resultant_force( elemental_force_vector , row_variable_set , num_var * n_bc_nodes ) ;
   end if
   

  end subroutine HW_model_reservior_truncation_boundary
   
   
 
 
 
   
   
  
  end submodule smod_res_trunc_boundary




