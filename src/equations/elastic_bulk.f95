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
!!  This submodule evaluates FEM calculation in bulk of elastic media, i.e. dam and foundation.
!!  No boundary condition is taken into account here.
!!
!!===============================================================================!!

  
  submodule ( mod_fem ) smod_elastic_bulk
  
    use  mod_utils , only : r_kind
    implicit none
    
  
    contains
    
    


!!====================================================================!!
!!
!!  The fem analysis in bulk of elastic media in linear model
!!  Here equations related to bulk of elements are evaluated.
!!  No boundary conditions is applied in this subroutine.
!!
!!  The equation is, here,
!!            M \ddot{ u } + K u = -\bar{M} * J
!!
!!  Mass matrix : M_{I,J} = \int( rho N_IN_J ) I_{num_dim) : submatrices are block-diagonal
!!  Stiffness   : K_{IJ}  = \int( (dN)_I * D * (dN)_J )
!!
!! Standard in defining variable set:
!! Every translational dof is described by a single state function in this model.
!! Therefore to every node num_dim variable is assigned and matrices are formed.
!! Then var_set is defined as follows:
!! Define:
!! v_1 : u_x_1 , u_y_1 , u_z_1 ,  ! translational dof of node 1
!! v_2 : u_x_2 , u_y_2 , u_z_2 ,  ! translational dof of node 2
!! v_k : u_x_k , u_y_k , u_z_k ,  ! translational dof of node k
!!
!! Therefore, variable set is V = [ v_1 , v_2 , v_3 , ... , v_k ]
!!
!! If num_dim == 2 then the components along z direction is not involved.
!! If a node i is absent then v_i = 0 = [0 , 0 , 0 ]
!! If r-th direction of a node is fixed then v_i = [ 0 , ...] , when r = x
!! If varset(i) = 0 then nothing is done.
!! if varset(i) \neq 0 then the information in M(i,:) and K(i,:) should
!! be added to global matrix that corresponds to global var set V( varset(i)).
!!
!!====================================================================!!
   
   subroutine basic_model_elastic_bulk( elem_id , is_drm_calc )
      
      use mod_physics  , only : element_property , material , rigidity_matrix, is_complete_mass_model  
      use mod_physics  , only : dam_jinertia , add_to_global_matrices , add_to_resultant_force
      use mod_physics  , only : dead_weight_force , is_time_domian , gravitational_constant, drm_layer_ndyn
      use mod_geometry , only : element_matrix, degree_of_freedom , num_dim 
      use mod_fem      , only : set_gauss_points , set_rigidity_matrix
      implicit none
      
      integer , intent( in ) :: elem_id 
      logical , intent( in ) :: is_drm_calc
      integer :: k , n, r , s , t , n_gauss_point , region_id , row , col , num_var , a , b
      real( kind = r_kind ) :: rho , d_volume 
        
                   
      !----------------------------------------------------------------!
      !! set initial values for elemental matrices
   
      elemental_mass_matrix(:,:)      = 0.0_r_kind
      
      elemental_stiffness_matrix(:,:) = 0.0_r_kind
      elemental_damping_matrix(:,:)   = 0.0_r_kind 
      
      elemental_force_vector(:,:)     = 0.0_r_kind

      !----------------------------------------------------------------!
      !! set physical property of the element
      
      rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
      
      region_id     = element_property( elem_id , 1 )
      n_gauss_point = element_property( elem_id , 3 )
      
      call  set_gauss_points( n_gauss_point )
      call  set_rigidity_matrix( elem_id , region_id ) 
      
      !---------------------------------------------------------------!
      ! set suitable variables set: every dof has 1 state function in this model
     
      row_variable_set = 0
      num_var = 0 
   
    do  k = 1 , the_element%num_node_in_element
      
        n = element_matrix( elem_id , k )
         
        if( n .eq. 0 ) cycle                                           ! if node is absent.
                                                
        if( degree_of_freedom( n , 1) .eq. 0 ) cycle                   ! if node has no dof.
            
        row = ( k - 1 ) * num_dim
                                                           
        do r = 1 , num_dim                                              ! get state functions of translational dof
        
            if( degree_of_freedom( n , 2 * r ) .eq. 0 ) cycle            ! r-th direction of translation vector is fixed.   
           
            row_variable_set( row + r ) = degree_of_freedom( n , 2 * r )
           
        end do ! over r           
         
    end do ! over k
     
     num_var = the_element%num_node_in_element * num_dim 
     
     col_variable_set = row_variable_set
          
     !!---------------------------------------------------------------!!
     !!  2D element : calculate elemental mass and stiffness matrices
      
      
      if( num_dim .eq. 2 ) then
                 
          do r = 1 , n_gauss_point
             do s = 1 , n_gauss_point
             
                call the_element%get_shape_function( elem_id , gauss_points(r) , gauss_points(s)  )
                
                call the_element%get_derivative_shape_function( elem_id , gauss_points(r) , gauss_points(s) )
                
                d_volume = weight_of_gauss_point(r) * weight_of_gauss_point(s) *  &
                        &  abs( the_element%determinant_of_jacobian )
                          
                
                do  k = 1 , size( the_element%shape_vector )            !! do sum over all nodes
                    
                    !! first form dN_I : derivative of shape function
                    dN = 0.0_r_kind
                    dN( 1 , 1 ) = the_element%derivative_shape_matrix( 1 , k )    ! D_x N_k  
                    dN( 2 , 2 ) = the_element%derivative_shape_matrix( 2 , k )    ! D_y N_k 
                    
                    
                    dN( 1 , 3 ) = dN( 2 , 2 )                           ! D_y N_k
                    dN( 2 , 3 ) = dN( 1 , 1 )                           ! D_x N_k
                    
                    
                    !! in 2D element, dN is 2*3 and rigidity matrix is 3 * 3, so dN_rigidity becomes 2 by 3
                    dN_rigidity = matmul( dN , rigidity_matrix ) 
                    
                    row = ( k - 1 ) * num_dim                            ! row index
                    
                    
                    !! calculate PBEL only for dam elements
                    if( is_time_domian .and. ( dam_jinertia .eq. 0 ) .and. (region_id .eq. 1) .and. &
                        &  ( .not.  is_drm_calc ) ) then
                        if( row_variable_set(2*k) /= 0 ) then
                            dead_weight_force(row_variable_set(2*k)) = dead_weight_force(row_variable_set(2*k)) &
                        &   - rho * gravitational_constant * the_element%shape_vector(k) * d_volume ;
                        end if
                    end if
                    
                    do n = 1 , size( the_element%shape_vector )
                       
                       col = ( n - 1 ) * num_dim
                       
                       
                       !! mass matrix : submatrix is diagonal: a multiple of I_{num_dim}
                       do  a = 1 , num_dim
                           elemental_mass_matrix( row + a , col + a ) =  &
                        &  elemental_mass_matrix( row + a , col + a ) +  &
                        &  the_element%shape_vector(k) * the_element%shape_vector(n) * d_volume * rho
                       end do
                      
                       
           
                       !! stiffness matrix , now form dN_J
                       dN = 0.0_r_kind
                       dN( 1 , 1 ) = the_element%derivative_shape_matrix( 1 , n )    ! D_x N_n   
                       dN( 2 , 2 ) = the_element%derivative_shape_matrix( 2 , n )    ! D_y N_n 
                       
                       
                       dN( 1 , 3 ) = dN( 2 , 2 )                               ! D_y N_n
                       dN( 2 , 3 ) = dN( 1 , 1 )                               ! D_x N_n
                                   
                       dN_D_dN =  matmul( dN_rigidity , transpose( dN ) )      ! a matrix of size num_dim \times num_dim
                       
                       !! stiffness matrix is not submatrix diagonal
                       do a = 1 , num_dim
                          do b = 1 , num_dim
                          
                             elemental_stiffness_matrix( row + a  , col + b  )  = &
                           & elemental_stiffness_matrix( row + a  , col + b  )  + &
                           & dN_D_dN( a , b ) * d_volume
                              
                          end do ! over b
                       end do ! over a
                      
                     
                     
                    end do ! over n
                end do ! over k
                  
             end do ! over s
          end do ! over r
          
          
      else  !! 3D problem
      
          do r = 1 , n_gauss_point
             do s = 1 ,  n_gauss_point
                do t = 1 , n_gauss_point
                
                   call the_element%get_shape_function( elem_id , gauss_points(r) , &
                                       &   gauss_points(s) , gauss_points(t) )
                   
                   call the_element%get_derivative_shape_function( elem_id , gauss_points(r) ,  &
                                                    & gauss_points(s) , gauss_points(t) )
                   
                   d_volume = weight_of_gauss_point(r) * weight_of_gauss_point(s) *  &
                            & weight_of_gauss_point(t) * abs( the_element%determinant_of_jacobian )  
                   
                    
                   !! do sum over all nodes
                
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
                       
                       !! in 3D element, dN is 3 by 6 and rigidity matrix is 6 by 6 so dN_rigidity is 3 by 6
                       dN_rigidity = matmul( dN , rigidity_matrix ) 
                       
                       row = (k-1) * num_dim 
                       
                       
                       !! calculate PBEL
                        if( is_time_domian .and. ( dam_jinertia .eq. 0 ) .and. (region_id .eq. 1) .and.  &
                            &  ( .not.  is_drm_calc ) ) then
                            if( row_variable_set(3*k) /= 0 ) then
                                dead_weight_force(row_variable_set(3*k)) = dead_weight_force(row_variable_set(3*k)) &
                            &   - rho * gravitational_constant * the_element%shape_vector(k) * d_volume ;
                            end if
                        end if
                    
                    
                        do n = 1 , size( the_element%shape_vector )
                        
                            col = (n-1) * num_dim

                          !! mass matrix : submatrix is diagonal: a multiple of I_{num_dim}
                          do  a = 1 , num_dim
                              elemental_mass_matrix( row + a , col + a ) =  &
                           &  elemental_mass_matrix( row + a , col + a )  +  &
                           &  the_element%shape_vector(k) * the_element%shape_vector(n) * d_volume * rho
                          end do
                         
                       
                     
                       
                          !! stiffness matrix, first form dN_J
                       
                          dN = 0.0_r_kind
                          dN( 1 , 1 ) = the_element%derivative_shape_matrix( 1 , n )  ! D_x N_k  
                          dN( 2 , 2 ) = the_element%derivative_shape_matrix( 2 , n )  ! D_y N_k
                          dN( 3 , 3 ) = the_element%derivative_shape_matrix( 3 , n )  ! D_z N_k
                       
                          dN( 1 , 4 ) = dN( 2 , 2 )                             ! D_y N_k 
                          dN( 2 , 4 ) = dN( 1 , 1 )                             ! D_x N_k 
                          dN( 2 , 5 ) = dN( 3 , 3 )                             ! D_z N_k 
                          dN( 3 , 5 ) = dN( 2 , 2 )                             ! D_y N_k 
                          dN( 1 , 6 ) = dN( 3 , 3 )                             ! D_z N_k 
                          dN( 3 , 6 ) = dN( 1 , 1 )                             ! D_x N_k 
                       
                            
                            
                          dN_D_dN =  matmul( dN_rigidity , transpose( dN ) ) 
                            
                          !! stiffness matrix is not submatrix diagonal
                          do a = 1 , num_dim
                             do b = 1 , num_dim
                          
                                elemental_stiffness_matrix( row + a  , col + b  )  = &
                             &  elemental_stiffness_matrix( row + a  , col + b  )  + &
                             &  dN_D_dN( a , b ) * d_volume
                              
                             end do ! over b
                          end do ! over a
                       
                       
                       end do ! over n
                   end do ! over k
                
                
                end do ! over t
             end do ! over s
          end do !  over r
          
          
    end if
     
    
    if( is_drm_calc ) return ; !! only elemental matrices is needed
     
    
    !! put elemental mass matrix in global mass matrix
    !! foundation mass involves only if the mass model is complete ( Jinertia = 1 )
     
    if( region_id .eq. 2 ) then
        if( is_complete_mass_model ) then
          
            call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                        & num_var , col_variable_set , num_var )
                                    
        end if
    else 
            call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , & 
                                        & num_var , col_variable_set , num_var )
                                    
    end if
     
     
    !! put elemental stiffness matrix in global stiffness matrix
     
    call  add_to_global_matrices( elemental_stiffness_matrix , 2 , row_variable_set , & 
                                & num_var , col_variable_set , num_var )
                                   
     
     
    ! the founation mass does not involve in force vector:
    ! Edit after DRM Model:
    ! If model is not DRM then M_bar * J of dam is calculated,
    ! Otherwise, M_bar * J for DRM layer is calculated.
    
    if( region_id .eq. 2 ) return ;
    
    if( drm_layer_ndyn /= 0 ) return  !! if there is DRM layer
        
     
    do k = 1 , the_element%num_node_in_element
        row = ( k - 1 ) * num_dim
         
        do n = 1 , num_dim
            do  r = 1 , the_element%num_node_in_element
                s = ( r - 1 ) * num_dim
                   elemental_force_vector( row + n , n ) = elemental_force_vector( row + n , n ) - &
                &  elemental_mass_matrix( row + n , s + n ) 
             
            end do ! over r
        end do ! over n
    end do ! over k
      
    call  add_to_resultant_force( elemental_force_vector , row_variable_set , num_var )
    
    
        

    
    
    
    
     
   end subroutine basic_model_elastic_bulk
   
 
 
!!==========================================================================================
!!
!! Surface pressure force in dam-reservoir interface: in dynamic analysis
!!
!!========================================================================================== 
   
    subroutine elastic_surface_pressure_force( elem_id  , bc_id )
      
        use mod_physics  , only : element_property , boundary_load , surface_pressure_force    
        use mod_geometry , only : element_matrix, degree_of_freedom , num_dim, node 
        use mod_fem      , only : set_gauss_points 
        use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face 
       
        implicit none
        integer , intent( in ) :: elem_id , bc_id
      
        integer :: k , j , n, r , s  , d  
        integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point, dir, ktype
        real( kind = r_kind ) ::  d_volume , fixed_value , pr , xiref , pr1 , xinew 
        
        
        !!----------------------------------------------------------------!
        !! set initial values
     
        indexed_set = 0

        !!----------------------------------------------------------------!
        !! get physical property of the element

        pr     = boundary_load( bc_id )%strength  ;                     ! PR
        xiref  = boundary_load( bc_id )%reference_height ;              ! XIREFF
        dir    = boundary_load( bc_id )%direction ;                     ! LDIR 
        n_face = boundary_load( bc_id )%face_id   ;                     ! NFACE
        ktype  = boundary_load( bc_id )%k_type    ;                     ! KTYPE 
        
        !!----------------------------------------------------------------!
        !! set gauss points
     
        n_gauss_point = element_property( elem_id , 3 )
      
        call  set_gauss_points( n_gauss_point )
     
        n_face = boundary_load( bc_id )%face_id
        n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
     
     
        call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
        fixed_value = real( direction_value(2) , kind = r_kind )           !! fixed_value = xi or eta or zeta at = +1 or -1
     
     
        
        !!----------------------------------------------------------------!
        !! set suitable variable set.
        row_variable_set = 0
        num_var = 0 
   
        do  k = 1 , the_element%num_node_in_element
            n = element_matrix( elem_id , k )
         
            if( n .eq. 0 ) cycle                                        ! if node is absent.
                                                
            if( degree_of_freedom( n , 1) .eq. 0 ) cycle                ! if node has no dof.
            
            j = ( k - 1 ) * num_dim
                                                           
            do  r = 1 , num_dim                                         ! get state functions of translational dof
                row_variable_set( j + r ) = degree_of_freedom( n , 2 * r )
            end do ! over r           
         
        end do ! over k
     
        num_var = the_element%num_node_in_element * num_dim
     
     
        !------------------------------------------------------------------!

        if( num_dim .eq. 2 ) then                                       ! there is only one direction to do integration
                 
            do  r = 1 ,  n_gauss_point 
              
                if( direction_value(1) .eq. 1 ) then                        !! find shape function on face
               
                    call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r)  )      ! xi = fixed
                    d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) ) 
                    d_volume = d_volume  * weight_of_gauss_point(r)  
                else
               
                    call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value  )      ! eta  = fixed
                    d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value ) 
                    d_volume =  d_volume * weight_of_gauss_point(r)   
                end if 
            
                call the_element%get_normal_at_face( elem_id , direction_value(1)  ,  &
                                                   & indexed_set(1:n_bc_nodes) , is_outward = .false. )     
            
           
                pr1 = pr ;
                if( ktype /= 1 ) then
                    xinew = 0.0_r_kind
                    do  k = 1 , the_element%num_node_in_element
                        n = element_matrix( elem_id , k )
                        if( n .eq. 0 ) cycle                            ! if node is absent.
                        xinew = xinew + node( n , dir ) * the_element%shape_vector( k ) ;
                    end do
                    pr1 = pr * (  xiref - xinew ) ;
                    if( pr1 .le. 0.0_r_kind ) pr1 = 0.0_r_kind ;
	
                end if   
            
            
           
                do  k = 1 , the_element%num_node_in_element  
                    do  j = 1 , num_dim
                        n = row_variable_set( (k-1)*num_dim + j ) ;
                        if( n .eq. 0 ) cycle ;
                        surface_pressure_force(n) = surface_pressure_force(n) +  &
                    &   the_element%shape_vector( k ) * the_element%normal_at_face(j) * pr1 *  d_volume ;
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
                
     call the_element%get_normal_at_face( elem_id, direction_value(1), indexed_set(1:n_bc_nodes), is_outward = .false. )
        
                pr1 = pr ;
                if( ktype /= 1 ) then
                    xinew = 0.0_r_kind
                    do  k = 1 , the_element%num_node_in_element
                        n = element_matrix( elem_id , k )
                        if( n .eq. 0 ) cycle                            ! if node is absent.
                        xinew = xinew + node( n , dir ) * the_element%shape_vector( k ) ;
                    end do
                    pr1 = pr * (  xiref - xinew ) ;
                    if( pr1 .le. 0.0_r_kind ) pr1 = 0.0_r_kind ;
                    	
                end if   
            
            
                do  k = 1 , the_element%num_node_in_element  
                    do  j = 1, num_dim
                        n = row_variable_set( (k-1)*num_dim + j ) ;
                        if( n .eq. 0 ) cycle ;
                        surface_pressure_force(n) = surface_pressure_force(n) +  &
                    &   the_element%shape_vector( k ) * the_element%normal_at_face(j) * pr1 *  d_volume ;
                    end do 
                    
                end do ! over k
                  
                
             end do ! over s
          end do ! over r
          
        end if
        
        
  
    end subroutine elastic_surface_pressure_force
   
   
   
   
   
  
  end submodule smod_elastic_bulk
