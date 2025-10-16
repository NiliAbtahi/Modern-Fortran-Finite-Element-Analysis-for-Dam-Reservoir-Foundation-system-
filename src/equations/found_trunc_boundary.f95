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
!!  This submodule evaluates FEM calculation in foundation truncation boundary.
!! 
!!
!!===============================================================================!!

  
  submodule ( mod_fem ) smod_found_trunc_boundary
  
    use  mod_utils , only : r_kind
    implicit none
    
    real( kind = r_kind )   , allocatable :: a(:) , b(:) 
  
    contains
    
    

  
!!====================================================================!!
!!
!!  The fem analysis in foundation_truncation_boundary in basic model:             Lysmer model.
!!  i.e. every component of translational dofs has only one state function :       Lysmer BC
!!  The equation is 
!!                 M \ddot{u} + C \dot{u} + K u = F
!!
!!  Here the damping term C \dot{u} is evaluated where
!!
!!    C = rho * \int{ N * C_{Lys} * N^T } , on boundary face.
!!    C_{I,J} = rho * \int{ N_I C_{Lys} N_J }
!!
!! The Lysmer rotation matrix is obtained from outward normal vector
!! as the angle of vectors is specified by outward normal.
!! Some conditions are imposed to get rotation matrix.
!!
!! 1. The rotation matrix should map the outward normal to e_x (unit vector along x-direction),
!!    that is R * n = e_x.
!!    This condition comes from the extension of Lysmer formulation from 1D to 2D or 3D.
!!
!! 2. The Lysmer rotation matrix is indeed a rotation, hence an orthogonal matrix with +1 determinant.
!!    This implies that R^T R = RR^T = I , and det(R) = + 1, as an element of orthogonal group.
!!
!! These two conditions are enough to specify the rotation matrix.
!!
!! To see this let n be the normal vector of the boundary face.
!! Also let R = [v1 ; v2 ; v3 ] be rotation matrix where v1 , v2 , v3 are rows of the matrix.
!!
!! Then we should have  R * n = e_x  or v1.n = 1 , v2.n = 0 , v3.n = 0 ,       (1)
!!
!! Also the condition R^T * R = I implies v1.v1 = v2.v2 = v3.v3 = 1,           (2)    
!!                                        v1.v2 = v1.v2 = v2.v3 = 0            (3)
!! 
!! The last condition is about determinant,   det(R) = + 1  ,                  (4).
!!
!! This means that v1 and v2 and v3 are orthonormal vectors ( orthogonal unit vectors )
!! 
!! Now from equation (2) we saw that \{ v1 , v2 , v3 \} are orthogonal vectors and 
!! equation (1) states that v2 and v3 are orthogonal to n, hence, in 3D space , v1 is parallel to n1.
!!
!! This means that v1 = k n for some scalar factor k. From equation (1) we have v1.n= 1 --> k = 1/norm(n1)^2
!! This means that v1 is identical to unit normal vector n.
!!
!! In 2D problem the normal unit vector can be written as, in polar coordinate, n = ( cos(theta) , sin(theta) ).
!! So R(1,:) = v1 = n1 = ( cos(theta) , sin(theta) )
!! The second condition is, from equation (1) , v2.v1 = v2.n = 0 
!! There is two possibilities: v2 = ( - sin(theta) , cos(theta) )  or v2 = ( + sin(theta) , - cos(theta) )
!! The second choice makes det(R) = -1 while the first choice implies det(R) = +1.
!! So the first possibility is chosen and we get the following matrix in 2D.
!! 
!!    R = [ v1 ; v2 ] = [ cos(theta) , sin(theta) ;
!!                       -sin(theta) , cos(theta) ];
!! or, R = [ v1 ; v2 ] = [ n(1) , n(2) ; 
!!                        -n(2) , n(1) ] where is unit outward normal vector to face.
!! This is a well-known result to you.
!!
!! In 3D again the equation (1) states that v1 = n1: the unit outward normal to the boundary face.
!! The boundary is 2D, and v2 and v3 are orthogonal to v1=n, hence v2 and v3 lie in boundary face.
!!
!! To get them we first find two tangent vectors in boundary face from a basic rule.
!!    1. we define position vector as r(xi , eta , zeta ) = x(xi , eta , zeta) e_x + y(xi , eta , zeta) e_y + z(xi , eta , zeta) e_z
!!
!!    2.1. If boundary face is xi = const, then we put v2 = ( dr / d eta) d eta , v3 = ( dr / d zeta) d zeta 
!!    2.2. If boundary face is eta = const, then we put v2 = ( dr / d xi) d xi , v3 = ( dr / d zeta) d zeta 
!!    2.3. If boundary face is zeta = const, then we put v2 = ( dr / d xi) d xi , v3 = ( dr / d eta) d eta
!!   
!!    The two vectors \{v2 , v3 \} should be 
!!    3.1. orthogonal to v1 = n , i.e. v2.n = v3.n = 0
!!    3.2. orthogonal to each other, i.e. v2.v3 = 0
!!    3.3. uint vectors, i.e. v2.v2 = v3.v3 = 1
!! 
!!   This step is done by Gram-Schmit orthogonalization process.
!!   With these conditions we have R^T * R = R * R^T = I, hence R belongs to orthogonal group O(3).
!!
!!   The last condition comes from positive determinant, i.e. det(R) = 1 , or R belongs to special orthogonal group SO(3).
!!   
!!   4. To impose the last condition we just find determinant of R. 
!!      4.1. If it is +1 there is nothing to do, i.e. R = [ n ; v2 ; v3]
!!      4.2. If it is -1 then we just need to replace the second and third rows, i.e. R = [n ; v3 ; v2]
!!      as changing the rows multiplies determinant by -1.
!!
!!   Therefore the Lysmer rotation matrix in 2D and 3D are found based on the 4 conditions above.
!! 
!!====================================================================!!
   
  module subroutine basic_model_foundation_truncation_boundary( elem_id , bc_id )
   
     use mod_physics  , only : element_property , material , boundary_load , add_to_global_matrices
     use mod_geometry , only : element_matrix, degree_of_freedom , num_dim 
     use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face 
      use  mod_fem     , only : set_gauss_points
       
     implicit none
     integer , intent( in ) :: elem_id , bc_id
      
     integer :: k , j , n, r , s  , d , q1 , q2  
     integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point
     real( kind = r_kind ) ::  d_volume , rho , nu , E , cp , cs , fixed_value , coeff 
     
        
     !!----------------------------------------------------------------!
     !! set initial values for damping matrix
     
     elemental_damping_matrix = 0.0_r_kind
     indexed_set = 0
     
     
     !!----------------------------------------------------------------!
     !! get physical property of the element
      
     E   = material( element_property( elem_id , 2) )%elastic_module  
     nu  = material( element_property( elem_id , 2) )%poisson_ratio  
     rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind
     
     
     
     cp = sqrt( E * ( 1.0_r_kind - nu ) / ( 1.0_r_kind + nu ) / ( 1.0_r_kind - 2.0_r_kind * nu ) / rho )
     cs = sqrt( E / ( 2.0_r_kind + 2.0_r_kind * nu ) / rho )
     
     
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
         
     do k = 1 , n_bc_nodes
        
        n = element_matrix( elem_id , indexed_set(k) ) 
            
        if( n .eq. 0 ) cycle                                            ! if node is absent or no dof
        if( degree_of_freedom( n , 1 ) .eq. 0 ) cycle                   ! node has no dof
   
        q1 = ( indexed_set(k) - 1 ) * num_dim
        do  d = 1 , num_dim
            row_variable_set( q1 + d ) = degree_of_freedom( n  , 2 * d ) 
        end do ! over d
        
     end do ! over k
     
     num_var = num_dim * the_element%num_node_in_element
     col_variable_set  = row_variable_set 
     
     !------------------------------------------------------------------!

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
            
            call the_element%get_normal_at_face( elem_id , direction_value(1)  ,  &
                                              &  indexed_set(1:n_bc_nodes) , is_outward = .true. )     
            
           
            !! form rotation matrix
            rotation = 0.0_r_kind
            rotation(1,1) = the_element%normal_at_face(1)
            rotation(1,2) = the_element%normal_at_face(2)
            rotation(2,1) = -rotation(1,2)
            rotation(2,2) =  rotation(1,1)
            
            
            !! form R^T * C_{Lysmer}  , resuse dN for this
            dN = 0.0_r_kind          
            
            do k = 1 , 2
               dN(k, 1 ) = cp * rotation(1 , k)                         ! dN is used tomporary for storing R^T * C
               dN(k, 2 ) = cs * rotation(2 , k)
            end do
            
            
            !! form whole Lysmer matrix : R^T * C * R = dN * R
            C_Lysmer = 0.0_r_kind
                
            do k = 1 , num_dim
               do n = 1 , num_dim
                  do d = 1 , num_dim
                      C_Lysmer(k , n ) = C_Lysmer(k , n ) + dN(k , d ) * rotation(d , n )
                  end do
               end do
            end do 
            
           
           
            do  k = 1 , n_bc_nodes     
                q1 = ( indexed_set(k) -1 ) * num_dim
                     
                do n = 1 , n_bc_nodes
                   q2 = ( indexed_set(n) -1) * num_dim
                   
                   coeff = the_element%shape_vector( indexed_set(k) ) *  d_volume * &   
                         & the_element%shape_vector( indexed_set(n) ) * rho
                   
                   
                   do  d = 1 , num_dim
                      do  j = 1 , num_dim
                          elemental_damping_matrix( q1 + d , q2 + j ) =  &
                       &  elemental_damping_matrix( q1 + d , q2 + j ) +  coeff * C_Lysmer( d , j )
                      end do
                   end do ! over d
                   
                   
                end do ! over n
             end do ! over k
                  
         end do ! over r
          
          
     else  !! 3D problem
      
          do r = 1 , n_gauss_point                                      
             do s = 1 ,  n_gauss_point
             
                rotation = 0.0_r_kind
                if( direction_value(1) .eq. 1 ) then                    !! xi = fixed_value
               
      call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r) , gauss_points(s) ) 
      d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) , gauss_points(s) ) 
      d_volume = abs( d_volume ) * weight_of_gauss_point(r) * weight_of_gauss_point(s)
                    
                   
                    !! first row of Jacobian matrix is zero 1-2-3
                    rotation(2,1:3) = the_element%jacobian_matrix(2,1:3)
                    rotation(3,1:3) = the_element%jacobian_matrix(3,1:3)
                                                
                elseif( direction_value(1) .eq. 2 ) then                !! eta = fixed_value
               
      call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value , gauss_points(s) )   
      d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value , gauss_points(s) ) 
      d_volume = abs( d_volume ) * weight_of_gauss_point(r) * weight_of_gauss_point(s) 
                             
                    !! second row of Jacobian matrix is zero 2-3-1
                    rotation(2,1:3) = the_element%jacobian_matrix(1,1:3)
                    rotation(3,1:3) = the_element%jacobian_matrix(3,1:3)
                    
                else                                                    !! zeta = fixed_value
                       
      call the_element%get_shape_function( elem_id , gauss_points(r) , gauss_points(s) , fixed_value  )
      d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , gauss_points(s) , fixed_value ) 
      d_volume = abs( d_volume ) * weight_of_gauss_point(r) * weight_of_gauss_point(s)
                                   
                    !! third row of Jacobian matrix is zero 3-1-2
                    rotation(2,1:3) = the_element%jacobian_matrix(1,1:3)
                    rotation(3,1:3) = the_element%jacobian_matrix(2,1:3)
                    
                    
                end if 
                
     call the_element%get_normal_at_face( elem_id, direction_value(1), indexed_set(1:n_bc_nodes), is_outward = .true. )
        
                rotation(1,1:3) = the_element%normal_at_face
                
                ! normalize normal vector
                coeff = sqrt( rotation(1,1) ** 2 +  rotation(1,2) ** 2 + rotation(1,3) ** 2 )
                rotation(1,1:3) = rotation(1,1:3) / coeff
                
                
                !! orthogonalize the second vector on the face with normal vector         
                coeff = 0.0_r_kind
                do d = 1 , 3
                   coeff = coeff + rotation( 2 , d) * rotation( 1 , d)
                end do
                rotation( 2 , 1:3) = rotation( 2 , 1:3) - coeff * rotation( 1 , 1:3) ! orthogonalize v2 and n = v1
                
                !! make second vector unit
                coeff = sqrt( rotation(2,1) ** 2 +  rotation(2,2) ** 2 + rotation(2,3) ** 2 )
                rotation(2,1:3) = rotation(2,1:3) / coeff
                
                
                !! orthogoanlize third vector with the two others  
                   
                do k = 1 , 2
                   coeff = 0.0_r_kind
                   do d = 1 , 3
                      coeff = coeff + rotation( 3 , d) * rotation( k , d)
                   end do
                   
                   rotation( 3 , 1:3) = rotation( 3 , 1:3) - coeff * rotation( k , 1:3)
                end do
                
                !! make third vector unit
                coeff = sqrt( rotation(3,1) ** 2 +  rotation(3,2) ** 2 + rotation(3,3) ** 2 )
                rotation(3,1:3) = rotation(3,1:3) / coeff
                
                !! the rotation matrix is now orthogonal, i.e. R^T * R = I , check if determinant is +1 or -1
                !! get determinant
                
                coeff = 0.0_r_kind
                coeff = coeff + rotation( 1 , 1 ) * rotation( 2 , 2 ) * rotation( 3 , 3 )
                coeff = coeff + rotation( 1 , 2 ) * rotation( 2 , 3 ) * rotation( 3 , 1 )
                coeff = coeff + rotation( 1 , 3 ) * rotation( 2 , 1 ) * rotation( 3 , 2 )
                coeff = coeff - rotation( 1 , 3 ) * rotation( 2 , 2 ) * rotation( 3 , 1 )
                coeff = coeff - rotation( 1 , 2 ) * rotation( 2 , 1 ) * rotation( 3 , 3 )
                coeff = coeff - rotation( 1 , 1 ) * rotation( 2 , 3 ) * rotation( 3 , 2 )
                
                if( coeff < 0.0_r_kind ) then                           ! if det(R) = -1 then replace v2 and v3.
                                          
                    dN( 1 , 1:3 ) = rotation( 2 , 1:3 )                 ! used as temp vector
                    rotation( 2 , 1:3 ) = rotation( 3 , 1:3 )
                    rotation( 3 , 1:3 ) = dN( 1 , 1:3 )
                    
                end if
                
                !! the 3D Lysmer rotation matrix is now avilable.
                
                
                !! form R^T * C_{Lysmer}  , resuse dN for this
                dN = 0.0_r_kind          
            
                do k = 1 , 3
                   dN(k, 1 ) = cp * rotation(1 , k)                     ! dN is used tomporary for storing R^T * C
                   dN(k, 2 ) = cs * rotation(2 , k)
                   dN(k, 3 ) = cs * rotation(3 , k)
                end do
            
                !! form whole Lysmer matrix : R^T * C * R = dN * R
                C_Lysmer = 0.0_r_kind
                
                do k = 1 , 3
                   do n = 1 , 3
                      do d = 1 , 3
                         C_Lysmer(k , n ) = C_Lysmer(k , n ) + dN(k , d ) * rotation(d , n )
                      end do
                   end do
                end do 
                
                
                !! form damping matrix from Lysmer boundary condition.
                
                do  k = 1 , n_bc_nodes                                  ! other terms are zero   
                    q1 =  ( indexed_set(k)  - 1 )  * num_dim
                                  
                    do n = 1 , n_bc_nodes                 
                       q2 =  ( indexed_set(n)  - 1 )  * num_dim
                       
                       coeff = the_element%shape_vector( indexed_set(k) ) *  d_volume * &   
                             & the_element%shape_vector( indexed_set(n) ) * rho
                        
                       do  d = 1 , num_dim
                          do  j = 1 , num_dim
                              elemental_damping_matrix( q1 + d , q2 + j ) =  &
                           &  elemental_damping_matrix( q1 + d , q2 + j ) +  coeff * C_Lysmer( d , j )
                          end do
                       end do ! over d
                   
             
                    end do ! over n
                end do ! over k
                  
                
             end do ! over s
          end do ! over r
          
     end if
       
    !! assemble to global matrix
    call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , & 
                                 & num_var , col_variable_set , num_var  ) 
  
     
     
  end subroutine basic_model_foundation_truncation_boundary
   
   
   
  
  
  
  
!!====================================================================!!
!!
!!  The fem analysis in foundation_truncation_boundary in HW model 
!!
!! 
!!====================================================================!!
   
 module subroutine HW_model_foundation_truncation_boundary( elem_id , bc_id )  
   
    use mod_physics  , only : element_property , material , boundary_load , element_type
    use mod_geometry , only : num_dim , degree_of_freedom , element_matrix
    use mod_utils    , only : get_nodes_on_element_face , local_coordinate_on_face
    use mod_physics  , only : the_model_description , add_to_resultant_force , add_to_global_matrices
    use mod_fem      , only : set_gauss_points
    implicit none
 
    integer , intent( in ) :: elem_id , bc_id
    integer :: num_var ,  n_bc_nodes , n_face , direction_value(2) , n_gauss_point , delta_corner(4)
    integer :: i , k , n , r, s , phi_id , special_node_index, var_id, row , col, num_row , num_col, dir
    real( kind = r_kind ) :: rho , E, alpha , mu , nu , kappa , fixed_value , d_volume, nx , ny , c
    real( kind = r_kind ) :: Mmatrix(3,3) , Kmatrix(3,3) , Smatrix(3,3)  , theta_coeff(3) , gamma_coeff(3)
    real( kind = r_kind ) :: B14matrix(3,3) , B23matrix(3,3), ellipitic_param = 0.0_r_kind
    logical  :: is_edge_vertical , auto_coeff = .false. 
    
    
    
    !!---------------------------------------------------
    !! step 0. for simplicty put incident and evanescent coefficients in a and b arrays
    
    
    !! 0.1. deallocate a_j if it is allocated
    if( .not. allocated( a ) )  then
        !! 0.1.1. allocate memory for a_j
        allocate( a( the_model_description%order_of_foundation_higher_order_bc(1) + 1 ) , stat = i )
        if( i /= 0 ) stop 'Error: failed in allocating a-array in HW foundation model.'
        
        !! 0.1.2. store foundation_incident_coefficient in a-array
        do  i = 1 , the_model_description%order_of_foundation_higher_order_bc(1) + 1
            a(i) = the_model_description%foundation_incident_coefficient(i) ;
        end do
    end if

    
    !! 0.2. allocate memory for b_j
    if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
        !! 0.2.1. deallocate b_j if it is allocated
        if( .not. allocated( b ) )  then
            allocate( b( the_model_description%order_of_foundation_higher_order_bc(2)  ) , stat = i )
            if( i /= 0 ) stop 'Error: failed in allocating b-array in HW foundation model.'
            
            !! 0.2.2. store foundation_evanescent_coefficient in b-array
            do i  = 1 , the_model_description%order_of_foundation_higher_order_bc(2)
                b(i) = the_model_description%foundation_evanescent_coefficient(i) ;
            end do
        end if
    end if
    
    
    !!----------------------------------------------------------------!
    !! step 1. get physical property of the element
    
    E   = material( element_property( elem_id , 2) )%elastic_module     ! elastic module
    nu  = material( element_property( elem_id , 2) )%poisson_ratio      ! Poisson ratio
    rho = material( element_property( elem_id , 2) )%special_weight / 9.810_r_kind   ! mass density
    
    mu  = E / ( 2.0_r_kind * ( 1.0_r_kind + nu ) ) ;                    ! shear module
    
    !! get element type of reagion 2 (foundation) : 5 is plane strain and 6 is plane stress
    if( element_type( 2 ) .eq. 5 ) then                                 ! plane strain 
        kappa = E / ( 3.0_r_kind * ( 1.0_r_kind - 2.0_r_kind* nu ) ) ;              
    elseif( element_type( 2 ) .eq. 6 ) then                             ! plane stress
        kappa =( 1.0_r_kind + 2.0_r_kind * nu )* E /( 3.0_r_kind * ( 1.0_r_kind - nu ** 2 ) ) ;   
    else
        stop 'The element type has not defined properly, plane stress or plane strain'
    end if
     
    alpha = kappa - 2.0_r_kind * mu / 3.0_r_kind ;
    
    !! speed of shear wave: cs = sqrt( G /rho )  where G =  E / ( 2  * (1+ nu) ) = mu 
    c   = sqrt(  mu / rho )  ; 
  
    
!!---------------------------------------------------------------------!!
!!
!!  Step 2. Form boundary_shape_matrix and boundary_der_shape_matrix.
!!
!!---------------------------------------------------------------------!!    
     
    !!---------------------------------------------------------------!!
    !! 2.1. get local node index of boundary nodes
     
    indexed_set = 0
     
    n_face     = boundary_load( bc_id )%face_id
    n_bc_nodes = get_nodes_on_element_face( n_face , the_element%node_parent , indexed_set  ) 
    
    
    
    call  local_coordinate_on_face( n_face , direction_value )         !! set value of local coordinates on boundary
    fixed_value = real( direction_value(2) , kind = r_kind )           !! either xi or eta or zeta is fixed at +1 or -1
     
     
    !! check if the current element contains corner nodes
    delta_corner(1:4) = -2 ;
    do  i = 1 , 4
        do  k = 1 , n_bc_nodes
            if( element_matrix( elem_id , indexed_set(k) ) .eq. found_inf_corner_node(i)  ) then
                delta_corner(i) = k ;
            end if
        end do  
    end do
    
    
    
    !! note that fixed_value is the same as normal n_x and n_y
    if( n_face .eq. 4 ) then
        is_edge_vertical = .false. 
        ny = -1.0_r_kind ;                                              !! the same as fixed-value: eta = -1.0
        dir = 1 ;                                                       !! parallel directon along boundary is x-axis
    elseif( ( n_face .eq. 1 ) .or. ( n_face .eq. 2 ) ) then
        is_edge_vertical = .true. 
        nx  = fixed_value   ;                                           !! noraml is either +1.0 or -1.0
        dir = 2 ;                                                       !! parallel directon along boundary is y-axis
    else 
        stop 'The boundary face for foundation is not valid in HW model for foundation'
    end if
    
    
    !!---------------------------------------------------------------!!
    !! set theta and gamma coefficients
    
    if( is_edge_vertical ) then
    
        theta_coeff(1) =   rho / ( alpha + 2.0_r_kind * mu ) ;
        theta_coeff(2) = - mu  / ( alpha + 2.0_r_kind * mu ) ;
        theta_coeff(3) = - nx *  ( alpha + mu ) / ( alpha + 2.0_r_kind * mu ) ;
        
        gamma_coeff(1) =    rho / mu ;
        gamma_coeff(2) =  - ( alpha + 2.0_r_kind * mu ) / mu ;
        gamma_coeff(3) =  - nx * ( alpha + mu ) /  mu  ;
        
    else 
        
        theta_coeff(1) =   rho / mu  ;
        theta_coeff(2) = - ( alpha + 2.0_r_kind * mu) / mu  ;
        theta_coeff(3) = - ny * ( alpha + mu ) /  mu  ;
        
        gamma_coeff(1) =    rho / ( alpha + 2.0_r_kind * mu ) ;
        gamma_coeff(2) =  - mu / ( alpha + 2.0_r_kind * mu ) ;
        gamma_coeff(3) =  - ny * ( alpha + mu ) / ( alpha + 2.0_r_kind * mu ) ;
   
    end if
     
    
     !!---------------------------------------------------------------!!
     !! 2.2.  calcuate local matrices for the three facial nodes
     
    Mmatrix(:,:)   = 0.0_r_kind  ;
    Kmatrix(:,:)   = 0.0_r_kind  ;
    Smatrix(:,:)   = 0.0_r_kind  ;
     
     
    n_gauss_point = element_property( elem_id , 3 ) 
    call  set_gauss_points( n_gauss_point )
    
    
     
     !! set normal and direction to be found
     
    do r = 1 ,  n_gauss_point 
              
        if( direction_value(1) .eq. 1 ) then                            !! find shape function on face
               
            call the_element%get_shape_function( elem_id , fixed_value , gauss_points(r)  )      ! xi = fixed
            call the_element%get_derivative_shape_function( elem_id , fixed_value , gauss_points(r) )
            d_volume = the_element%volume_element_at_face( elem_id , 1 , fixed_value , gauss_points(r) ) 
            d_volume = d_volume  * weight_of_gauss_point(r)  ; 
        else
               
            call the_element%get_shape_function( elem_id , gauss_points(r) , fixed_value  )      ! eta  = fixed
            call the_element%get_derivative_shape_function( elem_id , gauss_points(r) , fixed_value )
            d_volume = the_element%volume_element_at_face( elem_id , 2 , gauss_points(r) , fixed_value ) 
            d_volume =  d_volume * weight_of_gauss_point(r)  ;
        end if 
        
        
      
        do  k = 1 , n_bc_nodes   
            do n = 1 , n_bc_nodes 
              
               Mmatrix(k,n) = Mmatrix(k,n) +   &
            &  the_element%shape_vector(indexed_set(k)) * the_element%shape_vector(indexed_set(n)) * d_volume 
              
              
               Kmatrix(k , n) = Kmatrix(k,n) + &
            &  the_element%derivative_shape_matrix(dir,indexed_set(k)) *  &
            &  the_element%derivative_shape_matrix(dir,indexed_set(n)) * d_volume
            
            
               Smatrix(k , n) = Smatrix(k,n) + &
            &  the_element%shape_vector(indexed_set(k)) *  &
            &  the_element%derivative_shape_matrix(dir,indexed_set(n)) * d_volume
          
                          
            end do ! over n
        end do ! over k
                  
    end do ! over r
          
    
    !!--------------------------------------------------------------------------------------------------!!
    !!
    !! Boundary matrix that has one non-zero element: based on Samii BC
    !!
    !!---------------------------------------------------------------------------------------------------!!
    
    B14matrix(:,:) = 0.0_r_kind ; 
    B23matrix(:,:) = 0.0_r_kind ; 
    if( delta_corner(1) > 0 ) then
         n = delta_corner(1) ;
         B14matrix(n,n) = 1.0_r_kind ; 
    end if
    
    if( delta_corner(4) > 0 ) then
         n = delta_corner(4) ;
         B14matrix(n,n) = B14matrix(n,n) + 1.0_r_kind ; 
    end if
    
    if( delta_corner(2) > 0 ) then
         n = delta_corner(2) ;
         B23matrix(n,n) = 1.0_r_kind ; 
    end if
    
    if( delta_corner(3) > 0 ) then
         n = delta_corner(3) ;
         B23matrix(n,n) = B23matrix(n,n) + 1.0_r_kind ; 
    end if
    
    
     
!!---------------------------------------------------------------------!!
!!
!!  Step 3. Form row and col variable set index
!!
!!---------------------------------------------------------------------!!
    
     
    !! 3.1. get number of variables: N+M+2: get max num of variables for three boundary nodes
    !! To each displacement dof, there is N+M+2 variables, i.e. u(1), phi_i( N ) , phi_bar(1), phi_{N+i} (M) = N+M+2
     
    num_var = 1 + the_model_description%order_of_foundation_higher_order_bc(1) ;
    
    if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
        num_var = num_var + 1 + the_model_description%order_of_foundation_higher_order_bc(2) ; !! phi_bar and phi_{N+i}
    end if
               
    !! 3.2.1. There are num_dim of displacement dofs to each node, and also there are n_bc_nodes number of boundary nodes.
    num_row = n_bc_nodes * num_dim * num_var ;
    num_col = num_row ;
      
     
    
     
     ! ! 3.2. form index set for bigmatrices
     !! The order is 
     !!               ( U , V , Phi_1 , Psi_1 , Phi_2 , Psi_2 , ..., Phi_N , Psi_N , 
     !!                 PhiX  , PsiY , Phi_{N+1} , Psi_{N+1} , ... , Phi_{N+M} , Psi_{N+M} )
     !!
     !! where     
     !!      U = (u_{x,i1} , u_{x,i2} , u_{x,i3} )  : displacement along x-axis for three boundary nodes \{ i1 , i2 , i3 \}
     !!      V = (u_{y,i1} , u_{y,i2} , u_{y,i3} )  : displacement along y-axis for three boundary nodes \{ i1 , i2 , i3 \} 
     !!      Phi_k = ( phi_{k,i1} , phi_{k,i2} , phi_{k,i3} ) : the phi_vector along x-axis for boundary nodes
     !!      Psi_k = ( psi_{k,i1} , psi_{k,i2} , psi_{k,i3} ) : the psi_vector along y-axis for boundary nodes 
     !!      PhiX  = ( phi_{x,i1} , phi_{x,i2} , phi_{x,i3} ) : continuity function along x-axis for boundary nodes
     !!      PsiY  = ( psi_{y,i1} , psi_{y,i2} , psi_{y,i3} ) : continuity function along y-axis for boundary nodes
     !!
     
   
    row_variable_set(1:num_row) = 0 ;
             
    var_id = 0 ;
    do  i  = 1 , num_var
        do  r = 1 , num_dim    !! loop over dimension
            do  k = 1 , n_bc_nodes
                var_id  = var_id + 1 ;
                n = element_matrix( elem_id , indexed_set( k ) ) ;   !! each boundary node
             
                if( degree_of_freedom (n, 2 * r ) .eq. 0 ) cycle  ;  !! node may have no translational dof
             
                row_variable_set( var_id ) = degree_of_freedom (n, 2 * r ) + i -1 ;
            end do
        end do
    end do 
     
    
    col_variable_set(1:num_col ) = row_variable_set(1:num_row ) ;
    
   
    !! 3.3. prepare elemental matrices
    
    elemental_damping_matrix(:,:)    = 0.0_r_kind
    elemental_mass_matrix(:,:)       = 0.0_r_kind
    elemental_stiffness_matrix(:,:)  = 0.0_r_kind
    elemental_force_vector(:,:)      = 0.0_r_kind 
  
      
   
!!---------------------------------------------------------------------!! 
!!  
!! Step 4. Equation involving traction
!!   M \ddot{u} + K u + T = 0
!! where the traction is:
!!     T_x =  eta_0 M_3 \dot{U} - \eta_0 M_3 \dot{Phi}_1 + eta_2 S_3 V ;
!!     T_y =  eta_1 M_3 \dot{V} - eta_1 M_3  \dot{Psi}_1 + eta_3 S_3 V ;
!!     
!! The first two terms have already been inculded during bulk calculation.
!! The  traction term are now included in damping matrix and stiffness matrices.
!!
!! 1. The term phi_1 comes from first incident wave, so it enters when N + M > 0
!!
!!---------------------------------------------------------------------!!
    if( auto_coeff ) then
        do  i = 1 , the_model_description%order_of_foundation_higher_order_bc(1) + 1
            a(i) = theta_coeff(1) * ( c ** 2) + ellipitic_param ;
        end do
    end if 
    
    !! step 4.1. form coefficient.
    alpha_coeffs(:) = 0.0_r_kind
    if( is_edge_vertical ) then
        
        alpha_coeffs(1) = a(1) * ( alpha + 2.0_r_kind * mu ) / c ;      !! eta_0 = ( a0 / c ) * ( alpha + 2 mu )
        alpha_coeffs(2) = a(1) * mu  / c ;                              !! eta_1 = ( a0 / c ) * mu
        alpha_coeffs(3) = - alpha * nx ;                                !! eta_2 = -alpha * n_x
        alpha_coeffs(4) = - mu * nx ;                                   !! eta_3 = -mu * n_x
    else 
        
        alpha_coeffs(1) = a(1) * mu  / c ;                              !! eta_0 = ( a0 / c ) *  mu 
        alpha_coeffs(2) = a(1) * ( alpha + 2.0_r_kind * mu ) / c;       !! eta_1 = ( a0 / c ) * ( alpha + 2 mu ) 
        alpha_coeffs(3) = - mu * ny    ;                                !! eta_2 = -mu * n_y
        alpha_coeffs(4) = - alpha * ny ;                                !! eta_3 = -alpha * n_y
    end if
    
    !! 4.2. traction: stiffness: (row = U , col = V)  and (row = V , col = U)
    elemental_stiffness_matrix(1:3, 4:6) = alpha_coeffs(3) * Smatrix ;
    elemental_stiffness_matrix(4:6, 1:3) = alpha_coeffs(4) * Smatrix ;
    
    !! 4.3. traction:damping-1 (row = U , col = U)  and (row = V , col = V)
    elemental_damping_matrix(1:3, 1:3) = alpha_coeffs(1) * Mmatrix ;
    elemental_damping_matrix(4:6, 4:6) = alpha_coeffs(2) * Mmatrix ;
    
    !! 4.4. traction:damping-2 (row = U , col = phi1 )  and ( row = V , col = psi1 ) appears if N > 0
    if( the_model_description%order_of_foundation_higher_order_bc(1) + & 
      & the_model_description%order_of_foundation_higher_order_bc(2) > 0 )  then
      
        elemental_damping_matrix(1:3 , 7:9  ) = -alpha_coeffs(1) * Mmatrix ;
        elemental_damping_matrix(4:6 , 10:12) = -alpha_coeffs(2) * Mmatrix ;
    end if


    
    
    
!!!---------------------------------------------------------------------!! 
!!!  
!!! Step 5. Equation involving { phi_0, phi_1 , phi_2 }  and {psi_0 , psi_1 , psi_2 }
!!!   
!!!---------------------------------------------------------------------!!  

    if( the_model_description%order_of_foundation_higher_order_bc(1) > 0 ) then   !! if N > 0 
        
        !! 5.x. phi_1 equation or x-direction
        !! 5.x.0 set coefficients { d_i } used in this equation
       
        alpha_coeffs(:) = 0.0_r_kind ;   ! { d_i }
       
        alpha_coeffs(1) = 2.0_r_kind * a(2) * ( a(1) ** 2 - theta_coeff(1) * c ** 2 ) ; ! 2a_1 ( a_0^2 - theta_1 c^2 )
       
        !! d2 = -[ a0 * a1 ( a1 + 2 a0 ) + a0 * theta_1 * c^2 ]
        alpha_coeffs(2) = - a(1) * a(2) * ( a(2) + 2.0_r_kind * a(1) ) - a(1) * theta_coeff(1) * c ** 2 ;
        alpha_coeffs(3) = a(1) * ( a(2) ** 2 - theta_coeff(1) * c ** 2 ) ; !! d3 = a0 * ( a1^2 - theta_1 c^2 )
        alpha_coeffs(4) = - a(1) * a(2) * c * theta_coeff(3) ;             !! d4 = - a0 a1 c theta_3
        alpha_coeffs(5) = theta_coeff(2) * c ** 2 ;                        !! d5 = theta_2 * c^2 ;
        if( is_edge_vertical) then
            alpha_coeffs(6) = theta_coeff(2) * c * a(1) * a(2) * nx ;    !! 
        end if 
        alpha_coeffs(7) = theta_coeff(2) * c  ; 
        
        !! 5.x.1. add to mass matrix ( row = phi1 , col = phi0 , phi1 , phi2 )
        elemental_mass_matrix(7:9 , 1:3) = alpha_coeffs(1) * Mmatrix ;
        elemental_mass_matrix(7:9 , 7:9) = alpha_coeffs(2) * Mmatrix ;
        
        !! 5.x.1.1  row = phi_1 ,  col = phi_2 :  !! added if N+M+1 > 2 or N+M > 1 
        phi_id = the_model_description%order_of_foundation_higher_order_bc(1) + &
                 the_model_description%order_of_foundation_higher_order_bc(2) ;
                 
        if( phi_id  > 1 )  then
            elemental_mass_matrix( 7:9 , 13:15 ) = alpha_coeffs(3)  * Mmatrix ;
        end if
        
        
        !! 5.x.2. damping matrix row =phi1 , col = psi0 , psi1 , psi2
        elemental_damping_matrix(7:9 , 4:6    ) = -2.0_r_kind * ( alpha_coeffs(4) * Smatrix + alpha_coeffs(6) * B23matrix ) ;
        elemental_damping_matrix(7:9 , 10:12  ) =  alpha_coeffs(4) * Smatrix + alpha_coeffs(6) * B23matrix ;
        
        if( phi_id  > 1 )  then
           elemental_damping_matrix( 7:9 , 16:18 ) = alpha_coeffs(4)  * Smatrix + alpha_coeffs(6) * B23matrix ;
        end if
        
        !! 5.x.3. damping matrix row =phi1 , col = phi0 , phi1 , phi2
        elemental_damping_matrix(7:9 , 1:3    ) = 2.0_r_kind * a(2) * alpha_coeffs(7) * B14matrix  ;    !! phi_0
        elemental_damping_matrix(7:9 , 7:9  ) =  a(1) * alpha_coeffs(7) * B14matrix ;                   !! phi1
        
        if( phi_id  > 1 )  then
           elemental_damping_matrix( 7:9 , 13:15 ) = a(1) * alpha_coeffs(7) * B14matrix ;               !! phi2
        end if
        
        !! 5.x.3. stiffness matrix ( row = phi1 ,  col = phi0 , phi1 , phi2 ) plus BC of surface term of type = 1
        elemental_stiffness_matrix(7:9 , 1:3) = 2.0_r_kind * a(2) * alpha_coeffs(5) * Kmatrix ;
        elemental_stiffness_matrix(7:9 , 7:9) = a(1) * alpha_coeffs(5) * Kmatrix ;
                 
        if( phi_id  > 1 )  then
            elemental_stiffness_matrix( 7:9 , 13:15 ) = a(1) * alpha_coeffs(5) * Kmatrix ;
        end if
        
        
!!---------------------------------------------------------------------------------------------!!
        
       !! 5.y. psi_1 equation or y-direction
       !! 5.y.0 set coefficients used in this equation
        if( auto_coeff ) then
            do  i = 1 , the_model_description%order_of_foundation_higher_order_bc(1) + 1
                a(i) = gamma_coeff(1) * ( c ** 2) + ellipitic_param ;
            end do
        end if
    
       alpha_coeffs(:) = 0.0_r_kind ;   ! { e_i }
       
       alpha_coeffs(1) = 2.0_r_kind * a(2) * ( a(1) ** 2 - gamma_coeff(1) * c ** 2 ) ; ! e1 = 2a_1 ( a_0^2 - gamma_1 c^2 )
       
       !! e2 = -[ a0 * a1 ( a1 + 2 a0 ) + a0 * gamma_1 * c^2 ]
       alpha_coeffs(2) = - a(1) * a(2) * ( a(2) + 2.0_r_kind * a(1) ) - a(1) * gamma_coeff(1) * c ** 2 ;
       alpha_coeffs(3) = a(1) * ( a(2) ** 2 - gamma_coeff(1) * c ** 2 ) ; !! e3 = a0 * ( a1^2 - gamma_1 c^2 )
       alpha_coeffs(4) = - a(1) * a(2) * c * gamma_coeff(3) ;             !! e4 = - a0 a1 c gamma_3
       alpha_coeffs(5) = gamma_coeff(2) * c ** 2 ;                        !! e5 = gamma_2 * c^2 ;
       if( is_edge_vertical ) then
           alpha_coeffs(6) = gamma_coeff(2) * c * a(1) * a(2) * alpha * nx / ( alpha + 2.0_r_kind * mu ) ;   
       end if
       alpha_coeffs(7) = gamma_coeff(2) * c  ; 
       
       
       !! 5.y.1. add to mass matrix ( row = psi1 , col = psi0 , psi1 , psi2 )
        elemental_mass_matrix(10:12 , 4:6  ) = alpha_coeffs(1) * Mmatrix ;
        elemental_mass_matrix(10:12 , 10:12) = alpha_coeffs(2) * Mmatrix ;
        
        !! 5.y.1.1  row = psi_1 ,  col = psi_2 :  !! added if N+M+1 > 2 or N+M > 1 
        if( phi_id  > 1 )  then
           elemental_mass_matrix( 10:12 , 16:18 ) = alpha_coeffs(3)  * Mmatrix ;
        end if
        
        
        !! 5.y.2. damping matrix row =psi1 , col = phi0 , phi1 , phi2
        elemental_damping_matrix(10:12 , 1:3    ) = -2.0_r_kind * ( alpha_coeffs(4) * Smatrix + alpha_coeffs(6) * B23matrix ) ;
        elemental_damping_matrix(10:12 , 7:9    ) =  alpha_coeffs(4) * Smatrix + alpha_coeffs(6) * B23matrix ;
        
        if( phi_id  > 1 )  then
           elemental_damping_matrix( 10:12 , 13:15 ) = alpha_coeffs(4)  * Smatrix + alpha_coeffs(6) * B23matrix ;
        end if
        
        !! 5.y.3. damping matrix row =psi1 , col = psi0 , psi1 , psi2
        elemental_damping_matrix(10:12 , 4:6    ) = 2.0_r_kind * a(2) * alpha_coeffs(7) * B14matrix  ;
        elemental_damping_matrix(10:12 , 10:12  ) =  a(1) * alpha_coeffs(7) * B14matrix  ;
        
        if( phi_id  > 1 )  then
           elemental_damping_matrix( 10:12 , 16:18 ) = a(1) * alpha_coeffs(7) * B14matrix  ;
        end if
        
        !! 5.y.4. stiffness matrix ( row = psi1 ,  col = psi0 , psi1 , psi2 )
        elemental_stiffness_matrix(10:12 , 4:6  ) = 2.0_r_kind * a(2) * alpha_coeffs(5) *  Kmatrix ;
        elemental_stiffness_matrix(10:12 , 10:12) = a(1) * alpha_coeffs(5) *  Kmatrix  ;
                 
        if( phi_id  > 1 )  then
           elemental_stiffness_matrix( 10:12 , 16:18 ) = a(1) * alpha_coeffs(5) * Kmatrix ;
        end if
        
        
        

    end if ! of phi1-psi1 eqn
    
    

     
!!!---------------------------------------------------------------------!!
!!!
!!!  Step 6.x. Add equations governing auxilary functions 
!!!                     { \phi_{j-1} , \phi_j , \phi_{j+1} }  
!!!
!!!     
!!!   1. The first equation is added to phi_j  ( 2 <= j <= N )
!!!   2. The term phi_{N+1} does not appear if M = 0,
!!!          
!!!---------------------------------------------------------------------!!  

    if( auto_coeff ) then
        do  i = 1 , the_model_description%order_of_foundation_higher_order_bc(1) + 1
            a(i) = theta_coeff(1) * ( c ** 2) + ellipitic_param ;
        end do
    end if
    
    do phi_id = 2 , the_model_description%order_of_foundation_higher_order_bc(1)
        
        
       !! 6.x. phi_j equation or x-direction
       !! 6.x.0 set coefficients { d_i } used in this equation
       
       alpha_coeffs(:) = 0.0_r_kind ;   
       
       !! d1 = a_j ( a_{j-1}^2 - theta_1 c^2 )
       alpha_coeffs(1) =  a( phi_id+1 ) * ( a(phi_id) ** 2 - theta_coeff(1) * c ** 2 ) ; 
       
       !! d2 = -( a_j + a_{j-1} ) [ a_j * a_{j-1}  +  theta_1 * c^2 ]
       alpha_coeffs(2) = - ( a( phi_id+1 ) + a(phi_id) ) * ( a( phi_id+1 ) * a(phi_id) + theta_coeff(1) * c ** 2 ) ;
       
       !! d3 = a_{j-1} ( a_j^2 - theta_1 c^2 )
       alpha_coeffs(3) =  a( phi_id ) * ( a(phi_id + 1) ** 2 - theta_coeff(1) * c ** 2 ) ; 
       alpha_coeffs(4) = -a(phi_id + 1) * a(phi_id) * c * theta_coeff(3) ;          !! d4 = - a_j a_{j-1} c theta_3
       alpha_coeffs(5) = theta_coeff(2) * c ** 2 ;                                  !! d5 = theta_2 * c^2 ;
       if( is_edge_vertical ) then
           alpha_coeffs(6) = theta_coeff(2) * c * a(phi_id + 1) * a(phi_id) * nx ;  
       end if
       alpha_coeffs(7) = theta_coeff(2) * c ;
       
       !! 6.x. row = phi_j : adding to equation phi_j for  2 <= j <= N
       row = phi_id * 2 * n_bc_nodes ;
       col = row - 2 * n_bc_nodes ;
       
       !! 6.x.1. add to mass matrix ( row = phi_j , col = phi_{j-1} , phi_j , phi_{j+1} )
        elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(1) * Mmatrix ;
        elemental_mass_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(2) * Mmatrix ;
                 
        !! 6.x.1.2. if j < N then phi_{j+1} is at most phi_N, so add the term to incident-wave part         
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 2 * n_bc_nodes ;
            elemental_mass_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(3)  * Mmatrix ;
        else
            !! Then phi_{j+1} = phi_{N+1} that appears only if M > 0
            !! In this case phi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 4 * n_bc_nodes ;
                elemental_mass_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(3)  * Mmatrix ;
            end if
        end if
        
        
        !! 6.x.2. damping matrix row =phi_j , col = psi_{j-1} , psi_j , psi_{j+1}
        col = row - n_bc_nodes ;
        elemental_damping_matrix(row+1:row+3,col+1:col+3) = -alpha_coeffs(4) *Smatrix -alpha_coeffs(6)*B23matrix ;
        
        !! 6.x.2.1. if j < N then psi_{j+1} is at most psi_N, so add the term to incident-wave part         
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 3 * n_bc_nodes ;
       elemental_damping_matrix(row+1:row+3,col+1:col+3) = alpha_coeffs(4) * Smatrix + alpha_coeffs(6)* B23matrix ;
        else
            !! Then psi_{j+1} = psi_{N+1} that appears only if M > 0
            !! In this case psi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 5 * n_bc_nodes ;
           elemental_damping_matrix(row+1:row+3,col+1:col+3) = alpha_coeffs(4)*Smatrix + alpha_coeffs(6)*B23matrix ; 
            end if
        end if
        
        !! 6.x.3. damping matrix row =phi_j , col = phi_{j-1} , phi_j , phi_{j+1}
        col = row - 2 * n_bc_nodes ;
        elemental_damping_matrix(row+1:row+3,col+1:col+3) = a(phi_id+1) * alpha_coeffs(7) *B14matrix ;
        elemental_damping_matrix(row+1:row+3,row+1:row+3) = ( a(phi_id) + a(phi_id+1) ) * alpha_coeffs(7) *B14matrix ;
        
        !! 6.x.3.1. if j < N then phi_{j+1} is at most phi_N, so add the term to incident-wave part         
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 2 * n_bc_nodes ;
       elemental_damping_matrix(row+1:row+3,col+1:col+3) = a(phi_id) * alpha_coeffs(7) *B14matrix ;
        else
            !! Then phi_{j+1} = phi_{N+1} that appears only if M > 0
            !! In this case phi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 4 * n_bc_nodes ;
           elemental_damping_matrix(row+1:row+3,col+1:col+3) = a(phi_id) * alpha_coeffs(7) *B14matrix ; 
            end if
        end if
        
        
        !! 6.x.4. add to stiffness matrix ( row = phi_j , col = phi_{j-1} , phi_j , phi_{j+1} )
        col = row - 2 * n_bc_nodes ;
       
        elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(5) * a(phi_id+1) * Kmatrix  ;
        elemental_stiffness_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(5) *(a(phi_id+1) + a(phi_id)) * Kmatrix ;
                 
        !! 6.x.3.2. if j < N then phi_{j+1} is at most phi_N, so add the term to incident-wave part  
               
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 2 * n_bc_nodes ;
        elemental_stiffness_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(5) * a(phi_id) * Kmatrix ;
        else
            !! Then phi_{j+1} = phi_{N+1} that appears only if M > 0
            !! In this case phi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 4 * n_bc_nodes ;
        elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(5) * a(phi_id) * Kmatrix ;
            end if
        end if
        

    end do !! over phi_id 
    
    
    !!!---------------------------------------------------------------------!!
!!!
!!!  Step 6.y. Add equations governing auxilary functions 
!!!                     { \psi_{j-1} , \psi_j , \psi_{j+1} } 
!!!
!!!   1. The second equation is added to psi_j ( 2 <= j <= N )
!!!   2. The term phi_{N+1} does not appear if M = 0,
!!!          
!!! 
!!!---------------------------------------------------------------------!!  
    
    if( auto_coeff ) then
        do  i = 1 , the_model_description%order_of_foundation_higher_order_bc(1) + 1
            a(i) = gamma_coeff(1) * ( c ** 2) + ellipitic_param ;
        end do
    end if
    
    do phi_id = 2 , the_model_description%order_of_foundation_higher_order_bc(1)
        
        
       !! 6.y. psi_j equation or y-direction
       !! 6.y.0 set coefficients { d_i } used in this equation
       
       alpha_coeffs(:) = 0.0_r_kind ;   
       
       ! d1 = a_j ( a_{j-1}^2 - gamma_1 c^2 )
       alpha_coeffs(1) =  a( phi_id+1 ) * ( a(phi_id) ** 2 - gamma_coeff(1) * c ** 2 ) ; 
       
       !! d2 = -( a_j + a_{j-1} ) [ a_j * a_{j-1}  +  gamma_1 * c^2 ]
       alpha_coeffs(2) = -( a( phi_id+1 ) + a(phi_id) ) * ( a( phi_id+1 ) * a(phi_id) + gamma_coeff(1) * c ** 2 )
       
       ! d3 = a_{j-1} * ( a_j^2 - gamma_1 c^2 )
       alpha_coeffs(3) =  a( phi_id ) * ( a(phi_id + 1) ** 2 - gamma_coeff(1) * c ** 2 ) ; 
       alpha_coeffs(4) = -a(phi_id + 1) * a(phi_id) * c * gamma_coeff(3) ;         !! d4 = - a_j a_{j-1} c gamma_3
       alpha_coeffs(5) =  gamma_coeff(2) * c ** 2 ;                                !! d5 = gamma_2 * c^2 ;
       if( is_edge_vertical ) then
           alpha_coeffs(6) = gamma_coeff(2) * a(phi_id + 1) * a(phi_id) *c* nx * alpha /( alpha + 2.0_r_kind * mu) ;
       end if 
       alpha_coeffs(7) = gamma_coeff(2) * c ;
       
       !! 6.y. row = psi_j : adding to equation psi_j for  2 <= j <= N, row = psi_{j-1} , psi_j  , psi_{j+1}
       row = ( 2 * phi_id + 1 ) * n_bc_nodes ;
       
       !! 6.y.1:  row = psi_j ,  col = psi_{j-1} 
       col = row - 2 * n_bc_nodes ;
       
       !! 6.y.1. add to mass matrix ( row = psi_j , col = psi_{j-1} , psi_j , psi_{j+1} )
        elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(1) * Mmatrix ;
        elemental_mass_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(2) * Mmatrix ;
                 
        !! 6.y.1.2. if j < N then psi_{j+1} is at most psi_N, so add the term to incident-wave part         
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 2 * n_bc_nodes ;
            elemental_mass_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(3)  * Mmatrix ;
        else
            !! Then psi_{j+1} = psi_{N+1} that appears only if M > 0
            !! In this case psi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 4 * n_bc_nodes ;
                elemental_mass_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(3)  * Mmatrix ;
            end if
        end if
        
        
        !! 6.y.2. damping matrix row =psi_j , col = phi_{j-1} , phi_{j+1}
        col = row - 3 * n_bc_nodes ;
        elemental_damping_matrix(row+1:row+3 , col+1:col+3) = - alpha_coeffs(4) * Smatrix   & 
                                                            & - alpha_coeffs(6) * B23matrix ;
        
        !! 6.y.2.1. if j < N then phi_{j+1} is at most psi_N, so add the term to incident-wave part         
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + n_bc_nodes ;
            elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(4)  * Smatrix   &
                                                                & + alpha_coeffs(6)  * B23matrix ;
        else
            !! Then psi_{j+1} = psi_{N+1} that appears only if M > 0
            !! In this case psi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 3 * n_bc_nodes ;
                elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(4) * Smatrix   &
                                                                    & + alpha_coeffs(6) * B23matrix ;
            end if
        end if
        
        !! 6.y.3. damping matrix row =psi_j , col = psi_{j-1} , psi_j , psi_{j+1}
        col = row - 2 * n_bc_nodes ;
        elemental_damping_matrix(row+1:row+3 , col+1:col+3) = a(phi_id+1) * alpha_coeffs(7) * B14matrix ;
        elemental_damping_matrix(row+1:row+3 , row+1:row+3) = (a(phi_id) + a(phi_id+1)) * alpha_coeffs(7) * B14matrix ;
        
        !! 6.y.3.1. if j < N then psi_{j+1} is at most psi_N, so add the term to incident-wave part         
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 2 * n_bc_nodes ;
            elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = a(phi_id) * alpha_coeffs(7) * B14matrix ;
        else
            !! Then psi_{j+1} = psi_{N+1} that appears only if M > 0
            !! In this case psi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 4 * n_bc_nodes ;
                elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = a(phi_id) * alpha_coeffs(7) * B14matrix ;
            end if
        end if
        
        
        !! 6.y.4. add to stiffness matrix ( row = psi_j , col = psi_{j-1} , psi_j , psi_{j+1} )
        col = row - 2 * n_bc_nodes ;
       
        elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(5) * a(phi_id+1) * Kmatrix ;
        elemental_stiffness_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(5) *(a(phi_id+1) + a(phi_id)) * Kmatrix ;
                 
        !! 6.y.4.2. if j < N then psi_{j+1} is at most psi_N, so add the term to incident-wave part  
               
        if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(1) )  then
            col = row + 2 * n_bc_nodes ;
        elemental_stiffness_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(5) * a(phi_id) * Kmatrix  ;
        else
            !! Then psi_{j+1} = psi_{N+1} that appears only if M > 0
            !! In this case psi_{N+1} appears after phi_bar and psi_bar
            
            if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
                col = row + 4 * n_bc_nodes ;
            elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(5) * a(phi_id) *Kmatrix  ;
            end if
        end if
        
        

    end do !! over phi_id


!!---------------------------------------------------------------------!!
!!
!!  Step 7.x. Add equations containing  { D_t psi_bar ,    D_t psi_bar }
!!          All the remaining equations appears only if M > 0
!!  
!!          These equations are added to phi_bar and psi_bar equations
!!
!!---------------------------------------------------------------------!! 

    !! 7.0. check if M > 0 and N >= 1 as it used one of the equations phi_j where j \ne 0
    
    if ( ( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) .and. &
       & ( the_model_description%order_of_foundation_higher_order_bc(1) > 0 ) ) then
   
   
       phi_id = the_model_description%order_of_foundation_higher_order_bc(1)      ! take as N 
       
       !! 7.x.0. set coefficients used in this equation
       alpha_coeffs(:) = 0.0_r_kind ;   
       
       ! d1 =  a_N^2 - theta_1 c^2 
       alpha_coeffs(1) =  a( phi_id + 1 ) ** 2 - theta_coeff(1) * c ** 2  ; 
       
       !! d2 = -[ a_N^2 +  theta_1 * c^2 ]
       alpha_coeffs(2) = - ( a( phi_id+1 ) ** 2  + theta_coeff(1) * c ** 2 ) ;
       
       !! d3 = 2 a_N c n_n 
       if( is_edge_vertical ) then
           alpha_coeffs(3) = 2.0_r_kind * a(phi_id + 1) * c  * nx ;                 !! d3 = 2 a_N c nx
       else 
           alpha_coeffs(3) = 2.0_r_kind * a(phi_id + 1) * c  * ny ;                 !! d3 = 2 a_N c ny
       end if 
       
       alpha_coeffs(4) = theta_coeff(2) * c ** 2 ;                                  !! d4 = theta_2 * c^2 ;
       alpha_coeffs(5) = a(phi_id + 1)  * c * theta_coeff(3) ;                      !! d5 =  a_N c theta_3
       if( is_edge_vertical ) then
           alpha_coeffs(6) = theta_coeff(2) * c * a(phi_id + 1) * nx  ;             !! surface_bc_type = 2
       end if 
       alpha_coeffs(7) = theta_coeff(2) * c ;
       
       !! 7.x. row = phi_bar : ( N (phi) + 1 ( phi_0) ) * 2 * n_bc_nodes
       row = ( phi_id + 1 ) * 2 * n_bc_nodes ;
       
       !! 7.x.1:  mass matrix (row = phi_bar ,  col = ( phi_N , phi_{N+1} , 
       col = row - 2 * n_bc_nodes ;
       
       elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(1) * Mmatrix ;
       col = row + 2 * n_bc_nodes ;
       elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(2) * Mmatrix ;
                 
       
       !! 7.x.2. damping matrix
       !! 7.x.2.1. row = phi_bar , col = psi_N , phi_bar , psi_{N+1}
       col = row - n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =  alpha_coeffs(5) * Smatrix - alpha_coeffs(6) * B23matrix ;
       col = row + 3 * n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = -alpha_coeffs(5) * Smatrix + alpha_coeffs(6) * B23matrix ;
       elemental_damping_matrix( row+1:row+3 , row+1:row+3 ) =  alpha_coeffs(3)  * Mmatrix ;
       
       !! 7.x.3. damping matrix
       !! 7.x.3.1. row = phi_bar , col = phi_N , phi_{N+1}
       col = row - 2 * n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =  alpha_coeffs(7) * B14matrix ;
       col = row + 2 * n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =  alpha_coeffs(7) * B14matrix ;
       
       
       !! 7.x.4. stiffness matrix: row = phi_bar , col = phi_N , phi_{N+1}
       col = row - 2 * n_bc_nodes ;
       
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(4) * Kmatrix  ;
       col = row + 2 * n_bc_nodes ;
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(4) * Kmatrix  ;
       
       
       
!!-------------------------------------------------------------------------------------------!! 
       
       !! 7.y.0. set coefficients used in this equation
       alpha_coeffs(:) = 0.0_r_kind ;   
       
       ! d1 =  a_N^2 - gamma_1 c^2 
       alpha_coeffs(1) =  a( phi_id + 1 ) ** 2 - gamma_coeff(1) * c ** 2  ; 
       
       !! d2 = -[ a_N^2 +  gamma_1 * c^2 ]
       alpha_coeffs(2) = - ( a( phi_id+1 ) ** 2  + gamma_coeff(1) * c ** 2 ) ;
       
       !! d3 = 2 a_N c n_n 
       if( is_edge_vertical ) then
           alpha_coeffs(3) = 2.0_r_kind * a(phi_id + 1) * c  * nx ;                 !! d3 = 2 a_N c nx
       else 
           alpha_coeffs(3) = 2.0_r_kind * a(phi_id + 1) * c  * ny ;                 !! d3 = 2 a_N c ny
       end if 
       
       alpha_coeffs(4) = gamma_coeff(2) * c ** 2 ;                                  !! d4 = gamma_2 * c^2 ;
       alpha_coeffs(5) = a(phi_id + 1)  * c * gamma_coeff(3) ;                      !! d5 =  a_N c gamma_3
       if( is_edge_vertical ) then
           alpha_coeffs(6) = gamma_coeff(2) * c * a(phi_id + 1) * nx * alpha / ( alpha + 2.0_r_kind * mu ) ;   !! surface_bc_type = 2
       end if 
       alpha_coeffs(7) = gamma_coeff(2) * c ;
       
       !! 7.y. row = psi_bar : id( phi_bar ) + n_bc_nodes
       row = ( phi_id + 1 ) * 2 * n_bc_nodes + n_bc_nodes ;
       
       !! 7.y.1:  row = psi_bar ,  col = psi_N 
       col = row - 2 * n_bc_nodes ;
       
       !! 7.y.1:  mass matrix (row = psi_bar ,  col = ( psi_N , psi_{N+1} , 
       col = row - 2 * n_bc_nodes ;
       
       elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(1) * Mmatrix ;
       col = row + 2 * n_bc_nodes ;
       elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(2) * Mmatrix ;
                 
       
       !! 7.y.2. damping matrix
       !! 7.y.2.1. row = psi_bar , col = phi_N , psi_bar , phi_{N+1}
       col = row - 3 * n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =  alpha_coeffs(5)  * Smatrix - alpha_coeffs(6) * B23matrix ;
       col = row + n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = -alpha_coeffs(5)  * Smatrix + alpha_coeffs(6) * B23matrix ;
       elemental_damping_matrix( row+1:row+3 , row+1:row+3 ) =  alpha_coeffs(3)  * Mmatrix ;
       
       !! 7.y.3. damping matrix
       !! 7.y.3.1. row = psi_bar , col = psi_N , psi_{N+1}
       col = row - 2 * n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =  alpha_coeffs(7) * B14matrix ;
       col = row + 2 * n_bc_nodes ;
       elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(7) * B14matrix ;
      
       
       
       !! 7.y.4. stiffness matrix: row = psi_bar , col = psi_N , psi_{N+1}
       col = row - 2 * n_bc_nodes ;
       
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(4) * Kmatrix  ;
       col = row + 2 * n_bc_nodes ;
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(4) * Kmatrix  ;
       
        
       
    end if !! of first continuity eqn

    
!!---------------------------------------------------------------------!!
!!
!!  Step 8. Add equations containing  psi_bar and phi_bar
!!          This equation appears only if M > 0, as it uses first equation of \phi_{N+j},
!!          The term phi_{N+2} appears only if M > 1, i.e. if N+2 < N + M+ 1 --> M > 1
!!---------------------------------------------------------------------!!
   
   if( the_model_description%order_of_foundation_higher_order_bc(2) > 0 ) then
     
       phi_id = the_model_description%order_of_foundation_higher_order_bc(1)      ! take as N 
       
       !! 8.x.0. set coefficients used in this equation
       alpha_coeffs(:) = 0.0_r_kind ;   
       
       alpha_coeffs(1) =  theta_coeff(1) * c ** 2  ;                    !! d1 =  theta_1 c^2
       alpha_coeffs(2) = b(1) ** 2  ;                                   !! d2 = b_1^2

       !! d3 = 2 b(1) c n_n 
       if( is_edge_vertical ) then
           alpha_coeffs(3) = 2.0_r_kind * b(1) * c  * nx ;              !! d3 = 2 b_1 c nx
       else 
           alpha_coeffs(3) = 2.0_r_kind * b(1) * c  * ny ;              !! d3 = 2 b_1 c ny
       end if 
       
       alpha_coeffs(4) = -theta_coeff(2) * c ** 2 ;                     !! d4 = -theta_2 * c^2 ;
       alpha_coeffs(5) = -b(1)  * c * theta_coeff(3) ;                  !! d5 = -theta_3 * b_1 * c
       if( is_edge_vertical ) then
           alpha_coeffs(6) = -theta_coeff(2) * c * b(1) * nx ;          !! surface_bc_type = 2
       end if
       alpha_coeffs(7) = -theta_coeff(2) * c ;


       !! 8.x. row = phi_{N+1} : ( N (phi) + 1 ( phi_0) + 1( phi_bar) ) * 2 * n_bc_nodes
       row = ( phi_id + 2 ) * 2 * n_bc_nodes ;
       
       !! 8.x.1: mass matrix:  row = phi_{N+1} ,  col = phi_{N+1} , phi_{N+2} 
       elemental_mass_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(1) * Mmatrix ;
       
       if( the_model_description%order_of_foundation_higher_order_bc(2) > 1 ) then
           col = row + 2 * n_bc_nodes ;
           elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(1) * Mmatrix ;
       end if          
       
       !! 8.x.2. damping matrix: 
       elemental_damping_matrix( row+1:row+3 , row+1:row+3 ) =  alpha_coeffs(7) * B14matrix ;
       
       if( the_model_description%order_of_foundation_higher_order_bc(2) > 1 ) then
           col = row + 2 * n_bc_nodes ;
           elemental_damping_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(7) * B14matrix ;
       end if
       
       
       
       !! 8.x.3. stiffness matrix: row = phi__{N+1} , col = phi_bar , phi_{N+1} , psi_{N+1} , phi_{N+2} , psi_{N+2}
       col = row - 2 * n_bc_nodes ; !! for phi_bar
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(3) * Mmatrix ;  !! phi_bar
       
    elemental_stiffness_matrix(row+1:row+3 , row+1:row+3)= alpha_coeffs(2) *Mmatrix +alpha_coeffs(4) * Kmatrix ;  !! phi_{N+1}
 
     
       col = row + n_bc_nodes ;
  elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(5) * Smatrix - alpha_coeffs(6)* B23matrix ;  !! psi_{N+1}
       
       if( the_model_description%order_of_foundation_higher_order_bc(2) > 1 ) then
       
           col = row + 2 * n_bc_nodes ;
           
           elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = - alpha_coeffs(2) * Mmatrix    &   !! phi_{N+2}
                                                                 & + alpha_coeffs(4) * Kmatrix    ;   !! phi_{N+2}
           
           col = row + 3 * n_bc_nodes ;
  elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = -alpha_coeffs(5) * Smatrix + alpha_coeffs(6) * B23matrix ;  !! psi_{N+2}
       end if 
       
       
       
!!-------------------------------------------------------------------------------------------!! 


       !! 8.y.0. set coefficients used in this equation
       alpha_coeffs(:) = 0.0_r_kind ;   
       
       alpha_coeffs(1) =  gamma_coeff(1) * c ** 2  ;                    !! d1 =  gamma_1 c^2
       alpha_coeffs(2) = b(1) ** 2  ;                                   !! d2 = b_1^2

       !! d3 = 2 b(1) c n_n 
       if( is_edge_vertical ) then
           alpha_coeffs(3) = 2.0_r_kind * b(1) * c  * nx ;              !! d3 = 2 b_1 c nx
       else 
           alpha_coeffs(3) = 2.0_r_kind * b(1) * c  * ny ;              !! d3 = 2 b_1 c ny
       end if 
       
       alpha_coeffs(4) = -gamma_coeff(2) * c ** 2 ;                     !! d4 = -gamma_2 * c^2 ;
       alpha_coeffs(5) = -b(1)  * c * gamma_coeff(3) ;                  !! d5 = -gamma_3 * b_1 * c
       if( is_edge_vertical ) then
           alpha_coeffs(6) = -gamma_coeff(2) * c * b(1) * nx * alpha / ( alpha + 2.0_r_kind * mu ) ;    !! surface_bc_type = 2
       end if
       alpha_coeffs(7) = -gamma_coeff(2) * c ;

       !! 8.y. row = psi_{N+1} : ( N (phi) + 1 ( phi_0) + 1( phi_bar) ) * 2 * n_bc_nodes + n_bc_nodes
       row = ( phi_id + 2 ) * 2 * n_bc_nodes + n_bc_nodes ;
       
       !! 8.y.1: mass matrix:  row = psi_{N+1} ,  col = psi_{N+1} , psi_{N+2} 
       elemental_mass_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(1) * Mmatrix ;
       
       if( the_model_description%order_of_foundation_higher_order_bc(2) > 1 ) then
           col = row + 2 * n_bc_nodes ;
           elemental_mass_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(1) * Mmatrix ;
       end if          
       
       !! 8.y.2. damping matrix:
       elemental_damping_matrix( row+1:row+3 , row+1:row+3 ) =  alpha_coeffs(7) * B14matrix ;
       
       if( the_model_description%order_of_foundation_higher_order_bc(2) > 1 ) then
           col = row + 2 * n_bc_nodes ;
           elemental_damping_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(7) * B14matrix ;
       end if
       
       !! 8.y.3. stiffness matrix: row = psi_{N+1} , col = psi_bar , phi_{N+1} , psi_{N+1} , phi_{N+2} , psi_{N+2}
       col = row - 2 * n_bc_nodes ; !! for psi_bar
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(3) * Mmatrix ;  !! psi_bar
       
       col = row -  n_bc_nodes ; !! for phi_{N+1}
       elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = alpha_coeffs(5) * Smatrix  &
                                            & - alpha_coeffs(6) * B23matrix ;  !! phi_{N+1}
       
       elemental_stiffness_matrix(row+1:row+3 , row+1:row+3) = alpha_coeffs(2) * Mmatrix  &              !! psi_{N+1}
                                                           & + alpha_coeffs(4) * Kmatrix  ;  !! psi_{N+1}
     
      
       if( the_model_description%order_of_foundation_higher_order_bc(2) > 1 ) then
           col = row + n_bc_nodes ; !! for phi_{N+2}
   elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = -alpha_coeffs(5) * Smatrix  &
                                        & + alpha_coeffs(6) * B23matrix ;  !! phi_{N+2}
           
           col = row + 2 * n_bc_nodes ; !! for psi_{N+2}
   elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = -alpha_coeffs(2) * Mmatrix  &
                                        & + alpha_coeffs(4) *  Kmatrix ; !! psi_{N+2}
           
       end if 
       

    end if !! of second continuity equation
    
    
     
!!---------------------------------------------------------------------!!
!!
!!  Step 9.x Add equations governing: phi_{N+j , N+j + 1 , N+j+2} : for j = 1 , M-1
!!  
!!  Also N+j+2 <=  N+M+1 --> j <= M-1, and j>= 1 , so  1 <= j <= M-1
!!  For the last term, \phi_{N+j+2} does not appear.
!!
!!---------------------------------------------------------------------!!  

    do phi_id = 1 , the_model_description%order_of_foundation_higher_order_bc(2)-1
       
      
      !! 9.x. phi_{N+j+1} equation or x-direction
      !! 9.x.0 set coefficients { d_i } used in this equation
       
      alpha_coeffs(:) = 0.0_r_kind ;   
       
      
      alpha_coeffs(1) =   theta_coeff(1) * c ** 2  ;                    !! d1 =   theta_1 c^2 
      alpha_coeffs(2) = - b( phi_id ) * b( phi_id + 1 )  ;              !! d2 = - b_j * b_{j+1} 
      alpha_coeffs(3) = - theta_coeff(2) * c ** 2  ;                    !! d3 = - theta_2 c^2
      alpha_coeffs(4) =   b( phi_id ) * b( phi_id + 1 ) * c *  theta_coeff(3) ;  !! d4 = b_j * b_{j+1} c theta_3
      if( is_edge_vertical ) then
          alpha_coeffs(5) = - theta_coeff(2) * c * b( phi_id ) * b( phi_id + 1 ) * nx  ;  !! surface_bc_type = 2
      end if 
      alpha_coeffs(6) = - theta_coeff(2) * c 
       
      !! 9.x. row = phi_{N+j+1} : for  1 <= j <= M-1: 
      !! index: 1( phi_0 ) + N( phi_j ) + 1( phi_bar ) + 1( phi_{N+1} ) + j:  begins from ( N + 2 + j ) * 2 * n_bc_nodes
      
      row = ( the_model_description%order_of_foundation_higher_order_bc(1) + 2 + phi_id) * 2 * n_bc_nodes ;
       
      !! 9.x.1. add to mass matrix ( row = phi_{N+j+1} , col = phi_{N+j} , phi_{N+j+1} , phi_{N+j+2} )
      col = row - 2 * n_bc_nodes ; !! for phi_{N+j}
      elemental_mass_matrix(row+1:row+3 , col+1:col+3) = b( phi_id + 1 ) * alpha_coeffs(1) * Mmatrix ;
      elemental_mass_matrix(row+1:row+3 , row+1:row+3) = ( b(phi_id) + b( phi_id + 1 ) ) * & 
                                                       &   alpha_coeffs(1) * Mmatrix ;
                 
      !! 9.x.1.2. if N + j + 2 < N + M+1 or j < M-1 then phi_{N+j+2} appears.
      if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(2) -1 )  then
          col = row + 2 * n_bc_nodes ; 
          elemental_mass_matrix( row+1:row+3 , col+1:col+3 ) = b(phi_id) * alpha_coeffs(1)  * Mmatrix ; !!  phi_{N+j+2}
      end if
        
        
      !! 9.x.2. damping matrix:
      col = row - 2 * n_bc_nodes ; !! for phi_{N+j}
      elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =   b( phi_id + 1 ) * alpha_coeffs(6) * B14matrix ;
      
      elemental_damping_matrix( row+1:row+3 , row+1:row+3 ) = (b(phi_id) + b(phi_id+1)) * alpha_coeffs(6) * B14matrix ;
       
      if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(2) -1 )  then
          col = row + 2 * n_bc_nodes ; 
          elemental_damping_matrix(row+1:row+3 , col+1:col+3) =  b(phi_id) * alpha_coeffs(6) * B14matrix ;
      end if
       
      !! 9.x.3. add to stiffness matrix ( row = phi_{N+j+1} , col = phi_{N+j} , psi_{N+j} , phi_{N+j+1} , phi_{N+j+2} , psi )
      col = row - 2 * n_bc_nodes ;  !! for phi_{N+j}
       
      elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) =   &           
   &  b(phi_id) * alpha_coeffs(2) * Mmatrix + b(phi_id+1) * alpha_coeffs(3) * Kmatrix  ;     !! phi_{N+j}
     
   
      
      col = row - n_bc_nodes ;  !! for psi_{N+j}
      elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = - alpha_coeffs(4) * Smatrix &
                                                            & - alpha_coeffs(5) * B23matrix ;   !! psi_{N+j}
      
      elemental_stiffness_matrix(row+1:row+3 , row+1:row+3) = ( b(phi_id+1) + b(phi_id) ) *  & !! phi_{N+j+1}
    & ( -alpha_coeffs(2) * Mmatrix + alpha_coeffs(3) * Kmatrix  ) ; 

      !! row = phi_{N+j+1} , col = phi_{N+j+2} , psi_{N+j+2}
      if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(2) -1 )  then
      
          col = row + 2 * n_bc_nodes ;   !! for phi_{N+j+2}
          elemental_stiffness_matrix( row+1:row+3 , col+1:col+3 ) =  &
      &   ( b(phi_id+1) * alpha_coeffs(2)  * Mmatrix + b(phi_id) * alpha_coeffs(3)  * Kmatrix ) ; !!  phi_{N+j+2}
      
          col = row + 3 * n_bc_nodes ;   !! for psi_{N+j+2}
  elemental_stiffness_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(4)  * Smatrix &
                                                        & + alpha_coeffs(5) * B23matrix ; 
      end if
      
                
    
   end do !! of evanescence waves in x-dir
    
   
    
!!---------------------------------------------------------------------!!
!!
!!  Step 9.y Add equations governing: psi_{N+j , N+j + 1 , N+j+2} : for j = 1 , M-1
!!  
!!  Also N+j+2 <=  N+M+1 --> j <= M-1, and j>= 1 , so  1 <= j <= M-1
!!  For the last term, \phi_{N+j+2} does not appear.
!!
!!---------------------------------------------------------------------!!  

   do phi_id = 1 , the_model_description%order_of_foundation_higher_order_bc(2)-1
       
      
      !! 9.y. psi_{N+j+1} equation or x-direction
      !! 9.y.0 set coefficients { d_i } used in this equation
       
      alpha_coeffs(:) = 0.0_r_kind ;   
       
      
      alpha_coeffs(1) =   gamma_coeff(1) * c ** 2  ;                    !! d1 =   gamma_1 c^2 
      alpha_coeffs(2) = - b( phi_id ) * b( phi_id + 1 )  ;              !! d2 = - b_j * b_{j+1} 
      alpha_coeffs(3) = - gamma_coeff(2) * c ** 2  ;                    !! d3 = - gamma_2 c^2
      alpha_coeffs(4) =   b( phi_id ) * b( phi_id + 1 ) * c *  gamma_coeff(3) ;  !! d4 = b_j * b_{j+1} c gamma_3
      if( is_edge_vertical ) then                                                !! surface_bc_type = 2
          alpha_coeffs(5) = - gamma_coeff(2) * c * b( phi_id ) * b( phi_id + 1 ) * nx * alpha / ( alpha + 2.0_r_kind * mu)  ;  
      end if 
      alpha_coeffs(6) = - gamma_coeff(2) * c ;
      
      !! 9.y. row = psi_{N+j+1} : for  1 <= j <= M-1: 
      !! index: 1( phi_0 ) + N( phi_j ) + 1( phi_bar ) + 1( phi_{N+1} ) + j:  begins from ( N + 2 + j ) * 2 * n_bc_nodes
      
      row = (the_model_description%order_of_foundation_higher_order_bc(1) + 2 + phi_id) * 2 * n_bc_nodes ; !! phi_{N+j+1}
      row = row  + n_bc_nodes ; 
      
      
      !! 9.y.1. add to mass matrix ( row = psi_{N+j+1} , col = psi_{N+j} , psi_{N+j+1} , psi_{N+j+2} )
      col = row - 2 * n_bc_nodes ; !! for psi_{N+j}
      elemental_mass_matrix(row+1:row+3 , col+1:col+3) = b( phi_id + 1 ) * alpha_coeffs(1) * Mmatrix ;
      elemental_mass_matrix(row+1:row+3 , row+1:row+3) = ( b(phi_id) + b( phi_id + 1 ) ) * & 
                                                       &   alpha_coeffs(1) * Mmatrix ;
                 
      !! 9.y.1.2. if N + j + 2 < N + M+1 or j < M-1 then psi_{N+j+2} appears.
      if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(2) -1 )  then
          col = row + 2 * n_bc_nodes ; 
          elemental_mass_matrix( row+1:row+3 , col+1:col+3 ) = b(phi_id) * alpha_coeffs(1)  * Mmatrix ; !!  phi_{N+j+2}
      end if
        
        
      !! 9.y.2. damping matrix: 
      col = row - 2 * n_bc_nodes ; !! for psi_{N+j}
      elemental_damping_matrix( row+1:row+3 , col+1:col+3 ) =   b( phi_id + 1 ) * alpha_coeffs(6) * B14matrix ;
      
      elemental_damping_matrix( row+1:row+3 , row+1:row+3 ) = (b(phi_id) + b(phi_id+1)) * alpha_coeffs(6) * B14matrix ;
       
      if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(2) -1 )  then
          col = row + 2 * n_bc_nodes ; 
          elemental_damping_matrix(row+1:row+3 , col+1:col+3) =  b(phi_id) * alpha_coeffs(6) * B14matrix ;
      end if
      
      
      !! 9.y.3. add to stiffness matrix ( row = psi_{N+j+1} , col = phi_{N+j} , psi_{N+j} , phi_{N+j+1} , phi_{N+j+2} , psi )
      
      col = row - 3 * n_bc_nodes ;  !! for phi_{N+j}
  elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) = - alpha_coeffs(4) * Smatrix &
                                                        & - alpha_coeffs(5) * B23matrix ;   !! phi_{N+j}
       
      col = row - 2 * n_bc_nodes ;  !! for psi_{N+j}
      elemental_stiffness_matrix(row+1:row+3 , col+1:col+3) =   &           
   &  b(phi_id) * alpha_coeffs(2) * Mmatrix + b(phi_id+1) * alpha_coeffs(3) * Kmatrix  ;     !! psi_{N+j}
     
   
      
      elemental_stiffness_matrix(row+1:row+3 , row+1:row+3) = ( b(phi_id+1) + b(phi_id) ) *  & !! psi_{N+j+1}
    & ( -alpha_coeffs(2) * Mmatrix + alpha_coeffs(3) * Kmatrix ) ; 

      !! row = psi_{N+j+1} , col = phi_{N+j+2} , psi_{N+j+2}
      if( phi_id  < the_model_description%order_of_foundation_higher_order_bc(2) -1 )  then
          
          col = row +  n_bc_nodes ;   !! for phi_{N+j+2}
          elemental_stiffness_matrix( row+1:row+3 , col+1:col+3 ) = alpha_coeffs(4)  * Smatrix &
                                                                & + alpha_coeffs(5) * B23matrix ;  
          
          col = row + 2 * n_bc_nodes ;   !! for psi_{N+j+2}
          elemental_stiffness_matrix( row+1:row+3 , col+1:col+3 ) =  &
      &   ( b(phi_id+1) * alpha_coeffs(2)  * Mmatrix + b(phi_id) * alpha_coeffs(3)  * Kmatrix ) ; !!  phi_{N+j+2}
      
          
      end if
      
           
      
    
   end do !! over phi_id = 1:M-1 
   
 
    
  !! send to global matrices
  
    call  add_to_global_matrices( elemental_mass_matrix , 1 , row_variable_set , num_row , & 
                                & col_variable_set , num_col  ) ;
                               
    call  add_to_global_matrices( elemental_stiffness_matrix , 2 ,row_variable_set, num_row , & 
                                & col_variable_set , num_col  ) ;
   
    call  add_to_global_matrices( elemental_damping_matrix , 3 , row_variable_set , num_row , & 
                                & col_variable_set , num_col  ) ;
   
   
  

  end subroutine HW_model_foundation_truncation_boundary
   
   
 
 
 
 
 
  
  end submodule smod_found_trunc_boundary



