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
!Note that module USE is transitive, if A uses B and B uses C you don't 
!have to also declare that A uses C (though if you've renamed entities or 
!specified ONLY clauses you'll have to make sure what is transitive in a particular case).



!! Assumption 1: index of nodes begins from 1 
!! Assumption 2: an element can contain zero for absent nodes
!! Assumption 3: not all members of an element are zero, i.e. not all nodes are absent

program main

    use mod_analysis , only : start_analysis
    use mod_utils , only : r_kind 
    implicit none
    
    real(kind = r_kind) :: start_time, stop_time
                    
    call  cpu_time( start_time )
    
    call read_input_data() 
    
    
    call start_analysis() 
    
    call cpu_time( stop_time )
    
    write(*, *) 'Elapsed time, s : ',  (stop_time - start_time)
     
   
    
    

  contains
     
     
     

     
!===============================================================================!
!   Read input data                                                             !
!              
! 
! 
!      n_dyn = 1    "LINEAR STATIC ANALYSIS [LSA]" 
!      n_dyn = 2    "LINEAR DYNAMIC ANALYSIS [LDA]" 
!      n_dyn = 3    "NONLINEAR STATIC ANALYSIS [NLSA]" 
!      n_dyn = 4    "NONLINEAR DYNAMIC ANALYSIS [NLDA]"
!      n_dyn = 5    "EIGEN VALUES AND VECTORS ANALYSIS [EGVVA]"
!      n_dyn = 6    "SPECTRAL ANALYSIS [SPA]"
!      n_dyn = 7    "MODAL ANALYSIS IN TIME DOMAIN [MATD]"
!      n_dyn = 8    "MODAL ANALYSIS IN FREQUENCY DOMAIN [MAFD-COMPLETE]"
!      n_dyn = 9    "MODAL ANALYSIS IN FREQUENCY DOMAIN [MAFD-TRANSFER FUNCTIONS]"
!      n_dyn = 10   "DIRECT ANALYSIS IN FREQUENCY DOMAIN [DAFD-TRANSFER FUNCTIONS]"
!      n_dyn = 11   "DIRECT ANALYSIS IN FREQUENCY DOMAIN [DAFD-COMPLETE]" 
!
!      element type: IELTYPE , 1, 4, 14 : 3D Brick with 20 nodes
!      element type: IELTYPE , 5, 6, 12 : 2D quadliraterial with 8 nodes 
!      element type: IELTYPE , 3 (interface) or others    
!                                                                              !
!===============================================================================!



subroutine  read_input_data() 
   
    use mod_utils
    use mod_geometry
    use mod_physics
    use mod_analysis
             
    implicit none
        
    character(LEN=500) :: buffer , temp_text
    integer :: n_dis, n_spring, n_loadc , n_algo , n_accel , i_sym 
    integer :: n_freq , n_mmt , ifss , ifpr , k_opt , npl
          
    integer , parameter :: input_u = 80  
    real(kind = r_kind) :: t_step , r_tol , xx , yy , zz
    integer :: state , i , j , k , id_mat , id_elem, id_bc , n_mat, n_elem, n_bc, id_region, r  
    real(kind = r_kind )    :: array_real(12) = 0.0
    integer :: array_int(20), array_int_2(8), num_group_element(3)
    real(kind = r_kind) :: T1, F1, T2, F2, T3, F3, T4, F4, T5, F5, T6, F6
     
     
     
          

    ! Open input file 
    open ( unit   = input_u  , file    = "./input/INPUT.TXT"  , action  = 'read'   ,  &
            access = 'stream' , iostat  = state             , form    = 'formatted' )
                                  
    if ( state /= 0) stop 'Error: cannot open input file in preprocessing step. '
                
 
!!====================================================================================================! 
!!====================================================================================================!
       
    read( input_u    ,'(A)')  buffer
    
    !! check if it is TR dam or Pine Flat dam
    if( index( buffer, 'TRIANGULAR' ) /= 0 ) then
        is_TR_dam = .true.
    else if( (index( buffer, 'PINE' ) /= 0)  .or. (index( buffer, 'FLAT' ) /= 0) ) then
        is_TR_dam = .false.
    else
        close( input_u ) ;
        stop 'Error: input file should be either pine flat or triangular dam '
    end if
    
    
    
    read( input_u  ,'(A)')  buffer
    read( input_u  ,'(6I12)')  analysis_kind, num_dim , num_node , num_degree_of_freedom , &
                                num_substructure , n_dis
         
    
    
    
    read( input_u  ,'(A)')  buffer
    read( input_u  ,'(2I12,F12.7,3I12)') n_spring, n_loadc ,t_step , n_algo , n_accel , i_sym 
    time_step = t_step ;                                                !!  from analysis module
    num_load_cases = n_loadc  ;                                         !!  from physics module
    newmark_accel_type = n_accel ;
    
    if( ( n_accel /= 1 ) .and. ( n_accel /= 2 ) ) then 
            close( input_u )
            stop ' Error in input file: the NACCEL parameter should be either 1 or 2.' 
    end if
        
    read( input_u  ,'(A)')  buffer
    read( input_u  ,'(4I12,F12.7,I12)') n_freq , n_mmt , ifss , ifpr , r_tol, k_opt
    read( input_u  ,'(A)')  buffer
    read( input_u  ,'(F12.2,3I12,F12.3)') hysteresis_damping_coeff , report_node_index , &
                                 earthquake_direction  , reference_node_index , frequency_step


!!==================================================================================================!
!!
!! Set the physical model in bulk and truncation boundaries
!!
!!==================================================================================================!

    read( input_u  ,'(A)')  buffer   !! MODEL PRBLEM AT TRUNCATION BOUNDARIES[ ( i , j ) : i = 0 (SOMERFEL) , i = 1 (HW) , i = 2(GN) ; j = 0 (LYSMER) , j = 1 (HW) ]   
         
    array_int_2 = 0
    read( input_u  ,'(2I12)') array_int_2(1:2)
         
         
    the_model_description%accustic_bulk_model_id           = 0
    the_model_description%elastic_bulk_model_id            = 0
    the_model_description%reservior_inf_boundary_model_id  = array_int_2(1)
    the_model_description%foundation_inf_boundary_model_id = array_int_2(2)
    
         
    ! values up to this version: 0 for somerfeld, and 1  for HW model, and 2 for GN model 
    if( (the_model_description%reservior_inf_boundary_model_id > 2 ) .or.  &
      & (the_model_description%reservior_inf_boundary_model_id < 0 )  ) then  
        close( input_u )
        write(*,*) ' The demanded higher order bc model has not been covered uo to the current version'
        write(*,*) ' Up to now, the supported models for reservior truncation boundary are: '
        write(*,*) ' 0 (Somerfeld) , 1 (HW) , 2 (GN) '
        stop
    end if
         
         
    ! values up to this version: 0 for Lysmer. 
    if( ( the_model_description%foundation_inf_boundary_model_id > 1 ) .or.  &
    &   ( the_model_description%foundation_inf_boundary_model_id < 0 )  ) then 
        close( input_u )
        write(*,*) ' The demanded higher order bc model has not been covered up to the current version'
        write(*,*) ' Up to now, the supported models for foundation truncation boundary are: 0 (Lysmer)'
        stop 
    end if
         
         
    read( input_u  ,'(A)')  buffer !! RESERVIOR TRUNCATION BOUNDARY MODEL: order: O(N,M) , {a_j }, {b_j }
         
    !! check if reservior truncation boundary is higher order model
    if( the_model_description%reservior_inf_boundary_model_id > 0 ) then  ! higher order bc model
             
            !! read the order
            read( input_u  ,'(2I12)') array_int_2(1:2)
            the_model_description%order_of_reservior_higher_order_bc = array_int_2(1:2)
             
             
            !! assign memory for stroring coefficients \{ a_j \} and \{ b_j \}
            call  the_model_description%initialize_model_description( boundary_id = 1 )
             
            read( input_u  ,'(A)')  buffer    ! read line containing coefficient \{ a_j \} as string
             
            do  i = 1 , size( the_model_description%reservior_incident_coefficient )  ! N for GN and N+1 for HW
                read( buffer , * , iostat = state) the_model_description%reservior_incident_coefficient(1:i)
                if( state /= 0 )then
                     
                    if( the_model_description%reservior_inf_boundary_model_id .eq. 1 ) then    !! HW model
                        write(*,*) 'The order of HW model for reservior truncation boundary is:' ,  &
                            & the_model_description%order_of_reservior_higher_order_bc
                    elseif( the_model_description%reservior_inf_boundary_model_id .eq. 2 ) then    !! GN model
                        write(*,*) 'The order of GN model for reservior truncation boundary is:' ,  &
                            & the_model_description%order_of_reservior_higher_order_bc
                    end if
                     
                    write(*,*) 'The number of coefficients a_j in this model should be: ' , &
                               & size( the_model_description%reservior_incident_coefficient )
                    write(*,*) ' Bad Input Error.'
                     
                    close( input_u )
                    stop
                end if 
            end do
             
             
            if( the_model_description%order_of_reservior_higher_order_bc(2) > 0 ) then
                read( input_u  ,'(A)')  buffer    ! read line containing coefficient \{ b_j \} as string
                do  i = 1 , size( the_model_description%reservior_evanescent_coefficient )  ! N for GN and HW
                     read( buffer , * , iostat = state ) the_model_description%reservior_evanescent_coefficient(1:i)
                     if( state /= 0 )then
                         if( the_model_description%reservior_inf_boundary_model_id .eq. 1 ) then    !! HW model
                             write(*,*) 'The order of HW model for reservior truncation boundary is:' ,  &
                                      & the_model_description%order_of_reservior_higher_order_bc
                         elseif( the_model_description%reservior_inf_boundary_model_id .eq. 2 ) then    !! GN model
                             write(*,*) 'The order of GN model for reservior truncation boundary is:' ,  &
                                     & the_model_description%order_of_reservior_higher_order_bc
                         end if
                     
                         write(*,*) 'The number of coefficients b_j in this model should be: ' , &
                                   & size( the_model_description%reservior_evanescent_coefficient ) 
                         write(*,*) ' Bad Input Error.'
                     
                         close( input_u )
                         stop
                     end if 
                end do  ! over i
            end if
        
    end if
        
    read( input_u  ,'(A)')  buffer !! FOUNDATION TRUNCATION BOUNDARY MODEL: order: O(N,M) , {a_j }, {b_j } 
        
        !! check if foundation truncation boundary is higher order model
        if( the_model_description%foundation_inf_boundary_model_id > 0 ) then  ! higher order bc model
             
            !! read the order
            read( input_u  ,'(2I12)') array_int_2(1:2)
            the_model_description%order_of_foundation_higher_order_bc = array_int_2(1:2)
             
            !! assign memory for storing coefficients \{ a_j \} and \{ b_j \}
             
            call  the_model_description%initialize_model_description( boundary_id = 2 )
             
            read( input_u  ,'(A)')  buffer    ! read line containing coefficient \{ a_j \} as string
            do  i = 1 , size( the_model_description%foundation_incident_coefficient )  ! N for GN and N+1 for HW
                read( buffer , * ) the_model_description%foundation_incident_coefficient(1:i)
            end do
              
            if( array_int_2(2) > 0 ) then
                read( input_u  ,'(A)')  buffer    ! read line containing coefficient \{ b_j \} as string
                do  i = 1 , size( the_model_description%foundation_evanescent_coefficient )  ! N for GN and N+1 for HW
                    read( buffer , * ) the_model_description%foundation_evanescent_coefficient(1:i)
                end do
            end if 
             
        end if
       
         
    
         
    
   
     
     

       
!!====================================================================================================! 
!!      reading coordinate of nodes
!!====================================================================================================!
   
         
    if( num_dim .eq. 2 ) then
        max_num_node_in_element = 8    ! Q8 element : variables in geometry module
        num_element_around_node = 4
        
    else
        max_num_node_in_element  = 20  ! Brick element
        num_element_around_node  = 8
        
    end if
         
    call  initialize_node_matrix()
     
    if( num_dim .eq. 2 ) then
        temp_text = '(I12,2F12.4)'
    else
        temp_text = '(I12,3F12.4)'
    end if
         
    !!  read nodes
    read( input_u  ,'(A)')  buffer   !! reads: I     X   Y  Z  coordinates
        
   
    do  i = 1 , num_node
        read(  input_u    , temp_text ) k , node(k , 1:num_dim)        
    end do
         
    
         
!!!====================================================================================================! 
!!!       reading dof matrix
!!====================================================================================================!
         
    
     
         
    read(  input_u  ,'(A)') buffer
         
    do i = 1 , num_node
        read(  input_u   , '(8I12)')  k , array_int_2(1:7)   ! as node constraint
        
        do j = 1 , num_dim                                              ! first col of dof_matrix is class_id
           if( array_int_2(j ) .eq. 0 ) then
               degree_of_freedom(k , 2*j ) = 1                          !  translatioal dof
           end if                                                       
        end do 
                
        if( num_degree_of_freedom > num_dim ) then
            if( array_int_2(7) .eq. 0 ) then
                degree_of_freedom(k , 2 * num_degree_of_freedom ) = 1   !  pressure dof (if any)
            end if
        end if
    end do ! over i
         
           
!!====================================================================================================! 
!!      reading force: enforced displacement
!!====================================================================================================!  
! this is not used in freq-domain analysis     
     
    read(  input_u , '(A)')  buffer
    
    
    
    if( n_dis /= 0 ) then
        enforced_displacement%num_nodes =  n_dis ;                  !!  from physics module
        call   enforced_displacement%initialize_enforced_disp() ;
        do i = 1 , n_dis
            read(  input_u   , '(I12,3F12.3)') enforced_displacement%node_ids(i)  ,  &
                                           &  enforced_displacement%disp_vector(i,1:3)
        end do
    end if
         

!!====================================================================================================! 
!!     reading force : spring support
!!====================================================================================================!
! this is not used in freq-domain analysis  unless it is DRM layer

    read(  input_u   , '(A)') buffer
         
    
    if( n_spring /= 0 ) then
        boundary_springs%num_fixed_nodes = n_spring ;               !!  from physics module
        call   boundary_springs%initialize_boundary_spring() ;
        do  i = 1 , n_spring
            read(  input_u   , '(2I12,6E12.3)') k , j , array_real(1:6) ;
            boundary_springs%fixed_node_ids(i) = k ;
            boundary_springs%k_type(i) = j         ;
            boundary_springs%disp_vector(i,1:3)      = array_real(1:3)  ;
            boundary_springs%spring_stiffness(i,1:3) = array_real(4:6)  ;
        end do
    end if
    
    
    
!!====================================================================================================! 
!!        read physical info of substructures
!!====================================================================================================!
! first read size of matrices and arrays
       
    num_element         = 0
    num_material        = 0
    num_boundary_load   = 0
          
    inquire( input_u , pos = k )  ! get pointer to the current position
     
     
    n_bc = 0  ! used to specify max nuld among different substructures
        
    do  i = 1 , num_substructure
        read(  input_u   ,*  ) 
        read(  input_u   ,'(7I12)') array_int_2(1:7) ! IGROUP,IELTYPE,NEL,NCRIT,NMAT,NULD,JINERTIA
        
        num_element = num_element + array_int_2(3)  
        num_material = num_material + array_int_2(5)
        num_boundary_load = num_boundary_load + array_int_2(6)   

        
        n_bc = max( n_bc ,  array_int_2(6) )
        
        read(  input_u , * )
              
        do j = 1 , array_int_2(5)  ! n_mat
           read(  input_u  , *)    ! PMAT(IMAT,1:6),GRV
        end do ! over j
             
        read(  input_u   , * ) 
              
        do j = 1 , array_int_2(6) ! n_uld
           read(  input_u   , *)   
        end do 
             
        read(  input_u    , * ) 
        

              
        if( array_int_2(2) .eq. 3  ) stop 'Error: element type is not covered yet.'
                    
        do j = 1 , array_int_2(3)  ! NEL     ! IEL,IELCONNECT(1:20),IMAT,NINT,IULDEL(1:6)
           read(  input_u  , * ) 
        end do         
    end do ! over i
          
         
   
          
!! Second, assign suitable memory to matrices
          
    call  initialize_element_matrix()
    call  initialize_physics( n_bc )

    
!! Third, rewind to previously marked point and read again
      
    read(input_u ,'(A)', advance ='NO', pos = k )       ! return to pointer position 
          
    
    
    id_elem = 0
    id_mat  = 0
    id_bc   = 0
          
    num_group_element(:) = 0
          
    do  i = 1 , num_substructure
        
        read(  input_u   ,'(A)'   ) buffer
        read(  input_u   ,'(7I12)') array_int_2(1:7) ! IGROUP,IELTYPE,NEL,NCRIT,NMAT,NULD,JINERTIA
        
        
        n_elem   = array_int_2(3) ;
        n_mat    = array_int_2(5) ;
        n_bc     = array_int_2(6) ;
        
        if( i .eq. 1 ) then
            dam_jinertia = array_int_2(7) ;
        end if
                   
        read(  input_u   ,'(A)'   ) buffer   ! IMAT  E    POI    GAMMA    N   SIGY0    HPRIM  GRV 
               
        do  j = 1 , n_mat
            id_mat = id_mat + 1
            read(input_u , '(I12, E12.6,F12.3,F12.6,4F12.3)') k    , &
                                material( id_mat )%elastic_module  , & 
                                material( id_mat )%poisson_ratio   , &
                                material( id_mat )%special_weight  , &
                                material( id_mat )%n_factor        , &
                                material( id_mat )%sigma0_y        , &
                                material( id_mat )%h_prime         , &
                                gravitational_constant
        end do ! over j
               
        read(  input_u   , '(A)' ) buffer
                
        !  IULD    KTYPE      PR    NFACE    LDIR     XIREF
        do j = 1 , n_bc
           id_bc = id_bc + 1
           read(  input_u   , '(I12,F12.0,F12.3,2F12.0,F12.3)') k, array_real(1:5)  ! IULD,PULD(IULD,1:5)
                  
           boundary_load( id_bc )%k_type           =  int(  array_real(1)  )
           boundary_load( id_bc )%strength         =  real( array_real(2) , r_kind )
           boundary_load( id_bc )%face_id          =  int(  array_real(3)  )
           boundary_load( id_bc )%direction        =  int(  array_real(4)  )
           boundary_load( id_bc )%reference_height =  real( array_real(5) , r_kind )
        end do ! over j
            
       
             
        read(  input_u    , '(A)' ) buffer
        
        id_region  = i
        element_type(i) = array_int_2(2)
        
        num_group_element(i) = n_elem ;
        
        if( i .eq. 2 ) then  ! check if  foundation is absent
            if( ( num_substructure .eq. 2 ) .and. ( num_dim  /= num_degree_of_freedom ) ) then   
                  id_region = 3
                  element_type(3) = element_type(2)
                  element_type(2) = 0
                  num_group_element(3) = num_group_element(2) ;
                  num_group_element(2) = 0 ;
            else
               is_complete_mass_model = array_int_2(7) .eq. 1
            end if 
        end if    
        
        
        do j = 1 , n_elem                           ! IEL,IELCONNECT(1:20),IMAT,NINT,IULDEL(1:6)
            id_elem = id_elem + 1
            read(  input_u  , '(I12,22I6,6I3)' ) k  , array_int(1:20) , array_int_2(1:2), array_int_2(3:8)
           
            element_property( id_elem , 1 ) = id_region                       ! region_id
            element_property( id_elem , 2 ) = id_mat - n_mat + array_int_2(1) ! material_id  
            element_property( id_elem , 3 ) = array_int_2(2)                  ! num_gauss_points or nint
            
            do k = 1 , n_bc                                                   ! column 4 to 4+ nuld = boundary condition
               if( array_int_2( k + 2) /= 0 ) then
                   element_property( id_elem , 3 + k ) = id_bc - n_bc + k 
               end if
            end do
            
            element_matrix( id_elem , 1:max_num_node_in_element ) = array_int(1:max_num_node_in_element)
          
        end do ! loop over j
          
    end do ! over i

!!!====================================================================================================!
!!!     Read elements in DRM layer
!!!====================================================================================================!          
          
   
    
    read(input_u  , '(A)' ) buffer  !DRM ANALYSIS KIND ( 0 = NO DRM, 1 = FREE FIELD DATA IS USED , 2 = DIRECT DECONVOLUTION DATA IS USED=NO FREE FIELD DATA)  
    read(input_u,*) drm_layer_ndyn          
    read(input_u  , '(A)' ) buffer  !DRM LAYER ID , NUM ELEMENTS IN LAYER , NUM NODE IN LAYER
    read(input_u,*) drm_layer_info%layer_id , drm_layer_info%num_element_in_layer, drm_layer_info%num_node_in_layer 
    read(input_u  , '(A)' ) buffer  !	ID OF ELEMENTS IN DRM LAYER
    
    if( (drm_layer_ndyn < 0 ) .or. ( drm_layer_ndyn > 2 ) ) then
        stop 'Error: input file: DRM NDYN can be 0 (no DRM) ,1 (free field data) or 2 (Deconv data)'
    end if 
    
    if( (drm_layer_ndyn .eq. 0 ) .and. (drm_layer_info%layer_id /= 0 )  ) then
        stop 'DRM analysis kind is nonzero but Layer id is not specified'
    end if 
    
    if(   ( drm_layer_info%layer_id > 0 ) .and. ( drm_layer_info%num_element_in_layer > 0) .and. &
        & ( drm_layer_info%num_node_in_layer > 0 ) ) then
        
        call drm_layer_info%initialize_drm_layer() ;
        
        do  i = 1 ,  drm_layer_info%num_element_in_layer
            read(  input_u  , * ) j
            drm_layer_info%element_index(i) = num_group_element(1) + j
        end do
        
        read(input_u  , '(A)' ) buffer  !	NODE INDEX CORRESPONDENCE:
        
        do  i = 1 ,  drm_layer_info%num_node_in_layer 
            read(  input_u  , * ) drm_layer_info%paired_node_index(i,1:2)
        end do
        
    else
        read(input_u  , '(A)' ) buffer  !	NODE INDEX CORRESPONDENCE:
    end if 
   
    !! from physics module
    is_time_domian = (analysis_kind .eq. 2) ! .or. &
                !  & ((analysis_kind .eq. 10) .and. ( drm_layer_info%layer_id > 0 )) ;
    
    
    
!!!====================================================================================================! 
!!!        read forces: needed only for time domain
!!!====================================================================================================!

    
    
    if( .not. is_time_domian ) then    !! freq domain not drm analysis    
        close (input_u) ;
        return ;
    end if
    
    
    
    read(  input_u  , '(A)' ) buffer  ! CONCENTRATED NODAL LOADS IN THE STRUCTURE:
    
    read(  input_u  , '(A48,I12)' ) buffer , npl
          
    nodal_loads%num_point_load = npl ;                              !! from physics module
          
    read(  input_u  , '(A)'  )   buffer    
          
    if( npl /= 0 ) then
        call   nodal_loads%initialize_nodal_load() ;
            
        do  i = 1 , npl
            read(  input_u  , '(I12,3F12.3)' ) nodal_loads%node_ids(i) , nodal_loads%load_vector(i,1:3)
        end do     
    end if
        

!!!====================================================================================================! 
!!!        load conditions: needed only for time domain
!!!====================================================================================================!
    
    
        
        
    if( num_load_cases .eq. 0 ) then
        close (input_u) ;
        return ;
    end if 
    
    if( num_load_cases /= 1 ) then
        close (input_u) ;
        stop 'Warning: only one load case is implemented. The others are ignored.'
    end if 
              
    read(  input_u  , '(A)' ) buffer                    !   LOAD CONDITIONS FOR STRUCTURE:
    read(  input_u  , '(A)' ) buffer                    !! LOADC NU.      DEAD WEIGHTS         POINTLOAD 
    read(  input_u  , '(I18,7F18.3)') j , load_cases%load_coeffs(1:7)
    read(  input_u  , '(A)' ) buffer                    !!  LOAD CONDITIONS TIME DEPENDNCY:  TMIN  
    read(  input_u  , '(F54.2,F18.2,I18)' ) load_cases%tmin , load_cases%tmax , num_static_step

    if( time_step < 1.0e-9 ) then
        close (input_u) ;
        stop 'Warning: The time step is smaller than 1e-9.'
    end if 
    load_cases%num_rows = ceiling( 1e-9 + ( load_cases%tmax -load_cases%tmin )/time_step ) ;
    
    call load_cases%initialize_load_case() ;
    load_cases%time_load(1,1) = load_cases%tmin ;
    do  i = 2 , load_cases%num_rows              
        load_cases%time_load(i,1) = load_cases%time_load(i-1,1) + time_step ;  !! first row is time
    end do
	
	
	
    temp_text = '(F6.2,F6.3,F6.2,F6.3,F6.2,F6.3,F6.2,F6.3,F6.2,F6.3,F6.2,F6.3)' ;
	
    do  i = 1 , 7
        read(  input_u  , '(A)' ) buffer                            !! reads this !----------------!
        read(  input_u  , '(A)' ) buffer                            !!  IFUNC   NTIME
        read(  input_u  ,'(2I18)') array_int(1:2)                   !!  IFUNC , NTIME
        read(  input_u  , '(A)' ) buffer                            !!   [T    Fj     T    Fj ...]
	    
        n_bc = 0 ;
        array_real(1) = 0.0_r_kind ;
		
        if( array_int(2) > 0 ) then
            array_int(1) = ceiling( real( array_int(2) , kind = r_kind ) / 6.0_r_kind ) ;
            do  j = 1 , array_int(1)
                read(  input_u, temp_text ) T1, F1, T2, F2, T3, F3, T4, F4, T5, F5, T6, F6
                array_int_2(1) = int( abs( load_cases%tmin / time_step) ) + 1 ;
                array_int_2(2) = floor( 1e-10 + T1 / time_step)  + array_int_2(1) ;
                array_int_2(3) = floor( 1e-10 + T2 / time_step)  + array_int_2(1) ;
                array_int_2(4) = floor( 1e-10 + T3 / time_step)  + array_int_2(1) ;
                array_int_2(5) = floor( 1e-10 + T4 / time_step)  + array_int_2(1) ;
                array_int_2(6) = floor( 1e-10 + T5 / time_step)  + array_int_2(1) ;
                array_int_2(7) = floor( 1e-10 + T6 / time_step)  + array_int_2(1) ;
                load_cases%time_load( array_int_2(2) , i+1 ) = F1 ;
                load_cases%time_load( array_int_2(3) , i+1 ) = F2 ;
                load_cases%time_load( array_int_2(4) , i+1 ) = F3 ;
                load_cases%time_load( array_int_2(5) , i+1 ) = F4 ;
                load_cases%time_load( array_int_2(6) , i+1 ) = F5 ;
                load_cases%time_load( array_int_2(7) , i+1 ) = F6 ;
			    
                do  k = array_int_2(2) , array_int_2(7)
                    if ( (k > array_int_2(2)) .and. (k < array_int_2(3)) ) then
				    
                        array_real(2) = real( k - array_int_2(2) , kind = r_kind )  &
                                    & / real( array_int_2(3) - array_int_2(2) , kind = r_kind ) ;
				                    
                        load_cases%time_load( k , i+1 ) = F1 + (F2 - F1) * array_real(2) ;  
				                   
                    else if( (k > array_int_2(3)) .and. (k < array_int_2(4)) ) then
				        
                        array_real(2) = real( k - array_int_2(3) , kind = r_kind )  &
                                    & / real( array_int_2(4) - array_int_2(3) , kind = r_kind ) ;
				                    
                        load_cases%time_load( k , i+1 ) = F2 + (F3 - F2) * array_real(2) ;  
				    
                    else if( (k > array_int_2(4)) .and. (k < array_int_2(5)) ) then
				        
                        array_real(2) = real( k - array_int_2(4) , kind = r_kind )  &
                                    & / real( array_int_2(5) - array_int_2(4) , kind = r_kind ) ;
				                    
                        load_cases%time_load( k , i+1 ) = F3 + (F4 - F3) * array_real(2) ;  
				    
                    else if( (k > array_int_2(5)) .and. (k < array_int_2(6)) ) then
                        array_real(2) = real( k - array_int_2(5) , kind = r_kind )  &
                                    & / real( array_int_2(6) - array_int_2(5) , kind = r_kind ) ;
				                    
                        load_cases%time_load( k , i+1 ) = F4 + (F5 - F4) * array_real(2) ;  
				    
                    else if( (k > array_int_2(6)) .and. (k < array_int_2(7)) ) then
                        array_real(2) = real( k - array_int_2(6) , kind = r_kind )  &
                                    & / real( array_int_2(7) - array_int_2(6) , kind = r_kind ) ;
				                    
                        load_cases%time_load( k , i+1 ) = F5 + (F6 - F5) * array_real(2) ;  
				                                  
                    end if
                end do ! over k
				
                if( n_bc /= 0 ) then
                    do  k = n_bc , array_int_2(2)
                        array_real(2) = real( k - array_int_2(2) , kind = r_kind )  &
                                    & / real( array_int_2(2) - n_bc , kind = r_kind ) ;
				                    
                        load_cases%time_load(k,i+1) = F1 + (F1 - load_cases%time_load(n_bc,i+1))*array_real(2);  
                    end do
                end if
                n_bc = array_int_2(7) ;
				
            end do
        end if
    end do ! over i
	
    read(  input_u  , '(A)' ) buffer   !---------------------------------------  
	
    
    if( analysis_kind .eq. 10 ) return ; !! for freq domain no need to these nodes no matter if DRM exists or not
 
!!!====================================================================================================! 
!!!        read analysis
!!!====================================================================================================!     
          
    read(  input_u, '(A41,I3)' ) buffer , k    ! NUMBER OF NODS FOR TIME HISTORY RECORD:  
    
    if( k < 1 ) stop 'Error: input file: num of nodes for reporting history should be > 0.'
    i = 0 ;
    if( allocated( reproted_node_index_timeD ) ) deallocate( reproted_node_index_timeD, stat = i ) ;  !! from analysis module
    if( i /= 0 ) stop 'Failed in deallocating memory for reported_node_index_timeD in main' ;
    
    allocate( reproted_node_index_timeD( k ) , stat = i )   
    if( i /= 0 ) stop 'Failed in allocating memory for reported_node_index_timeD in main' ; 
     
    do  i = 1, k
        read(  input_u, '(I12)' ) reproted_node_index_timeD(i)
    end do
          
         
    read(  input_u, '(A41,I3)' ) buffer , num_node_report_time_data    !! TIMES IN WHICH THE OUTPUT IS PRINTED  :  ...
    
    if( num_node_report_time_data < 0 ) stop 'Error: input file: num of time instances for reporting history should be >= 0.'
    
    if( num_node_report_time_data .eq. 0 ) then
        close (input_u) ;
        return ;
    end if

    
    
    i = 0 ;
    if( allocated( reported_time ) ) deallocate( reported_time, stat = i ) ;  !! from analysis module
    if( i /= 0 ) stop 'Failed in deallocating memory for reported_time in main' ;
    
    allocate( reported_time( num_node_report_time_data ) , stat = i )   
    if( i /= 0 ) stop 'Failed in allocating memory for reported_time in main' ; 
    
    
    do  i = 1 , num_node_report_time_data
        read(input_u,*) reported_time(i)
    end do
          
          
    close (input_u) 
    
   
   
      
     
end subroutine read_input_data






end program

