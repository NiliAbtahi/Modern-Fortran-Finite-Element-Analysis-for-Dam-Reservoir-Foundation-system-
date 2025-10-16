!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!

    module  mod_analysis

        use mod_utils , only :  r_kind
        
        implicit none
        
        ! parameters used in frequency-domain analysis
        real( kind = r_kind ) , public  :: frequency_step              ! DOMEGA
        real( kind = r_kind ) , public  :: hysteresis_damping_coeff    !HYSDP
        
        
        integer , public :: report_node_index         ! NJTF
        integer , public :: reference_node_index      ! NRTF
        integer , public :: num_node_report_time_data = -1 ;
        integer , allocatable , public  :: reproted_node_index_timeD(:)  !! Node at which time domain data is reported
        real( kind = r_kind ) , allocatable, public  :: reported_time(:) !! Time instances that data is reported 
         
        
        ! parameters used in time-domain analysis
        real( kind = r_kind ) , public  :: time_step              ! Timestep
        
        !! functions and subroutines
        public   :: start_analysis
       
        private  :: free_unused_memory
        private  :: freq_domain_analysis
        private  :: time_domain_analysis
        private  :: store_in_file
        private  :: find_variable_identifier
        private  :: fft_transformation
        
        private  :: clean_up_analysis
                
            
  
    contains
  


!!=================================================================
!!
!!   start_analysis subroutine
!!
!!================================================================= 
  
    subroutine start_analysis()

        use mod_physics  , only : prepare_physical_data , analysis_kind ,drm_layer_info
        use mod_fem      , only : do_fem_analysis 
        implicit none   

        
        
        call  prepare_physical_data()  

        call  do_fem_analysis()
        
       
        ! call  store_in_file()
       
       
        if( analysis_kind .eq. 10) then
            if( drm_layer_info%layer_id .eq. 0 ) then
                call  free_unused_memory()  ;
            end if
            call  freq_domain_analysis() ;
        elseif( analysis_kind .eq. 2 ) then
            call  time_domain_analysis() ;
        else
            stop 'Error: only NDYN=2 (Linear Dynamic) and NDYN=10 (Freq) analysis are allowed' ; 
        end if
      
        ! call  store_in_file()
       
        call  clean_up_analysis()
       
       
    end subroutine start_analysis
  
  


!!====================================================================!!
!!
!!  time_domain_analysis
!!
!!
!!====================================================================!!

    subroutine  time_domain_analysis()
       
        use  mod_physics  , only : mass , stiffness , damping , force , dictionary , dependency 
        use  mod_physics  , only : translational_variable, enforced_displacement, boundary_springs
        use  mod_physics  , only : nodal_loads , load_cases, get_variable_location , newmark_accel_type
        use  mod_physics  , only : dead_weight_force , surface_pressure_force , num_static_step
        use  mod_physics  , only : tension, compression, element_property , the_model_description
        use  mod_physics  , only : is_dam_present, is_reservior_present, is_foundation_present
        use  mod_physics  , only : drm_layer_info , is_TR_dam , analysis_kind , u_uddot_njtf_drm
        use  mod_physics  , only : find_DRM_internal_face , drm_layer_ndyn
        use  mod_geometry , only : num_dim , degree_of_freedom , node, num_degree_of_freedom, num_node
        use  mod_solver   , only : initialize_real_solver , clean_up_solver , solve_real_system, static_solution
        use  mod_solver   , only : matrix , response , solution, u_vec, udot_vec, uddot_vec, temp_real_vec
        use  mod_fem      , only : find_strain_stress , set_gauss_points, gauss_points , force_in_drm_layer
        implicit none
       
        integer :: ffa_file_id , deconv_data_x_file_id , deconv_data_y_file_id , max_step 
        integer :: i , j , k , n, d, r , s , istep, num_var, num_region_element(3)
        real( kind = r_kind ) :: tmp_re  
        real( kind = r_kind ) :: alpha_k  , alpha_m  ;  !! damping coefficients
        real( kind = r_kind ) :: gamma_newmark  , beta_newmark , point_data_vector(8) ; 
        real( kind = r_kind ) :: alpha0, alpha1, alpha2, alpha3, alpha4, alpha5 , large_stiff ;  
        logical ::  is_data_printed ; 
        character(LEN=512) :: buffer, ff_file_name
        
        
        
        
        !! set Rayley damping coefficients
        if( (.not. is_reservior_present ) .and. (.not. is_foundation_present )) then  !! only dam
            alpha_k = 0.0013390_r_kind
            alpha_m = 1.4535750_r_kind
        elseif( (.not. is_reservior_present ) .and. is_foundation_present ) then  !! dam & foundation
            alpha_k = 0.0013390_r_kind
            alpha_m = 1.4535750_r_kind
        elseif( is_reservior_present  .and. (.not. is_foundation_present ) ) then  !! dam and reservoir
            alpha_k = 0.0012790_r_kind
            alpha_m = 1.2719280_r_kind
        elseif( is_reservior_present  .and.  is_foundation_present ) then   !! dam-foundation reservoir
            alpha_k = 0.0012790_r_kind
            alpha_m = 1.2719280_r_kind
        else
            stop 'Solver: structure does not lie in assumed possibilities. '
        end if
        
        max_step = 1310  ;  !load_cases%num_rows ; ! 1310
        
        !! num elements in each substructure
        num_region_element(:) = 0 ;
        do i = 1 , size( element_property ,1 )
           n = element_property(i,1) ;
           num_region_element(n) = num_region_element(n) + 1 ;
        end do
        
        
        if( newmark_accel_type .eq. 1 ) then
            gamma_newmark = 0.50_r_kind ;
            beta_newmark  = 0.250_r_kind ;
        else if( newmark_accel_type .eq. 2 ) then
            gamma_newmark = 0.50_r_kind ;
            beta_newmark  = 1.0_r_kind / 6.0_r_kind ;
        else
            stop 'Bad input: Newmark Naccel should be either 1 or 2. '
        end if 
        
        alpha0 = 1.0_r_kind      / ( beta_newmark * ( time_step ** 2 ) ) ;
        alpha1 = gamma_newmark   / ( beta_newmark * time_step ) ;
        alpha2 = 1.0_r_kind      / ( beta_newmark * time_step ) ;
        alpha3 = 0.50_r_kind     / beta_newmark   - 1.0_r_kind  ;
        alpha4 = gamma_newmark   / beta_newmark   - 1.0_r_kind  ;
        alpha5 = ( gamma_newmark / (2.0_r_kind * beta_newmark ) - 1.0_r_kind ) * time_step ;
        
        
        large_stiff = 100.00000D0
        large_stiff = large_stiff **3
        large_stiff =9999999999.9999999999999999999999999999999D0*(large_stiff**3)
        
        
        
        !! step 0. Add Rayleigh damping to damping matrix
        !! 0.1. specify nodes with translational dof
        
        allocate( translational_variable( size( dictionary ) - 1 ) , stat = k )
        if( k /=0 ) stop 'Error in allocating variable class array in analysis module '
     
        translational_variable(:) = 0                                   !! non-zero value for vars with translational dof
     
        do  i = 1 , size( degree_of_freedom , 1 )                       !! over nodes 
            if( degree_of_freedom( i , 1 ) .eq. 0 ) cycle
            do  j = 1 , num_dim
                if( degree_of_freedom( i , 2 * j ) .eq. 0 ) cycle
                translational_variable( degree_of_freedom( i , 2 * j )) = 1 ;
            end do ! over j
        end do ! over i
     
        !! 0.2. Form Rayleigh damping
        do  k = 1 , size( translational_variable )
            if( translational_variable( k ) .eq. 0 ) cycle              !! 1.2. The variable = eq has no translational dof
            do  j = dictionary( k ) + 1 , dictionary( k + 1 )           !! take k-th equation
                if ( translational_variable( dependency(j) ) .eq. 0 )   cycle
                   damping(j) = damping(j) + alpha_m * mass(j) + alpha_k * stiffness(j) ; !! 1.3. Add Rayleigh damping
            end do ! over j
        end do ! over k
                
        
        if( allocated(  translational_variable )) then
            deallocate( translational_variable , stat = k )
            if( k /= 0 ) stop ' Error in deallocating arrays in fem module'
        end if
        
    
    !! if DRM model, then open free field analysis (ffa) solution
    ffa_file_id = 0 
    deconv_data_x_file_id  = 0 ;
    deconv_data_y_file_id  = 0 ;
    
    if( drm_layer_ndyn > 0 ) then
        
        !! find internal face of each element
        !! find min and max x with max height
        call   find_DRM_internal_face()
        
        
        !! open file containing free field data
        if( is_TR_dam ) then
            
            if( drm_layer_ndyn .eq. 2 ) then
                stop 'Error: DRM from Deconv data is designed only for PF DAM: analysis suroutine. '
            end if 
            
            write(buffer,'(I2.2)') drm_layer_info%layer_id
            ff_file_name = './Data/TR10B-LAYER='//trim(buffer)//'.TXT'
                
            ffa_file_id = 5 
            
            open( unit = ffa_file_id , file = ff_file_name , action = 'read', access='stream',  &
                & iostat = k, form = 'formatted')
                
            if ( k /= 0) stop 'Error: cannot open free field data in freeFieldData analysis suroutine. ' 
        else
            
            if( drm_layer_ndyn .eq. 1 ) then
                write(buffer,'(I2.2)') drm_layer_info%layer_id
                ff_file_name = './Data/PF8B-LAYER='//trim(buffer)//'.TXT'
                
                ffa_file_id = 5 
                open( unit = ffa_file_id , file = ff_file_name , action = 'read', access='stream',  &
                    & iostat = k, form = 'formatted')
                
                if ( k /= 0) stop 'Error: cannot open free field data in freeFieldData analysis suroutine. ' 
                
                read( ffa_file_id , *) buffer ! Layer-ID    ,    Num-DRM-Nodes    ,   Num-Steps   ,   DeltaT
                read( ffa_file_id , *) i,j,k, tmp_re
                if( (i /= drm_layer_info%layer_id) .or. (j /= drm_layer_info%num_node_in_layer) .or. &
                    & ( abs(tmp_re - time_step ) > 1.0e-10 ) ) then
                stop 'Error: incompatible Layer-Id or numNodeInLayer or time-step inn freeFieldData in TimeDSolver. '
                end if
        
                max_step = min( max_step , k )  !! there is no more data in FreeFieldData
                read( ffa_file_id , *)   !! DRM Node info:
                do  k = 1 , drm_layer_info%num_node_in_layer
                    read( ffa_file_id , *)   !! DRM Node info:
                end do 
        
                read( ffa_file_id , * )  !! ID  ,   UX  ,   UY  ,   dotUX   ,   dotUY   ,   ddotUX   ,  ddotUY
        
            else
                
                !! open X dir if earthquake is in X dir 
                if( abs( load_cases%load_coeffs(5) ) > 1.0e-10 ) then  ! there is EQ in X direction
                    
                    deconv_data_x_file_id = 8 ;
                    open( unit = deconv_data_x_file_id , file = './Data/PF-EQ-X-DECONV-DIFF-LEVELS.TXT' ,  &
                        &    action = 'read', access='stream', iostat = k, form = 'formatted')
                
                    if ( k /= 0) stop 'Error: cannot open Deconv-x-data in time domain analysis. '
                    
                    !! read levels
                    read( deconv_data_x_file_id , *)   !! NUM LEVELS AT WHICH EQ WAS DECONVOLUTED:
                    read( deconv_data_x_file_id , *)   drm_layer_info%num_levels , tmp_re
                    
                    if( abs( tmp_re - time_step ) > 1.0e-6 ) then
                        stop 'Error: Inconsistent time steps in Deconv-x-data and input:time domain analysis. '
                    end if
                    
                    drm_layer_info%num_levels = drm_layer_info%num_levels + 1 ! to include zero level
                    
                    if( drm_layer_info%num_levels .eq. 1 ) then
                        stop 'Error: Num levels in DRM-X-Deconv_Data should be >0 time domain analysis. '
                    end if
                    
                    allocate( drm_layer_info%levels( drm_layer_info%num_levels ) , stat = k )
                    if ( k /= 0) stop 'Error: in allocating memory for DRM levels in time domain analysis. '
                    
                    drm_layer_info%levels(1:drm_layer_info%num_levels)  = 0 ! zero value is free surface
                    
                    read( deconv_data_x_file_id , *)   ! LEVELS AT WHICH EQ WAS DECONVOLUTED:
                    read( deconv_data_x_file_id , *) drm_layer_info%levels(2:drm_layer_info%num_levels)
                    read( deconv_data_x_file_id , *)   ! NUM ROWS AND NUM COLUMNS:
                    read( deconv_data_x_file_id , *) drm_layer_info%num_rows , drm_layer_info%num_cols
                    read( deconv_data_x_file_id , *)   ! TIME    A-L1   V-L1   U-L1    A-L2   V-L2    U-L2
                    max_step = min( max_step , drm_layer_info%num_rows )  !! there is no more data in DeconvX-file
                end if
                
                if( abs( load_cases%load_coeffs(6) ) > 1.0e-10 ) then  ! there is EQ in Y direction
                    
                    deconv_data_y_file_id = 9 
                    open( unit = deconv_data_y_file_id , file = './Data/PF-EQ-Y-DECONV-DIFF-LEVELS.TXT' ,  &
                        &    action = 'read', access='stream', iostat = k, form = 'formatted')
                
                    if ( k /= 0) stop 'Error: cannot open Deconv-y-data in time domain analysis. '
                    
                    !! read levels
                    read( deconv_data_y_file_id , *)   !! NUM LEVELS AT WHICH EQ WAS DECONVOLUTED:
                    read( deconv_data_y_file_id , *)   k , tmp_re 
                    
                    if( abs( tmp_re - time_step ) > 1.0e-6 ) then
                        stop 'Error: Inconsistent time steps in Deconv-Y-data and input:time domain analysis. '
                    end if
                    
                    if( drm_layer_info%num_levels .eq. 0 ) then ! EQ was not in X -dir
                        drm_layer_info%num_levels = k + 1 ;
                        
                        if( drm_layer_info%num_levels .eq. 1 ) then
                            stop 'Error: Num levels in DRM-X-Deconv_Data should be >0 time domain analysis.'
                        end if
                        allocate( drm_layer_info%levels( drm_layer_info%num_levels ) , stat = k )
                        if ( k /= 0) stop 'Error: in allocating memory for DRM levels in time domain analysis. '
                    
                        drm_layer_info%levels(1:drm_layer_info%num_levels)  = 0 ! zero value is free surface
                        read( deconv_data_y_file_id , *)   ! LEVELS AT WHICH EQ WAS DECONVOLUTED:
                        read( deconv_data_y_file_id , *) drm_layer_info%levels(2:drm_layer_info%num_levels)
                        read( deconv_data_y_file_id , *)   ! NUM ROWS AND NUM COLUMNS:
                        read( deconv_data_y_file_id , *) drm_layer_info%num_rows , drm_layer_info%num_cols
                        read( deconv_data_y_file_id , *)   ! TIME    A-L1   V-L1   U-L1    A-L2   V-L2    U-L2
                        max_step = min( max_step , drm_layer_info%num_rows )  !! there is no more data in DeconvY-file
                    else
                        read( deconv_data_y_file_id , *)   ! LEVELS AT WHICH EQ WAS DECONVOLUTED:
                        read( deconv_data_y_file_id , *) 
                        read( deconv_data_y_file_id , *)   ! NUM ROWS AND NUM COLUMNS:
                        read( deconv_data_y_file_id , *)   r , s
                        read( deconv_data_y_file_id , *)   ! TIME    A-L1   V-L1   U-L1    A-L2   V-L2    U-L2
                        max_step = min( max_step , r )  !! there is no more data in DeconvY-file
                    end if 
                end if
            end if  
        end if    
    end if    
    
        
        !! fem.str    unit = 21 ( only for dam foundation element at given times)
        !! fem.ten    unit = 19
        !! fem.cmp    unit = 20
        !! output.txt unit = 2
        !! FEM.DIS    unit = 24
       
        
    if( analysis_kind .eq. 2 ) then
        open( unit = 21, file = './output/FEM.STR', action = 'write', access='stream',  &
            & iostat = k, form = 'formatted')
            
        if ( k /= 0) stop 'Error: cannot write FEM.STR file in Time domain analysis suroutine. '
    
               
        !! create file FEM.DIS
        open( unit = 24 , file = './output/FEM.DIS', action = 'write', access='stream', &
            & iostat = k, form = 'formatted')
            
        if ( k /= 0) stop 'Error: cannot write FEM.HIS file in Time domain analysis suroutine. '
    
    
        open( unit = 12, file = './output/FEM.HIS', action = 'write', access='stream', &
            & iostat = k, form = 'formatted')
            
        if ( k /= 0) stop 'Error: cannot write FEM.STR file in Time domain analysis suroutine. '
    
        write(12,'(A10,9A14)') 'TIME','U-NJTF','V-NJTF','U-NRTF','V-NRTF','Uddot-NJTF','Vddot-NJTF' , &
                             & 'Uddot-NRTF','Vddot-NRTF'
    end if 
    
    
                                              
        !! step 1. initialize real solver for Newmark method : 
        !!         tension and compression matrices are also allocated   
            
    call  initialize_real_solver()
    
    !! step 2. Initial values
    if( analysis_kind .eq. 2) then    
        
        tension(:,:)     = 0.0_r_kind ;  !! zero matrix
        compression(:,:) = 0.0_r_kind ;  !! zero matrix
        tension(:,7)     = -10.0e20   ;  !! very small value
        compression(:,9) =  10.0e20   ;  !! very large values
    end if     
         
        !! step 3. Do time domain analysis
        u_vec(:)      = 0.0_r_kind ;
        udot_vec(:)   = 0.0_r_kind ;
        uddot_vec(:)  = 0.0_r_kind ;
        static_solution(:) = 0.0_r_kind ;
     
       
        !!==============================================================
        write(*,*) '       istep , max solver error ' 
        
        do  istep = 1 ,   max_step + num_static_step
            
            !!--------------------------------------------------------------------!!
            !! step 3.1. form equivalent stiffness matrix: K_equiv = K + alpha0 * M + alpha1 * C
            
            
            matrix = stiffness  ;                                   !! 3.1.1. For statis step K_equiv = K 
                
            if( istep > num_static_step ) then                      !! 3.1.2. for dynamic step, other terms come 2 play
                matrix = matrix + alpha0 * mass    ;
                matrix = matrix + alpha1 * damping ;
            end if

            
            
            !! step 3.2. add 7 types of forces: 
            response(:) = 0.0_r_kind ;
            
            !! 3.2.1. body force: dead weigth: when its coeff is not zero !! evaluated in elastic_bulk module for NDYN = 2
            
            tmp_re = load_cases%load_coeffs(1) * load_cases%time_load( istep , 2 ) ;
            if( abs(tmp_re) /= 0.0_r_kind ) then
                response = response + tmp_re * dead_weight_force ;  
            end if
            
            
            !! 3.2.2. concentrainted nodal load = point load
            
            tmp_re = load_cases%load_coeffs(2) * load_cases%time_load(istep,3) ;
            
            if( (nodal_loads%num_point_load /= 0 ) .and. (abs(tmp_re) /= 0.0_r_kind ) ) then
                do  i = 1 , nodal_loads%num_point_load
                    if( degree_of_freedom(nodal_loads%node_ids(i) , 1 ) .eq. 0 ) cycle ! node has no dof
                    do  d = 1 , num_dim
                        k = degree_of_freedom( nodal_loads%node_ids(i) , 2*d )
                        if( k .eq. 0 ) cycle
                        response(k) = response(k) + tmp_re * nodal_loads%load_vector(i,d) ;
                    end do
                end do
            end if
           

            
            !! 3.2.3. surface pressure load !! calculated in elastic bulk module for NDYN = 2
            tmp_re = load_cases%load_coeffs(3) * load_cases%time_load(istep,4) ;
            if( abs( tmp_re ) /= 0.0_r_kind ) then
                response = response + tmp_re * surface_pressure_force ;
            end if
            
            
            !! 3.2.4. enforced displacement
            tmp_re = load_cases%load_coeffs(4) * load_cases%time_load(istep,5) ;
            
            if( enforced_displacement%num_nodes /= 0) then  
                do i = 1 , enforced_displacement%num_nodes
                    if( degree_of_freedom(enforced_displacement%node_ids(i), 1) .eq. 0 ) cycle ! node has no dof
                    do  d = 1 , num_dim
                        if( enforced_displacement%disp_vector(i,d) .eq. 0.0_r_kind ) cycle
                        j = degree_of_freedom( enforced_displacement%node_ids(i) , 2*d )
                        if( j .eq. 0 ) cycle
                        response(j) = response(j) + tmp_re * enforced_displacement%disp_vector(i,d) * large_stiff
                        
                        s = 0 ;
                        call get_variable_location(j , j , s , .false. )
                        matrix(s) = matrix(s) + large_stiff ;
                    end do
                end do
            
            end if
            
            
            !! 3.2.5-7. EQ IN different directions: couped with M_bar_J: not included in DRM analysis
            !! For DRM method force vector only contains effect of reservior-foundation term
            !! and Mbar * J for dam is not included.
            
            
              
            if( istep > num_static_step ) then 
                do  i = 1 , num_dim
                    tmp_re = load_cases%load_coeffs(4+i) * load_cases%time_load(istep,5+i) ;
                    if( abs(tmp_re) /= 0.0_r_kind ) then
                        response(:) = response(:) + tmp_re * force(:,i) ;
                    end if  
                end do
            
                if( drm_layer_ndyn .eq. 1 ) then
                    call  force_in_drm_layer(istep, ffa_file_id , 0 , alpha_m , alpha_k ) 
                else if( drm_layer_ndyn .eq. 2 ) then
                    call  force_in_drm_layer(istep, deconv_data_x_file_id , deconv_data_y_file_id , alpha_m,alpha_k) 
                end if
            end if
            

           
            
            !! step 3.2.8  spring support
            !! k_type =0 (all steps), k_type = 1 (only static steps) , k_type =2 ( only dynamic steps)
            
            do  i = 1 , boundary_springs%num_fixed_nodes
            
                if( degree_of_freedom( boundary_springs%fixed_node_ids(i) , 1 ) .eq. 0) cycle ; !! node has no dof
                
                if( boundary_springs%k_type(i) .eq. 0 ) then
                
                    do  d = 1 , num_dim
                        n = degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d ) ;
                        if( n .eq. 0) cycle ;                           !! node has no dof in d-th direction
                        do k = n , degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d + 1) 
                            if( boundary_springs%disp_vector(i,d) /= 0.0_r_kind ) then
                        
                                response(k) = response(k) + boundary_springs%disp_vector(i,d) *  &
                                            & boundary_springs%spring_stiffness(i,d) ;
                            end if 
                            
                            s = 0 ;
                            call get_variable_location( k , k , s , .false. )
                            matrix(s) = matrix(s) + boundary_springs%spring_stiffness(i,d) ;
                        end do ! over k
                    end do
                    
                else if( boundary_springs%k_type(i) .eq. 1 ) then
                    if((analysis_kind /= 1) .and. (analysis_kind /= 3) .and. (istep > num_static_step)) cycle
                    do  d = 1 , num_dim
                        n = degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d ) ;
                        if( n .eq. 0) cycle ;                           !! node has no dof in d-th direction
                        do k = n , degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d + 1) 
                            if( boundary_springs%disp_vector(i,d) /= 0.0_r_kind ) then
                        
                                response(k) = response(k) + boundary_springs%disp_vector(i,d) *  &
                                            & boundary_springs%spring_stiffness(i,d) ;
                            end if 
                            
                            s = 0 ;
                            call get_variable_location( k , k , s , .false. )
                            matrix(s) = matrix(s) + boundary_springs%spring_stiffness(i,d) ;
                        end do ! over k
                    end do
                    
                else if( boundary_springs%k_type(i) .eq. 2 ) then
                
                    if((analysis_kind /= 2) .and. (analysis_kind /= 4) .and. (istep <= num_static_step)) cycle
                    do  d = 1 , num_dim
                        n = degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d ) ;
                        if( n .eq. 0) cycle ;                           !! node has no dof in d-th direction
                        do k = n , degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d + 1) 
                            if( boundary_springs%disp_vector(i,d) /= 0.0_r_kind ) then
                        
                                response(k) = response(k) + boundary_springs%disp_vector(i,d) *  &
                                            & boundary_springs%spring_stiffness(i,d) ;
                            end if 
                            
                            s = 0 ;
                            call get_variable_location( k , k , s , .false. )
                            matrix(s) = matrix(s) + boundary_springs%spring_stiffness(i,d) ;
                        end do ! over k
                    end do
                    
                end if
            end do
            
            
            !! setting some springs in reservoir truncation boundary for HW model
            if( (the_model_description%reservior_inf_boundary_model_id > 0 ) .and.  &
            &   ( istep .le. num_static_step ) ) then  
                
                !! loop over all nodes for finding nodes at reservoir truncation boundary
                do i = 1 , num_node
                    if( degree_of_freedom(i,1) /= 3 ) cycle             !! nodes at reservior truncation boundary
                    do  j = num_dim + 1 , num_degree_of_freedom
                        if( degree_of_freedom(i,2*j) .eq. 0 ) cycle ;
                        do  k = degree_of_freedom(i,2*j) , degree_of_freedom(i,2*j+1)
                            s = 0 ;
                            call get_variable_location( k , k , s , .false. )
                            matrix(s) = matrix(s) + 1.0e20 ;
                        end do
                    end do
                    
                end do
            end if
            
            
            
            !! 3.2.9. add damping inertia to response vector for dynamic steps
            if( istep > num_static_step ) then
                
                !! 3.2.9.1. form vector T = a0 * U + a2 * Udot + a3 * Uddot
                temp_real_vec(:) = alpha0 * u_vec(:) + alpha2 * udot_vec(:) + alpha3 * uddot_vec(:) ;
                
                !! 3.2.9.2. Calculate Response = Response +  M * T
                do  i = 1 , size( dictionary ) -1
                    tmp_re = 0.0_r_kind ;
                    do j = dictionary(i) + 1 , dictionary(i+1)
                        tmp_re = tmp_re + mass(j) * temp_real_vec( dependency(j) ) ;
                    end do
                    response(i) = response(i) + tmp_re ;
                end do ! over i
                
                
                !! 3.2.9.3. form vector T = a1 * U + a4 * Udot + a5 * Uddot
                temp_real_vec(:) = alpha1 * u_vec(:) + alpha4 * udot_vec(:) + alpha5 * uddot_vec(:) ;
                
                !! 3.2.9.4. Calculate Response = Response +  C * T
                
                do  i = 1 , size( dictionary ) -1
                    tmp_re = 0.0_r_kind ;
                    do j = dictionary(i) + 1 , dictionary(i+1)
                        tmp_re = tmp_re + damping(j) * temp_real_vec( dependency(j) ) ;
                    end do
                    response(i) = response(i) + tmp_re ;
                end do ! over i
                
                
                !! subtraction of spring force from static solution
                do  i = 1 , boundary_springs%num_fixed_nodes
            
                    if( degree_of_freedom( boundary_springs%fixed_node_ids(i) , 1 ) .eq. 0) cycle ; !! node has no dof
                
                    if( boundary_springs%k_type(i) /=1 ) cycle
                
                    do  d = 1 , num_dim
                        n = degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d ) ;
                        if( n .eq. 0) cycle ;                           !! node has no dof in d-th direction
                        do k = n , degree_of_freedom( boundary_springs%fixed_node_ids(i) , 2*d + 1) 
                    response(k) = response(k) - static_solution(k) * boundary_springs%spring_stiffness(i,d) ;
                        end do ! over k
                    end do ! over d
                end do ! over i
                
                
                !! subtraction of spring force in reservoir truncation boundary
                if( the_model_description%reservior_inf_boundary_model_id > 0 ) then  
                
                    !! loop over all nodes for finding nodes at reservoir truncation boundary
                    do  i = 1 , num_node
                        if( degree_of_freedom(i,1) /= 3 ) cycle             !! nodes at reservior truncation boundary
                        do  j = num_dim + 1 , num_degree_of_freedom
                            if( degree_of_freedom(i,2*j) .eq. 0 ) cycle ;
                            do  k = degree_of_freedom(i,2*j) , degree_of_freedom(i,2*j+1)
                                response(k) = response(k) - static_solution(k) * 1.0e20 ;
                            end do
                        end do
                    
                    end do
                end if
            
                
            end if
            
            
            !! 3.2.10. Solve system
            call  solve_real_system(istep) ;
            
            
            if( analysis_kind .eq. 2 ) then
                !! print data for njtf node
                
                if( istep > num_static_step ) then
                    point_data_vector(:) = 0.0_r_kind ;
                    do  j = 1 , 2 !num_dim
                        n = degree_of_freedom(report_node_index , 2*j) ;
                        if( n /= 0 ) then
                            point_data_vector(j)    = u_vec(n)     ! u
                            point_data_vector(j+ 4) = uddot_vec(n) ! uddot
                        end if 
                    end do
                    
                    do  j = 1 , 2 !num_dim
                        n = degree_of_freedom(reference_node_index , 2*j) ;
                        if( n /= 0 ) then
                            point_data_vector(j+2)   = u_vec(n)     ! u
                            point_data_vector(j+ 6)  = uddot_vec(n) ! uddot
                        end if 
                    end do
                    write(12,'(E15.6,8E14.6)') load_cases%time_load(istep,1),point_data_vector(1:8)
                end if
            
                !! 3.2.11. Strain-stress analysis
                !! 3.2.11.1. REASSEMBLE SUBROUTINE: print some data
            
                is_data_printed = .false.
            
                tmp_re = real( istep - num_static_step - 1, kind = r_kind) * time_step  ;
            
                do  i = 1 , num_node_report_time_data
                    if( abs( tmp_re - reported_time(i) ) < epsilon( tmp_re) ) then
                        is_data_printed = .true.
                        write(24,'(A21,F8.2,A12,I5,A1)') 'DISPLACEMENTS AT TIME', tmp_re ,'  :  (ISTEP=',istep,')'
                        write(24,'(I10)') num_node
                    end if 
                end do ! over i
            
                if( is_data_printed ) then
                    do i = num_node , 1 , -1
                
                        point_data_vector(:) = 0.0_r_kind ;
                        do j = 1 , num_dim
                            n = degree_of_freedom(i , 2*j) ;
                            if( n /= 0 ) then
                                point_data_vector(j) = solution(n)
                            end if 
                        end do
                    
                        if( num_degree_of_freedom > num_dim ) then   !! there is pressure dof
                            n = degree_of_freedom(i , 2 * num_degree_of_freedom ) ;
                            if( n /= 0 ) then
                                point_data_vector(4) = solution(n) ;
                            end if
                        end if
                    
                    write(24,'(I6,4X,7ES14.6)') i,point_data_vector(1:3),0.0D0,0.0D0,0.0D0,point_data_vector(4)
                    end do
                end if
            
                !! end of REASSEMBLE SUBROUTINE
           
                call  find_strain_stress( istep , tmp_re , is_data_printed , analysis_kind, num_region_element ) ;
            end if 
            
            !! for freq domain in time domain analysis
            if( (analysis_kind .eq. 10) .and. ( istep > num_static_step) ) then
                do  j = 1 , num_dim
                    k = degree_of_freedom( report_node_index , 2 * j )  !! NJTF node
                    if( k .eq. 0 ) cycle                                ! node has no translational dof in d-th direction
                    u_uddot_njtf_drm(istep - num_static_step , j)   = solution(k)  ;      !! insert u vector
                    u_uddot_njtf_drm(istep - num_static_step , num_dim +j) = uddot_vec(k) ;      !! insert u_ddot vector
                end do
            end if 
            
            
                !! 3.2.12. update velocity and acceleration and displacement
            if(  istep > num_static_step ) then
                
                u_vec = solution - u_vec ;    !! Ubar = U_{n+1} - U_n
            
                !! store Udot_{n+1} in temp_vec = Udot_{n+1} = alpha1 * Ubar - alpha2 * Udot_n - alpha3 * Uddot_n 
                temp_real_vec = alpha1 * u_vec ;   
                temp_real_vec = temp_real_vec - alpha4 * udot_vec ;
                temp_real_vec = temp_real_vec - alpha5 * uddot_vec ;
            
                !! update Uddot_{n+1}  =  alpha1 * Ubar - alpha4 * Udot_n - alpha5 * Uddot_n
                uddot_vec = -alpha3 * uddot_vec ;  
                uddot_vec = uddot_vec + alpha0 * u_vec ;
                uddot_vec = uddot_vec - alpha2 * udot_vec ;
            
                !! update Udot_{n+1}
                udot_vec = temp_real_vec  ;
                
                !! update U_vec
                u_vec = solution ;
                
            end if
            
            
            
            !! store static step solution
            if( istep .eq. num_static_step ) then
                static_solution = solution ;
                !! update U_vec
                u_vec = solution ;
            end if
            
            
            !! if ndyn = 2 then tension and compression are satisfied at istep = 1310
      
        end do ! over istep
      
      
        if( analysis_kind .eq. 2) then
            close(12)  !! close file FEM.HIS
            close(21)  !! close file FEM.STR
            close(24)  !! close file FEM.DIS
        end if 
        
        if( ffa_file_id > 0           ) close( ffa_file_id          )  !! close file FreeFieldData
        if( deconv_data_x_file_id > 0 ) close( deconv_data_x_file_id)  !! close file for DRM-DeconvX-data
        if( deconv_data_y_file_id > 0 ) close( deconv_data_y_file_id)  !! close file for DRM-DeconvY-data
           
        
       
 
        
        !! write tension and compression    
        open(unit = 19, file = './output/FEM.TEN', action = 'write', access='stream', iostat = k, form = 'formatted')
        if ( k /= 0) stop 'Error: cannot write FEM.TEN file in Time domain analysis suroutine. '
        open(unit = 20, file = './output/FEM.CMP', action = 'write', access='stream', iostat = k, form = 'formatted')
        if ( k /= 0) stop 'Error: cannot write FEM.CMP file in Time domain analysis suroutine. '
    
        
        tmp_re =  real( max_step - num_static_step - 1, kind = r_kind) * time_step  ;
        
        d = 0 ;                                                         !! Igauss
        
        
        
        
        write(19,'(A18,F7.2,A14)') 'SIG-MAX UP TO TIME',tmp_re ,'   :   (*.TEN)'
        write(20,'(A18,F7.2,A14)') 'SIG-MIN UP TO TIME',tmp_re ,'   :   (*.CMP)'
        tmp_re = 0.0_r_kind ;
        do  i = 1 , size( element_property , 1)
            
            !! get local element id
            istep = element_property(i , 1 )                            !! igroup
            if( istep .eq. 1 ) then
                j = i
            elseif( istep .eq. 2 ) then
                j = i - num_region_element(1)
            else
                j = i - num_region_element(1) - num_region_element(2)
            end if
            
            num_var = element_property(i , 3 )                           !! n_gauss_point
            call  set_gauss_points( num_var ) ;
            
            
            
            n = 0                                                       !! num gauss point
            
            if( element_property(i,1) < 3 ) then                        !! elastic element
            
                if( num_dim .eq. 2 ) then
                    do k = 1 , num_var
                        do  r = 1 , num_var
                            n = n + 1;
                            d = d + 1;
                            
                            if( n .eq. 1 ) then
                                write(19,'(I12,I10)') num_var**2,1
                                write(20,'(I12,I10)') num_var**2,1
                            end if
                            
                            write(19,'(3I4,I2,1X,20E11.4)') j,istep,n,1,tension(d,1:9),tmp_re,tmp_re,  &
                            & gauss_points(r),gauss_points(k),tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            
                            write(20,'(3I4,I2,1X,20E11.4)') j,istep,n,1,compression(d,1:9),tmp_re,tmp_re,  &
                            & gauss_points(r),gauss_points(k),tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                        end do ! over r
                    end do !over k
                else                         !! dim = 3
                    do k = 1 , num_var
                        do  r = 1 , num_var
                            do  s = 1 , num_var
                            n = n + 1;
                            d = d + 1;
                    write(19,'(3I4,I2,1X,20E11.4)') j,istep,n,1,tension(d,1:9),tmp_re,tmp_re,gauss_points(k),&
                    & gauss_points(r),gauss_points(s),tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            
                    write(20,'(3I4,I2,1X,20E11.4)') j,istep,n,1,compression(d,1:9),tmp_re,tmp_re,gauss_points(k),&
                    & gauss_points(r),gauss_points(s),tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            end do ! over s
                        end do ! over r
                    end do !over k
                
                end if
            
            else  !! fluid elements
            
                if( num_dim .eq. 2 ) then
                    do k = 1 , num_var
                        do  r = 1 , num_var
                            n = n + 1;
                            d = d + 1;
                            if( n .eq. 1 ) then
                                write(19,'(I12,I10)') num_var**2,1
                                write(20,'(I12,I10)') num_var**2,1
                            end if
                            write(19,'(3I4,I2,1X,20E11.4)') j,istep,n,1,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, &
                            & tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, gauss_points(k),gauss_points(r),&
                            & tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            
                            write(20,'(3I4,I2,1X,20E11.4)') j,istep,n,1,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, &
                            & tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, gauss_points(k),gauss_points(r),&
                            & tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            
                        end do ! over r
                    end do !over k
                else                         !! dim = 3
                    do k = 1 , num_var
                        do  r = 1 , num_var
                            do  s = 1 , num_var
                            n = n + 1;
                            d = d + 1;
                            
                            if( n .eq. 1 ) then
                                write(19,'(I12,I10)') num_var**2,1
                                write(20,'(I12,I10)') num_var**2,1
                            end if
                            
                    write(19,'(3I4,I2,1X,20E11.4)') j,istep,n,1,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, &
                            & tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, gauss_points(k),gauss_points(r),&
                            & gauss_points(s),tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            
                    write(20,'(3I4,I2,1X,20E11.4)') j,istep,n,1,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, &
                            & tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re, gauss_points(k),gauss_points(r),&
                            & gauss_points(s),tmp_re,tmp_re,tmp_re,tmp_re,tmp_re,tmp_re
                            
                            end do ! over s
                        end do ! over r
                    end do !over k      
                end if 
            end if
         
        end do ! over i
    
        close(19)  !! close file FEM.TEN
        close(20)  !! close file FEM.CMP
        
    
    end subroutine time_domain_analysis
    
    
!!====================================================================!!
!!
!!
!!
!!====================================================================!!
    
    subroutine fft_transformation()
        use  mod_physics  , only : u_uddot_njtf_drm
        use  mod_geometry , only : num_dim 
        implicit none
        
        integer :: i,j  , n ,d , num_omega 
        real( kind = r_kind) :: freq , theta,  pi = acos(-1.0_r_kind) , max_omega  = 100.0_r_kind  ;  
        complex( kind = r_kind) :: ux_cplx , uy_cplx , ar_cplx , hx_cplx, hy_cplx  !! u_x, u_y , ddot_u, Hx, Hy, Hz
        
        n = size( u_uddot_njtf_drm , 1)
        
        !! find radial accelration
        u_uddot_njtf_drm(:,5) = sqrt( u_uddot_njtf_drm(:,3) ** 2 + u_uddot_njtf_drm(:,4) ** 2) ;
       
        
        num_omega = ceiling( max_omega / frequency_step ) ;
        
        ! Open input file 
        open ( unit   = 99  , file    = "./output/OUTPUT.TXT"  , action  = 'write'   ,  &
            access = 'stream' , iostat  = i             , form    = 'formatted' )
                                  
        if ( i /= 0) stop 'Error: cannot write output.txt file. '
     
        write(99,'(A)') 'TRANSFER FUNCTION DETERMINATION FOR ARBITRARY FREQUENCIES:'
!        write(99,'(A29,F10.5)') 'VALUE OF DAMPING COEFFICIENT:', hysteresis_damping_coeff
     
!        write(99,'(A40,I6,A32,I6,A32,I6)') 'ACCEL. TRANS. FUNC. FOR THE JOINT NUM.', report_node_index , &
!                 'RELATIVE TO THE JOINT NUMABER', reference_node_index ,'BY GROUND MOTION IN DIRECTION',&
!                 earthquake_direction
        write(99,'(A19,6A50)') 'FREQUENCY VALUE','ACCEL. TRANS. FUNC. IN DIRECTION 1',&
            &  'ACCEL. TRANS. FUNC. IN DIRECTION 2','ACCEL. TRANS. FUNC. IN DIRECTION 3',&
            &  'ABSOLUTE VALUE IN DIRECTION 1','ABSOLUTE VALUE IN DIRECTION 2',&
            & 'ABSOLUTE VALUE IN DIRECTION 3'
        
        
        
        !! do freq domain integral
        do  i = 1 , num_omega
            freq   = real( i - 1 , kind = r_kind ) * frequency_step ;
            
            !! take integral by trapezoid role
            theta  = real( n - 1 , kind = r_kind ) * time_step * pi * freq / 180.0_r_kind     ;
            ux_cplx = (u_uddot_njtf_drm(1,1) + u_uddot_njtf_drm(n,1) * cmplx( cos(theta) , sin(theta)) ) /2.0_r_kind ;
            uy_cplx = (u_uddot_njtf_drm(1,2) + u_uddot_njtf_drm(n,2) * cmplx( cos(theta) , sin(theta)) ) /2.0_r_kind ;
            ar_cplx = (u_uddot_njtf_drm(1,5) + u_uddot_njtf_drm(n,5) * cmplx( cos(theta) , sin(theta)) ) /2.0_r_kind ;
            
            do  j = 2  , n-1
                theta  = real( j - 1 , kind = r_kind ) * time_step * pi * freq / 180.0_r_kind     ;
                ux_cplx  = ux_cplx  + u_uddot_njtf_drm(j,1) * cmplx( cos(theta) , sin(theta)) ;
                uy_cplx  = uy_cplx  + u_uddot_njtf_drm(j,2) * cmplx( cos(theta) , sin(theta)) ;
                ar_cplx  = ar_cplx  + u_uddot_njtf_drm(j,5) * cmplx( cos(theta) , sin(theta)) ;
            end do 
            
            hx_cplx = ux_cplx / ar_cplx ;
            hy_cplx = uy_cplx / ar_cplx ;
            
!            write(99,*) freq , real(hx_cplx), aimag(hx_cplx),'i', &
!            &   real(hy_cplx), aimag(hy_cplx),'i',&
!            &   0.0_r_kind, 0.0_r_kind,'i', abs(hx_cplx ) , abs(hy_cplx ) , 0.0_r_kind 
                  
            write(99,*) freq , abs(hx_cplx ) , abs(hy_cplx ) , 0.0_r_kind 

        end do 
        
        
!        do  i = 1 , n   
            
!            write(99,*) u_uddot_njtf_drm(i,:)
!!            freq_data(1:9) = cmplx( 0.0_r_kind , 0.0_r_kind)  
            
!!            do  d = 1 , 2*num_dim  !! ndim for u and ndim for ddot_u
!!                tmp_cplx = cmplx( 0.0_r_kind , 0.0_r_kind )
                
!!                !! fourier transform
!!                do  j = 1 , n
!!                    theta = real( 2 *(j-1)*i , kind = r_kind) * freq0 !/ real( n , kind = r_kind) 
!!                    tmp_cplx%re = tmp_cplx%re + u_uddot_njtf_drm(j,d  ) * cos( pi * theta ) ;
!!                    tmp_cplx%im = tmp_cplx%im - u_uddot_njtf_drm(j,d+2) * sin( pi * theta ) ;
!!                end do 
                
!!                if( d <= num_dim ) then
!!                    freq_data(d) = tmp_cplx ;
!!                else
!!                    freq_data(d-num_dim+3) = tmp_cplx ;
!!                end if
!!            end do
            
!!            !! calculate H 
!!            do  d = 1 , num_dim
!!                if( abs( freq_data( d+3) ) > 1.0e-14) then
!!                    freq_data(6 + d) = freq_data(d) / freq_data(d+3)
!!                end if
!!            end do 
            
!!            freq = real( i-1 , kind = r_kind) * freq0 ;
            
!!!            write(99,*) i , real( freq_data(7)), aimag( freq_data(7) ),'i', &
!!!            &   real(freq_data(8)), aimag(freq_data(8)),'i',&
!!!            &   real(freq_data(9)), aimag(freq_data(9)),'i', abs(freq_data(7:9) )
                  
!!            write(99,*) freq , abs(freq_data(7:9) )
                 
            
!        end do 
        
        
        close(99)   
        
    end subroutine fft_transformation 
  
!!====================================================================!!
!!
!!  free unused memory before doing main analysis
!!
!!====================================================================!!
    
    subroutine  free_unused_memory()
      
        use  mod_geometry , only : node , element_matrix
        use  mod_physics  , only : material , boundary_load , element_property
      
        implicit none
        integer :: state , state1 = 0
      
        if( allocated( node )) then
            deallocate( node , stat = state)
            state1 = state1 + state
        end if
       
        if( allocated( element_matrix ) ) then
            deallocate( element_matrix , stat = state )
            state1 = state1 + state
        end if
       
        if( allocated( material ) ) then
            deallocate( material , stat = state )
            state1 = state1 + state
        end if
       
        if( allocated( boundary_load ) ) then
            deallocate( boundary_load , stat = state )
            state1 = state1 + state
        end if
       
       
        if( allocated( element_property ) ) then
            deallocate( element_property , stat = state )
            state1 = state1 + state
        end if
       
        ! check if things are done correctly
        if( state1 /= 0 ) stop 'Error: failed in deallocating unused arrays in analysis module.'  
                  
    
    end subroutine  free_unused_memory
 
  
!!====================================================================!!
!!
!!  frequency domain analysis
!!
!!====================================================================!!
   
    subroutine  freq_domain_analysis()
      
        use  mod_physics  , only : mass , stiffness , damping , force , dictionary , dependency 
        use  mod_physics  , only : earthquake_direction , translational_variable, drm_layer_info
        use  mod_physics  , only : find_DRM_internal_face
        use  mod_geometry , only : num_dim , degree_of_freedom , node
        use  mod_solver   , only : initialize_complex_solver , clean_up_solver , solve_complex_system
        use  mod_solver   , only : cplx_matrix , cplx_response , cplx_solution
        use  mod_fem      , only : freq_domain_force_in_drm_layer
        implicit none
      
        integer :: i , k , num_steps , d  , j
        real( kind = r_kind ) :: a_g = 1.0_r_kind                         ! magnitude of earthquake in direction: ngdir
        real( kind = r_kind ) :: omega = 0.0_r_kind , final_omega = 100.0_r_kind       ! initial and final value of omega.
        real( kind = r_kind ) ::  tmp_re , tmp_im 
        complex( kind = r_kind ) :: relative_displacement(3) , cpl_test(9)
        integer :: evaluation_nodes(2)
        logical :: found ;
        
        
        
        !! specify nodes with translational dof
        allocate( translational_variable( size( dictionary ) - 1 ) , stat = k )
        if( k /=0 ) stop 'Error in allocating variable class array in physics module '
     
      
     
        translational_variable(:) = 0                                   !! non-zero value for vars with translational dof
     
        do  i = 1 , size( degree_of_freedom , 1 )                          !! over nodes 
            if( degree_of_freedom( i , 1 ) .eq. 0 ) cycle
            do  j = 1 , num_dim
                if( degree_of_freedom( i , 2 * j ) .eq. 0 ) cycle
                translational_variable( degree_of_freedom( i , 2 * j ) ) = 1
            end do ! over j
        end do ! over i
     
      
        evaluation_nodes(1) = report_node_index         ! NJTF
        evaluation_nodes(2) = reference_node_index      ! NRTF
     

        num_steps = floor( ( final_omega + 0.010_r_kind * frequency_step ) / frequency_step )
      
       
        call  initialize_complex_solver()

    
        if( drm_layer_info%layer_id > 0 ) then
            call   find_DRM_internal_face()
        end if
    
      
        ! Open input file 
        open ( unit   = 101  , file    = "./output/OUTPUT.TXT"  , action  = 'write'   ,  &
            access = 'stream' , iostat  = k             , form    = 'formatted' )
                                  
        if ( k /= 0) stop 'Error: cannot write output.txt file. '
     
        write(101,'(A)') 'TRANSFER FUNCTION DETERMINATION FOR ARBITRARY FREQUENCIES:'
        write(101,'(A29,F10.5)') 'VALUE OF DAMPING COEFFICIENT:', hysteresis_damping_coeff
     
        write(101,'(A40,I6,A32,I6,A32,I6)') 'ACCEL. TRANS. FUNC. FOR THE JOINT NUM.', report_node_index , &
                 'RELATIVE TO THE JOINT NUMABER', reference_node_index ,'BY GROUND MOTION IN DIRECTION',&
                 earthquake_direction
        write(101,'(A19,6A50)') 'FREQUENCY VALUE','ACCEL. TRANS. FUNC. IN DIRECTION 1',&
            &  'ACCEL. TRANS. FUNC. IN DIRECTION 2','ACCEL. TRANS. FUNC. IN DIRECTION 3',&
            &  'ABSOLUTE VALUE IN DIRECTION 1','ABSOLUTE VALUE IN DIRECTION 2',&
            & 'ABSOLUTE VALUE IN DIRECTION 3'
           
      
      
        do  i = 1 ,   num_steps
      
            ! step 1. increment frequency
          
            omega = omega + frequency_step    
          
!!--------------------------------------------------------------------!!
        !! step 2. form coefficient matrix
           
            cplx_matrix  = cmplx( 0.0_r_kind , 0.0_r_kind )
           
            do  k = 1 , size( mass )
                tmp_re = stiffness(k) - ( omega ** 2 ) * mass(k)
                tmp_im = omega *  damping(k) 
                cplx_matrix(k) = cmplx( tmp_re , tmp_im )
            end do
          
          
            !! step 2.3. add hysteresis damping : C_elastic = 2 \beta * K_elastic  to all variables describing translational dof
          
            do  k = 1 , size( translational_variable )
          
                if( translational_variable( k ) .eq. 0 ) cycle    ! the variable = eq has no translational dof
                do  j = dictionary( k ) + 1 , dictionary( k + 1 )          !! take k-th equation
                    if ( translational_variable( dependency(j) ) .eq. 0 ) cycle
                    tmp_im = 2.0_r_kind * hysteresis_damping_coeff * stiffness(j)
                    cplx_matrix(j) = cplx_matrix(j) + cmplx( 0.0_r_kind , tmp_im )
                end do ! over j
          
            end do ! over k
          

!!!--------------------------------------------------------------------!! 
    !!! step 3. form response vector
           
            cplx_response = cmplx( 0.0_r_kind , 0.0_r_kind )
            do k = 1 , size( force , 1 )
                tmp_re = a_g * force( k , earthquake_direction )
                cplx_response(k) =  cmplx( tmp_re , 0.0_r_kind )
            end do

!!!-------------------------------------------------------------------!!
    !! step 3.1. force in DRM layer
    if( drm_layer_info%layer_id > 0 ) then
        call  freq_domain_force_in_drm_layer( omega , a_g , hysteresis_damping_coeff , earthquake_direction )
    end if 
            
!!!--------------------------------------------------------------------!!
    !!! step 4. solve equations locally
          
            call  solve_complex_system( omega )
            
            
!!--------------------------------------------------------------------!!
!! step 5. evaluate relative acceleration
          
            relative_displacement = cmplx( 0.0_r_kind , 0.0_r_kind )
            do  d = 1 , num_dim
                k = degree_of_freedom( evaluation_nodes(1) , 2 * d )
                if( k .eq. 0 ) cycle                                    ! node has no translational dof in d-th direction
                relative_displacement(d) = cplx_solution( k )
            end do
            
            
            do  d = 1 , num_dim
                k = degree_of_freedom( evaluation_nodes(2) , 2 * d )
                if( k .eq. 0 ) cycle                                    ! node has no translational dof in d-th direction
                relative_displacement(d) = relative_displacement(d) - cplx_solution( k )
            end do
          
          
            relative_displacement = relative_displacement * ( omega ** 2 ) / a_g
          
!            write(101,*) omega , real( relative_displacement(1)), aimag( relative_displacement(1) ),'i', &
!            &   real(relative_displacement(2)), aimag(relative_displacement(2)),'i',&
!            &   real(relative_displacement(3)), aimag(relative_displacement(3)),'i', abs(relative_displacement(1:3) )
                  
              write(101,*) omega , abs(relative_displacement(1:3) )
                  
   
       
!!--------------------------------------------------------------------!!        
!! step 6. write in output file
          
        
      
        end do ! over i
      
      
        close( 101 )
    
      
        call  clean_up_solver()
      
     
   
    end subroutine freq_domain_analysis
   
   
   
   
  
   
   
    
   
!!====================================================================!!
!!
!!  store data in file
!!
!!====================================================================!!
   
    subroutine  store_in_file()
     
        use mod_physics, only: mass, stiffness, damping, force, dictionary, dependency, num_max_variables_in_eqns
        use mod_geometry, only: degree_of_freedom, num_dim, num_degree_of_freedom
        use mod_physics, only :  translational_variable
    
        implicit none
     
        integer :: i , state , k  , j , num_var , node_id , dof_id , nodal_var_id, r
        integer , allocatable :: non_zero_index(:)
        integer , allocatable :: convert_index(:)
        real( kind = r_kind ) :: epsi


        !! for index convertor and storing stiffness matrix
        !! index of array is var index in my notation, the value is that of hesam
        allocate( convert_index( size(dictionary)-1 ) , stat = state )
        if( state /= 0 ) stop ' failed in allocating non-zero_index in analysis module '
        
        convert_index(:) = 0 ;
        
        !! first iterate over translational dof
        num_var = 0 ;
        do i = 1 , size( degree_of_freedom , 1 )
            do  j = 1 , num_dim
                k = degree_of_freedom(i,2*j) ;
                if( k .eq. 0 ) cycle
                num_var = num_var + 1 ;
                convert_index(k) = num_var ; 
            end do 
        end do
     
        !! then iterate over pressure dof
        if( num_degree_of_freedom > num_dim ) then !! there is pressure
            do i = 1 , size( degree_of_freedom , 1 )
                k = degree_of_freedom(i,2*num_degree_of_freedom ) ;
                if( k .eq. 0 ) cycle
                num_var = num_var + 1 ;
                convert_index(k) = num_var ; 
            end do
        end if 
        
        
        
        
        !! print data in new index set
        open ( unit   = 88  , file    = "./output/FEMData.txt"  , action  = 'write'   ,  &
               access = 'stream' , iostat  = state             , form    = 'formatted' )

        if ( state /= 0) stop 'Error: cannot write output file in analysis module. '
        
        !! write global matrices
        write(88, '(A)') 'Hesam-I, Hesam J, mass(i,j) , damping(i,j),stiffness(i,j)'
        write(88, '(I12)') size(stiffness)
        do i = 1 , size( dictionary ) - 1
            k = convert_index(i) ;
            do  j = dictionary(i)+1 , dictionary(i+1)
                r = convert_index( dependency(j) ) ;
                write(88,'(2I12,3F25.15)') k,r, mass(j) , damping(j) , stiffness(j) 
            end do
        end do 
        
        write(88, '(A)') 'Hesam-I, force(I,J) '
        write(88,'(2I12)')  size(force,1) , size(force , 2)
        if( num_dim .eq. 2 ) then
            do  i = 1 , size( force ,1)
                k = convert_index(i) ;
                write(88,'(I12,2F25.15)') k , force(i,1:2)
            end do 
        else
            do  i = 1 , size( force ,1 )
                k = convert_index(i) ;
                write(88,'(I12,3F25.15)') k , force(i,1:3)
            end do 
        end if
        close(88) ;
        
        deallocate( convert_index , stat = state )
        if ( state /= 0) stop 'failed in deallocating convert_index vector in analysis module. '

         
!        write(*,*) 'degree of freedom '
!        do i = 1 , size( degree_of_freedom , 1 )
!           write(*,*) i, degree_of_freedom(i,:)
!        end do
        
!        write(*,*) 'var index '
!        do i = 1 , size( convert_index )
!           write(*,*) i , convert_index(i)
!        end do
        
        
       
        
!        open ( unit   = 88  , file    = "./output/nnz_column.txt"  , action  = 'write'   ,  &
!            access = 'stream' , iostat  = state             , form    = 'formatted' )
            
            
!        epsi = 1.0E-14
     
     
!        open ( unit   = 88  , file    = "./output/nnz_column.txt"  , action  = 'write'   ,  &
!            access = 'stream' , iostat  = state             , form    = 'formatted' )
     
!        if ( state /= 0) stop 'Error: cannot write output file in analysis module. '
     
     
!        do  i = 1 , size( dictionary ) -1  
!            write(88 , * ) 'row = ' ,  i , ':' 
!            j = dictionary(i) + 1 ;
!            k = dictionary(i+ 1) ;
!            write(88 , * )  dependency(j:k)
!            write(88 , * ) 
!        end do
     
!        close(88)  
    
    
!        open ( unit   = 87  , file    = "./output/nnz_column_beg_end.txt"  , action  = 'write'   ,  &
!               access = 'stream' , iostat  = state             , form    = 'formatted' )
     
!        if ( state /= 0) stop 'Error: cannot write output file in analysis module. '
     
     
!        do  i = 1 , size( dictionary ) -1  
        
!            j = dictionary(i) + 1 ;
!            k = dictionary(i+ 1) ;
!            write(87 , * ) i , dependency(j) , dependency(k)
!        end do
     
!        close(87)  
     
     
!        allocate( non_zero_index( 3 * num_max_variables_in_eqns ) , stat = state )
!        if( state /= 0 ) stop ' failed in allocating non-zero_index in analysis module '
     
!        open ( unit   = 101  , file    = "./output/mass_matrix.txt"  , action  = 'write'   ,  &
!            access = 'stream' , iostat  = state             , form    = 'formatted' )
     
!        if ( state /= 0) stop 'Error: cannot write output file in analysis module. '
     
!        write( 101 , * ) '!! ********************   Explanation of the file  ************************ !!'
!        write( 101 , * ) '!!'
!        write( 101 , * ) '!! To every node (given by node_id) some dofs ( given by dof_id ) is assigned. '
!        write( 101 , * ) '!!        dof_id = 1 : displacement in X direction ( u_x )'
!        write( 101 , * ) '!!        dof_id = 2 : displacement in Y direction ( u_y )'
!        write( 101 , * ) '!!        dof_id = 3 : in (3D) displacement in Z direction ( u_z )'
!        write( 101 , * ) '!!                   : in 2D problem pressure ( P )  '
!        write( 101 , * ) '!!        dof_id = 4 : in 3D problem is pressure (P) '
!        write( 101 , * ) '!!'
!        write( 101 , * ) '!! Every dof is desribed by a set of variables ( given by dof_var_id ) '
!        write( 101 , * ) '!! Values of dof_var_id: '
!        write( 101 , * ) '!!        dof_id = 1 , dof_var_id = 1  means function u_x. '
!        write( 101 , * ) '!!        dof_id = 2 , dof_var_id = 1  means function u_y. '
!        write( 101 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1  means function P. '
!        write( 101 , * ) '!!'
!        write( 101 , * ) '!! For nodes at reservior truncation boundary: HW model of order( N , M ): '
!        write( 101 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2        means function phi_1. '
!        write( 101 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1+j      means function phi_j. '
!        write( 101 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1+N+1    means function psi. '
!        write( 101 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2+N+j    means function phi_{N+j}. '
!        write( 101 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2+N+M+1  means function P^i ( induced pressure). '
!        write( 101 , * ) '!!'
!        write( 101 , * ) '!! Row_id = index of row in global matrix. '
!        write( 101 , * ) '!! col_id = index of col in global matrix that possibly has non-zero value ( up to 1E-13 ).'
!        write( 101 , * ) '!!'
!        write( 101 , * ) '!! ******************************************************************************* !!'
!        write( 101 , * )
!        write( 101 , * )
!        write( 101 , * )
     
     
!        do  i = 1 , size( dictionary ) -1 
!            non_zero_index = 0
!            num_var = 0
!            do  k  = dictionary( i ) + 1 , dictionary( i + 1 )
!                if( abs( mass( k ) ) > epsi ) then
!                    num_var = num_var + 1
!                    non_zero_index(num_var) = k 
!                end if
!            end do
         
!            call  find_variable_identifier( i , node_id , dof_id , nodal_var_id )
!            write( 101 , * ) 
!            write( 101 , * ) ' ***************************** '
!            write( 101 , * ) '  Row id ,    Node_id ,     dof_id ,   dof_var_id   ' 
!            write( 101 , * ) i , node_id , dof_id , nodal_var_id 
          
!            write( 101 , * ) 'Maximum number of non-zero columns in this row:  ' , num_var
         
!            if( num_var .eq. 0 ) cycle 
         
!            write( 101 , * ) 'List of columns with potentially non-zero values: '
!            write( 101 , * ) '        col_id  ,  node_id  ,  dof_id , dof_var_id  ,  mass(row_id,col_id) '
!            do k = 1 , num_var
!            j = dependency( non_zero_index(k) )
!                call find_variable_identifier( j , node_id , dof_id , nodal_var_id )
!                write( 101 , * ) j , node_id , dof_id , nodal_var_id , '    '  , mass( non_zero_index(k) )
!            end do
     
!        end do
     
!        close( 101 )
     
     
!        open ( unit   = 102  , file    = "./output/stiffness_matrix.txt"  , action  = 'write'   ,  &
!               access = 'stream' , iostat  = state             , form    = 'formatted' )
     
!        if ( state /= 0) stop 'Error: cannot write output file in analysis module. '
     
!        write( 102 , * ) '!! **********************  Explanation of the file:  *********************** !!'
!        write( 102 , * ) '!!'
!        write( 102 , * ) '!! To every node (given by node_id) some dofs ( given by dof_id ) is assigned. '
!        write( 102 , * ) '!!        dof_id = 1 : displacement in X direction ( u_x )'
!        write( 102 , * ) '!!        dof_id = 2 : displacement in Y direction ( u_y )'
!        write( 102 , * ) '!!        dof_id = 3 : in (3D) displacement in Z direction ( u_z )'
!        write( 102 , * ) '!!                   : in 2D problem pressure ( P )  '
!        write( 102 , * ) '!!        dof_id = 4 : in 3D problem is pressure (P) '
!        write( 102 , * ) '!!'
!        write( 102 , * ) '!! Every dof is desribed by a set of variables ( given by dof_var_id ) '
!        write( 102 , * ) '!! Values of dof_var_id: '
!        write( 102 , * ) '!!        dof_id = 1 , dof_var_id = 1  means function u_x. '
!        write( 102 , * ) '!!        dof_id = 2 , dof_var_id = 1  means function u_y. '
!        write( 102 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1  means function P. '
!        write( 102 , * ) '!!'
!        write( 102 , * ) '!! For nodes at reservior truncation boundary: HW model of order( N , M ): '
!        write( 102 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2        means function phi_1. '
!        write( 102 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1+j      means function phi_j. '
!        write( 102 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1+N+1    means function psi. '
!        write( 102 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2+N+j    means function phi_{N+j}. '
!        write( 102 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2+N+M+1  means function P^i ( induced pressure). '
!        write( 102 , * ) '!!'
!        write( 102 , * ) '!! Row_id = index of row in global matrix. '
!        write( 102 , * ) '!! col_id = index of col in global matrix that possibly has non-zero value ( up to 1E-13 ).'
!        write( 102 , * ) '!!'
!        write( 102 , * ) '!! ************************************************************************** !!'
!        write( 102 , * )
!        write( 102 , * )
!        write( 102 , * )
     
!        do  i = 1 , size( dictionary ) -1 
!            non_zero_index = 0
!            num_var = 0
!            do  k  = dictionary( i ) + 1 , dictionary( i + 1 )
!                if( abs( stiffness( k ) ) > epsi ) then
!                    num_var = num_var + 1
!                    non_zero_index(num_var) = k 
!                end if
!            end do
         
!            call  find_variable_identifier( i , node_id , dof_id , nodal_var_id )
!            write( 102 , * ) 
!            write( 102 , * ) ' ***************************** '
!            write( 102 , * ) '  Row id ,    Node_id ,     dof_id ,   dof_var_id   ' 
!            write( 102 , * ) i , node_id , dof_id , nodal_var_id 
          
!            write( 102 , * ) 'Maximum number of non-zero columns in this row:  ' , num_var
         
!            if( num_var .eq. 0 ) cycle 
         
!            write( 102 , * ) 'List of columns with potentially non-zero values: '
!            write( 102 , * ) '        col_id  ,  node_id  ,  dof_id , dof_var_id  ,  stiffness(row_id,col_id)  '
         
!            do k = 1 , num_var
!                j = dependency( non_zero_index(k) )
!                call  find_variable_identifier( j , node_id , dof_id , nodal_var_id )
!                write( 102 , * ) j , node_id , dof_id , nodal_var_id , '    '  , stiffness( non_zero_index(k) )
!            end do

 
     
!        end do
     
!        close( 102 )
     
!        open ( unit   = 103  , file    = "./output/damping_matrix.txt"  , action  = 'write'   ,  &
!               access = 'stream' , iostat  = state             , form    = 'formatted' )
     
!        if ( state /= 0) stop 'Error: cannot write output file in analysis module. '
     
!        write( 103 , * ) '!! **********************  Explanation of the file:  *********************** !!'
!        write( 103 , * ) '!!'
!        write( 103 , * ) '!! To every node (given by node_id) some dofs ( given by dof_id ) is assigned. '
!        write( 103 , * ) '!!        dof_id = 1 : displacement in X direction ( u_x )'
!        write( 103 , * ) '!!        dof_id = 2 : displacement in Y direction ( u_y )'
!        write( 103 , * ) '!!        dof_id = 3 : in (3D) displacement in Z direction ( u_z )'
!        write( 103 , * ) '!!                   : in 2D problem pressure ( P )  '
!        write( 103 , * ) '!!        dof_id = 4 : in 3D problem is pressure (P) '
!        write( 103 , * ) '!!'
!        write( 103 , * ) '!! Every dof is desribed by a set of variables ( given by dof_var_id ) '
!        write( 103 , * ) '!! Values of dof_var_id: '
!        write( 103 , * ) '!!        dof_id = 1 , dof_var_id = 1  means function u_x. '
!        write( 103 , * ) '!!        dof_id = 2 , dof_var_id = 1  means function u_y. '
!        write( 103 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1  means function P. '
!        write( 103 , * ) '!!'
!        write( 103 , * ) '!! For nodes at reservior truncation boundary: HW model of order( N , M ): '
!        write( 103 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2        means function phi_1. '
!        write( 103 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1+j      means function phi_j. '
!        write( 103 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 1+N+1    means function psi. '
!        write( 103 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2+N+j    means function phi_{N+j}. '
!        write( 103 , * ) '!!        dof_id = 3 (in 2D) , dof_var_id = 2+N+M+1  means function P^i ( induced pressure). '
!        write( 103 , * ) '!!'
!        write( 103 , * ) '!! Row_id = index of row in global matrix. '
!        write( 103 , * ) '!! col_id = index of col in global matrix that possibly has non-zero value ( up to 1E-13 ).'
!        write( 103 , * ) '!!'
!        write( 103 , * ) '!! ************************************************************************** !!'
!        write( 103 , * )
!        write( 103 , * )
!        write( 103 , * )
     
!        do  i = 1 , size( dictionary ) -1 
!            non_zero_index = 0
!            num_var = 0
!            do  k  = dictionary( i ) + 1 , dictionary( i + 1 )
!                if( abs( damping( k ) ) > epsi ) then
!                    num_var = num_var + 1
!                    non_zero_index(num_var) = k 
!                end if
!            end do
         
!            call  find_variable_identifier( i , node_id , dof_id , nodal_var_id )
!            write( 103 , * ) 
!            write( 103 , * ) ' ***************************** '
!            write( 103 , * ) '  Row id ,    Node_id ,     dof_id ,   dof_var_id   ' 
!            write( 103 , * ) i , node_id , dof_id , nodal_var_id 
          
!            write( 103 , * ) 'Maximum number of non-zero columns in this row:  ' , num_var
         
!            if( num_var .eq. 0 ) cycle 
         
!            write( 103 , * ) 'List of columns with potentially non-zero values: '
!            write( 103 , * ) '        col_id  ,  node_id  ,  dof_id , dof_var_id  ,  damping(row_id,col_id) '
!            if( num_var .eq. 0 ) cycle
         
!            do k = 1 , num_var
!                j = dependency( non_zero_index(k) )
!                call find_variable_identifier( j , node_id , dof_id , nodal_var_id )
!                write( 103 , * ) j , node_id , dof_id , nodal_var_id , '    ' , damping( non_zero_index(k) )
!            end do
!        end do
     
!        close( 103 )
     
!        deallocate( non_zero_index , stat = state )
!        if ( state /= 0) stop 'failed in deallocating non-zero_index vector in analysis module. '

    end subroutine store_in_file
 
   

!!================================================================
!!
!! For given variable find the node and dof to which the variable is assigned.
!!
!!
!!======================================================================
   
    subroutine find_variable_identifier( var_id ,  node_id , dof_id , nodal_var_id )
        use mod_geometry , only : degree_of_freedom , num_dim , num_degree_of_freedom
        use mod_physics  , only : the_model_description
      
        implicit none
        integer, intent( in ) :: var_id
        integer, intent( inout ) ::  node_id , dof_id ,  nodal_var_id
        integer :: i , d , k

      
      
        do  i = 1, size( degree_of_freedom  , 1 )
            if( degree_of_freedom(i , 1 ) .eq. 0 ) cycle        !! node has no dof
          
            do  d = 1 , num_degree_of_freedom
                do k = degree_of_freedom( i , 2*d ) , degree_of_freedom( i , 2*d + 1 )
                    if( var_id .eq. k ) then
                        node_id   = i 
                        dof_id = d 
                        nodal_var_id = k - degree_of_freedom( i , 2*d ) + 1
                        return                                             
                    end if 
                end do
            end do 
        end do ! over i
   
        stop 'failed in finding node and dof for given variable in analysis module'
      
     
      
   
    end subroutine find_variable_identifier
   
!!=================================================================
!!
!!   clean up analysis
!!
!!================================================================= 

 
    subroutine  clean_up_analysis()
     
        use mod_geometry , only : clean_up_geometry
        use mod_physics  , only : clean_up_physics
        use mod_fem      , only : destruct_fem_model
        use mod_solver   , only : clean_up_solver
        implicit none
        integer :: state = 0 ;
        
       
        call  clean_up_geometry()
        call  clean_up_physics() 
        call  destruct_fem_model()                                      !! free memory
        call  clean_up_solver() ;
   
        ! if analysis is in time domain then clear force data
        if( allocated( reproted_node_index_timeD ) ) then
            deallocate( reproted_node_index_timeD, stat = state ) ;
            if( state /= 0 ) stop 'Failed in deallocating memory in clean_up_analysis in analysis module' ;
        end if
        
        
        if( allocated( reported_time ) ) then
            deallocate( reported_time , stat = state ) ;
            if( state /= 0 ) stop 'Failed in deallocating memory in clean_up_analysis in analysis module' ;
        end if
        
        
        
    end subroutine clean_up_analysis
  
  


   


    
 
  
 end module mod_analysis
  
