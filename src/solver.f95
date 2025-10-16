!------------------------------------------------------------------------
!   Created by: Nili Abtahi
!
!   Laboratory for Computational Sensing and Robotics, John Hopkins University
!
!   Contact: Nili Abtahi (nabtahi1@jhu.edu)
!
!----------------------------------------------------------------------!



 module   mod_solver
      
     use  mod_utils , only : r_kind , qp_kind
     implicit none
     
     
            
!!====================================================================!!
!!  variables for complex linear system
      
    complex( kind = r_kind ) , dimension(:) , allocatable , public  ::  cplx_matrix   ! matrix A in equation A * x = b
    complex( kind = r_kind ) , dimension(:) , allocatable , public  ::  cplx_response ! vector b in equation A * x = b
    complex( kind = r_kind ) , dimension(:) , allocatable , public  ::  cplx_solution ! vector x in equation A * x = b
    
    
!!====================================================================!!
!!  variables for real linear system
    
    real( kind = r_kind ) , dimension(:) , allocatable , public  ::  matrix   ! matrix A in equation A * x = b
    real( kind = r_kind ) , dimension(:) , allocatable , public  ::  response ! vector b in equation A * x = b
    real( kind = r_kind ) , dimension(:) , allocatable , public  ::  solution ! vector x in equation A * x = b
     
     
    !! used for direct substitution method
    
    integer, dimension(:) , allocatable , private :: solver_dictionary
    integer, dimension(:) , allocatable , private :: index_set 
    integer, dimension(:) , allocatable , private :: solver_dependency       !! main dependency matrix 
    integer, dimension(:) , allocatable , private :: auxiliary_dependency    !! auxiliary dependency matrix
    integer, dimension(:) , allocatable , private :: permutation
    integer, private :: num_row_in_main_matrix
   

    !!================================================================!!
    !! complex solver
    
    complex( kind = r_kind ) , dimension(:) , allocatable , private   :: cplx_upper_matrix       !! main upper matrix 
    complex( kind = r_kind ) , dimension(:) , allocatable , private   :: cplx_auxiliary_matrix   !! auxiliary upper matrix
    complex( kind = r_kind ) , dimension(:) , allocatable , private   :: cplx_temp_row
    complex( kind = r_kind ) , dimension(:) , allocatable , private   :: cplx_transformed_response
   
    !!==================================================================!!
    !! real solver
    
    real( kind = r_kind ) , dimension(:) , allocatable , private   ::   upper_matrix       !! main upper matrix 
    real( kind = r_kind ) , dimension(:) , allocatable , private   ::   auxiliary_matrix   !! auxiliary upper matrix
    real( kind = r_kind ) , dimension(:) , allocatable , private   ::   temp_row
    real( kind = r_kind ) , dimension(:) , allocatable , private   ::   transformed_response
    real( kind = r_kind ) , allocatable  , dimension(:), public    ::   u_vec     !! displacement vector
    real( kind = r_kind ) , allocatable  , dimension(:), public    ::   udot_vec  !! velocity vector
    real( kind = r_kind ) , allocatable  , dimension(:), public    ::   uddot_vec !! acceleration vector
    real( kind = r_kind ) , allocatable  , dimension(:), public    ::   temp_real_vec !! a temp vector
    real( kind = r_kind ) , allocatable  , dimension(:), public    ::   static_solution  !! solution at static step
    
    
    real( kind = r_kind ) , parameter :: effective_zero   = 1.0E-13  ;
    
     
!!====================================================================!!
!! subroutines and functions

   !! public functions
   public :: initialize_complex_solver
   public :: clean_up_solver
   public :: solve_complex_system
   
   public :: initialize_real_solver
   public :: solve_real_system
   
   

   private :: add_new_row_to_cplx_matrix 
   private :: add_real_row_to_matrix
   private :: get_location
   
  contains 
  
    
 



!!======================================================================
!!
!!   solve_real_system by direct substitution method
!!   
!!
!!====================================================================== 
    
    
     
subroutine  solve_real_system(istep)
    use  mod_physics  , only : dictionary, dependency
    use  mod_geometry , only : degree_of_freedom 
    implicit none  
    integer, intent(in) :: istep
    integer :: i , j , k, row_id
    real( kind = r_kind ) ::  max_error , tmp , max_val 
    logical :: is_printed = .false.
    
    !if( istep .eq. 8 ) is_printed = .true. 
    
    if( is_printed ) then
        open ( unit   = 100  , file    = "./output/debug.txt"  , action  = 'write'   ,  &
            access = 'stream' , iostat  = i             , form    = 'formatted' )
    end if
    
    !! 0. normalize system
    do  i = 1 , size( dictionary ) -1
        tmp = 0.0_r_kind  ;
        do  j = dictionary(i)+1 , dictionary(i+1)
            tmp = tmp + matrix(j) ** 2  ;
        end do
        
        tmp =  sqrt( tmp ) ;
        
        if( tmp < effective_zero ) then
            write(*,*) 'The matrix is singular as norm of row id = ' , i , ' is = ' , tmp  ;
            stop 
        end if
        
        if( is_printed ) then
            write(100 , *) i , tmp
        end if
        
        do  j = dictionary(i)+1 , dictionary(i+1)
            matrix(j) = matrix(j) / tmp  ;
        end do
        response(i) = response(i) / tmp ;
    end do
    
    
    
    !! 1. set initial values
    !! Row permutation vectors are initially identity matrix, store index of non-zero column at each row
    do i = 1 , size( permutation ) 
       permutation(i) = i ;
    end do
    
    
    num_row_in_main_matrix  = 0 ;                                        !! Is the main matrix enough for storing transformed base matrix?
    solver_dependency(:)    = 0 ;                                        !! Initial value of dependency matrix is zero
    upper_matrix(:)         = 0.0_r_kind  ;                              !! This is also zero vector
    solver_dictionary(:)    = 0 ;                                        !! The same for solver dictionary
    transformed_response(:) = 0.0_r_kind  ;    
   

    !! 2. trasnform matrix to an upper triangular one
    !! 2.1. find the k-th row that has max x_i value:
    
    max_val = 0.0_r_kind ;
    row_id = 0 ;
    do  i = 1 , size( dictionary ) - 1                                  !! Loop over all rows
        if( dependency( dictionary(i)+1 ) > 1 ) cycle                   !! 2.1.1. If x1 is not a member of this row just ignore the row.    
                
        if( abs(matrix( dictionary(i)+1)) > max_val )  then             !! 2.1.2. check if the matrix A(i,1) is > max_val
            max_val = abs(matrix( dictionary(i)+1)) ;
            row_id  = i ;
            if( max_val > 0.10_r_kind ) exit ;
        end if 
    end do
    
    
    if( row_id .eq. 0 ) stop 'Real Solver: Singular matrix first column has no non-zero value.' ;
    
    permutation(row_id) = 1 ;                                       !! 2.1.4. permute i-th row with first one
    permutation(1) = row_id ;                                       !! 2.1.5. permute first row with i-th one
    
    !! 2.1.8. First row may have been replaced by permutation(1)-th row, so the true first row is j = permutation(1).
    !!        Put this row in the first row of upper triangular matrix.
    
   
    solver_dictionary(2) = dictionary( permutation(1) + 1) -dictionary( permutation(1) ) ;  !! note that solver_dictionary(1) = 0
    solver_dependency(1:solver_dictionary(2))   = dependency(   dictionary( permutation(1) ) + 1  :  &
                                                            &   dictionary( permutation(1) + 1 )) ;
                                                            
    upper_matrix(1:solver_dictionary(2)) = matrix(  dictionary( permutation(1) ) + 1  :  &
                                                 &  dictionary( permutation(1) + 1 )) ;
    
    transformed_response(1) = response( permutation(1) ) ;      !! put response vector to transformed response without permutation
    
    num_row_in_main_matrix  = 1 ;                               !! Is the main matrix enough for storing transformed base matrix? 


    if( is_printed ) then
        write(100 , *) '============================ row =' , 1  
        write(100 , *) 1 , permutation(1) , solver_dictionary(1:2)
        do  i = solver_dictionary(1) + 1 , solver_dictionary(2) 
            write(100 , *) solver_dependency(i) , upper_matrix(i)
        end do  
    end if 
        
    !! 2.2. add remaining rows to the upper matrix
    
    temp_row(:)   = 0.0_r_kind  ;
    index_set(:)  = 0 ;
    
    
    do  row_id = 2 , size( dictionary )-1           !! 2.2.1. Loop over all remaining rows: from row= 2 to row = n
        call   add_real_row_to_matrix( row_id ,is_printed ) ;   !! 2.2.2. transform and add rows to the upper matrix 
    end do
    
    
!    write(*,*) 'num_row_in_main_matrix = ' , num_row_in_main_matrix , size(permutation)
    
!    write(*,*) 'size(dep) ' , size( dependency ) 
!    write(*,*) 'size solv dep = ' , size( solver_dependency )
!    if( allocated( auxiliary_dependency ) ) then
!        write(*,*) 'size aux dep = ' , size( auxiliary_dependency )
!    else
!        write(*,*) ' auxilary arrays are not allocated' 
!    end if
!    write(*,*) 'max vals in dict = ' , dictionary( size(dictionary) )
!    write(*,*) 'max vals in solver dic = ' , solver_dictionary( size(solver_dictionary) )
    

    
!    open ( unit   = 103  , file    = "./output/upperMatrix.txt"  , action  = 'write'   ,  &
!            access = 'stream' , iostat  = i             , form    = 'formatted' )
    
    
!    j = size( solver_dictionary ) -1
!    do i = 1 , j
!        temp_row(1:j) = 0.0_r_kind ;
!        do  k = solver_dictionary(i) + 1 , solver_dictionary(i+1)
!            temp_row(solver_dependency(k)) = upper_matrix(k) ;
!        end do
!        write(103,*) temp_row(1:j)
!    end do
    
!    close(103) ;
    
!    open ( unit   = 104  , file    = "./output/upperResponse.txt"  , action  = 'write'   ,  &
!            access = 'stream' , iostat  = i             , form    = 'formatted' )
    
!    write(104,*) transformed_response
!    close(104) ;
    
     
    
    solution(:) =  0.0_r_kind  ;
    solution = transformed_response ;
    
    do  i = size( permutation ) , num_row_in_main_matrix  + 2 , -1
        solution(i) = transformed_response(i) ;
        do  j = solver_dictionary(i)+ 2  , solver_dictionary(i+1)
            solution(i) = solution(i) - auxiliary_matrix(j) * solution(auxiliary_dependency(j)) ;
        end do 
        solution(i) = solution(i)/ auxiliary_matrix( solver_dictionary(i) + 1 ) ;
    end do 
    
    !! 3.2. the first item in new database begins from id = 1
    if( size( permutation ) > num_row_in_main_matrix ) then
        i = num_row_in_main_matrix  + 1 ;
        solution(i) = transformed_response(i) ;
        do  j = 2  , solver_dictionary(i+1)
            solution(i) = solution(i) - auxiliary_matrix(j) * solution(auxiliary_dependency(j)) ;
        end do 
        solution(i) = solution(i)/ auxiliary_matrix(1 ) ;
    end if 
    
    
    ! 3.2. second for remaining variables
    do  i = num_row_in_main_matrix , 1 , -1
        
        do  j = solver_dictionary(i)+ 2  , solver_dictionary(i+1)
            solution(i) = solution(i) - upper_matrix(j) * solution( solver_dependency(j) ) ;
        end do 
        solution(i) = solution(i) / upper_matrix( solver_dictionary(i) + 1 ) ;
    end do
    
    
    if( is_printed ) then
        write(100 , *) '============================  solution '  
    end if
        
    ! 3.3 Testing solution
    temp_row(:)   = 0.0_r_kind  ;
    index_set(:)  = 0 ;
    max_error     = 0.0_r_kind ;
    
    do i = 1 , size( permutation )
        do j = dictionary(i) + 1 , dictionary(i+1)
           temp_row(i) = temp_row(i) + matrix(j) * solution( dependency(j) ) ;
        end do
        temp_row(i) = temp_row(i) - response(i) ;
        max_error  = max( max_error  , abs( temp_row(i) ) )
        
        if( is_printed ) then
            write(100 , *) i , temp_row(i) 
        end if
    end do
    
    if( is_printed ) then
        close(100)
    end if
    
    write(*,*)  istep  , max_error 

    
end subroutine  solve_real_system
  

!!=================================================================
!!
!!   for variable x_k (k = row_id), choose j-th rows of matrix A, such that A(k,k) /= 0 after substitution variables x_i, i = 1 , ..., k.
!!
!!================================================================= 


    subroutine  add_real_row_to_matrix( row_id , is_printed )  
    
        use  mod_physics , only : dictionary , dependency
        implicit none
        integer , intent( in  ) :: row_id 
        logical , intent( in ) :: is_printed
        integer :: i, j, k, k_beg , k_end , num_cols, beg_id , end_id , id_min , id_max, first_var
        logical :: found_any  ;
        real( kind = r_kind ) :: response_value, tmp 
      
        
        !! 1. loop over all unused rows that contain x_k, the unused rows are permutation( row_id:end)
        
        found_any = .false. ;
        
        do  i = row_id ,  size( dictionary ) - 1 
        
            !! 2. put i-th row in the index set and matrix row
            
            beg_id = dictionary( permutation( i ) ) + 1 ;
            end_id = dictionary( permutation( i ) + 1 ) ;
        
            !! 3. if the row does not contain x_k variable, just ignore it
            
            if( ( dependency(beg_id) > row_id ) .or. ( dependency(end_id) < row_id ) ) cycle ;  !! simple necessary condition
            

            if( get_location( permutation( i ) , row_id ) .eq. 0 ) cycle ;                          !! sufficient condition for checking if the variable is in the row
            
            id_min = dependency(beg_id) ;
            id_max = dependency(end_id) ;
            
            
            do j = beg_id , end_id
               index_set( dependency(j) ) = 1 ;                         !! these variables are present in this row
               temp_row( dependency(j) )  = matrix(j) ;                 !! put matrix row in the row-vector
            end do 
            
            first_var = id_min ;
            
            response_value = response( permutation( i ) ) ;
            
            
            !! 3. substitute all previously singled out variables, i.e. x_1 from first row, x_2 from second row and so on.
            !! Some variables are in the base matrix and others are stored in the auxiliary matrix, 
            !! So the substitution is separated to two parts, i.e. variables storing in base and auxiliary matrices.
            !! The whole variables in this row are dependency(beg_id:end_id), all variables before diagonal one, i.e. x_k, k= row_id
            !! should be substituted. Therefore, the following loop begins from dependency(beg_id) and runs to dependency( varlocRowId - 1 )
            !! Some are in the base matrix the others are in the auxiliary matrix.
            
            !! substituttion loop 
            !! There are three possibilitieis
            !!      1. if  num_row_in_main_matrix = row_id -1, the base matrix is under considration
            !!      2. if  num_row_in_main_matrix = row_id -2, the auxiliary matrix is chosen for the first time
            !!      3. otherwise the auxiliary matrix is under consideration
            
            !! the substitution loop should be over index set itself as substituting x_1 may introduce x_3 that has not been in the
            !! main index set. The beginning value in the index set is id_min and it runs to row_id -1 in the index set.
            !! if some index was zero then it is not calculated.
            
            
            do  j = first_var , row_id - 1
               
                if( index_set(j) /= 1 ) cycle ;
                 
                if( j .le. num_row_in_main_matrix ) then                !! 3.1. the j-th row is in base matrix
                
                    k_beg = solver_dictionary( j ) + 1 ;
                    k_end = solver_dictionary( j + 1 ) ;
                    
                    id_min = min( id_min , solver_dependency(k_beg) ) ;
                    id_max = max( id_max , solver_dependency(k_end) ) ;
                    
                    
                    tmp =   temp_row( j )  /  upper_matrix( k_beg )  ;
                    
                    
                    do k = k_beg , k_end
                        index_set( solver_dependency(k) ) = 1 ;
                        temp_row(solver_dependency(k)) = temp_row(solver_dependency(k))- tmp * upper_matrix(k)
                    end do  
                    
               
                else    !! 3.2. The j-th row is stored in auxiliary matrix
                    
                    !! 3.2.1. In the second case mentioned above, k_beg =1 , i.e. the first item in auxiliary matrix
                    if( j .eq. num_row_in_main_matrix + 1 ) then
                        k_beg = 1 ;
                    else
                        k_beg = solver_dictionary( j ) + 1 ;
                    end if 
                    
                    k_end = solver_dictionary( j + 1 )  ;
                    
                    id_min = min( id_min , auxiliary_dependency(k_beg) ) ;
                    id_max = max( id_max , auxiliary_dependency(k_end) ) ;
                    
                    tmp =  temp_row( j ) / auxiliary_matrix( k_beg );
                             
                    
                    do k = k_beg , k_end
                        index_set( auxiliary_dependency(k) ) = 1 ;
                        temp_row(auxiliary_dependency(k)) =   temp_row(auxiliary_dependency(k)) &
                                                          & - tmp * auxiliary_matrix( k ) ;
                    end do   
                end if
            
                !! 3.3. change response value 
                response_value = response_value - tmp * transformed_response(j) ;
                
            end do !! over j
            

            
            
            !! 4. check if A(k,k) /= 0, otherwise set vectors to default values
            if( abs( temp_row( row_id ) ) /= 0.0_r_kind ) then ! effective_zero ) then
            
                found_any = .true.
                permutation( row_id ) = i ;                                  !! 4.1. do suitable permutaion of rows
                permutation( i ) = row_id ;
                transformed_response(row_id) = response_value ;         !! 4.2. add transformed response to the vector
                exit ;                                                       !! 4.3. exit the i-loop
                
            else
            
                index_set( id_min : id_max ) = 0 ;
                temp_row( id_min:id_max ) = 0.0_r_kind  ;
            
            end if 
        
        end do ! over i
        
         
        !! 5. check if the variable x_k was singled out
        if( .not. found_any ) then
            write(*,*) 'Error: solver: the variable x_k for k = ' , row_id , ' cannot be singled out, singular matrix' ;
            stop 'Error: solver: Matrix is singular' 
        end if 
       
        
        
        !! 6. add new row to upper triangular matrix if there is enough memory
        !!  if num_row_in_main_matrix = RowId - 1 then it is checked whether there is room to add new elements.
        !!  If so, num_row is incremnted and the difference remains the same for next iteration.
        !!  Otherwise, the new auxiliary matrices are allocated and num_row is not changed. 
        !!  For other iterations, RowId - numRow > 1 and the new matrices are used for storing data
        
        if( row_id .eq. num_row_in_main_matrix + 1 ) then
            
            !! 6.1. calculate required space
            num_cols = 0 ;
            do  i = row_id , id_max
                if( index_set(i) .eq. 1 ) then
                    num_cols = num_cols + 1 ;
                end if
            end do 
            
            
            !! 6.2. calculate available space in main matrix, put it to beg_id
            beg_id = size( solver_dependency ) - solver_dictionary( row_id ) ;
            
            !! 6.3. if there is enough room, just add this row to main matrix
            if( beg_id .ge. num_cols ) then
                
                k = solver_dictionary( row_id ) ;                       !! 6.3.1. get location of new variables in solver dictionary
                
                
                do  i = row_id , id_max
                    if( index_set(i) .eq. 1 ) then
                        k = k + 1 ;
                        solver_dependency( k ) = i ;
                        upper_matrix( k ) = temp_row(i) ;
                    end if
                end do 
                
                
                num_row_in_main_matrix = row_id ;                       !! 6.3.2. increment num rows in main matrix
                solver_dictionary( row_id + 1 )  = k ;                  !! 6.3.3. set next value of dictionary
                
                if( is_printed ) then
                    write(100 , *) '============================ row =' ,row_id  
                    write(100 , *) row_id , permutation(row_id) , solver_dictionary(row_id:row_id+1)
                    do  i = solver_dictionary(row_id) + 1 , solver_dictionary(row_id+1) 
                        write(100 , *) solver_dependency(i) , upper_matrix(i)
                    end do  
                end if
                 
            else  
                !! 6.4. otherwise, allocated memory for new matrices
                !! 6.4.1. Calculate required memory for next matrix to be fully filled by other rows
                !!        It is (N - m )*( N -m + 1 ) / 2 where m = num_row_in_main_matrix
                
                k=  size( permutation ) - num_row_in_main_matrix  ;
                
                num_cols =  ( k * ( k + 1 ) ) / 2 ;                         !! find required memory
                
                !! 6.4.2. If auxilary ariables are not allocated just allocate them.
                
                if( .not. allocated( auxiliary_matrix ) ) then
                
                    allocate( auxiliary_matrix( int( 1.2 * num_cols ) ) , stat = k )   !! rescale memory by an offset
                    if( k /= 0 ) stop 'Error: failed in allocating cplx_auxiliary_matrix in solver module.' 
                
                    allocate( auxiliary_dependency( int( 1.2 * num_cols ) ) , stat = k )
                    if( k /= 0 ) stop 'Error: failed in allocating auxiliary_dependency in solver module.' 
                    
                else  !! 6.4.3. If auxiliary matrices are allocated just check if their size is enough large, as required.
                   
                    if( size( auxiliary_matrix ) < num_cols ) then
                    
                        deallocate( auxiliary_matrix , stat = k )
                        deallocate( auxiliary_dependency , stat = k )
                        if( k /= 0 ) stop ' Error in deallocating auxiliary arrays in solver module'
                        
                        allocate( auxiliary_matrix( int( 1.2 * num_cols ) ) , stat = k )
                        if( k /= 0 ) stop 'Error: failed in allocating cplx_auxiliary_matrix in solver module.' 
                
                        allocate( auxiliary_dependency( int( 1.2 * num_cols ) ) , stat = k )
                        if( k /= 0 ) stop 'Error: failed in allocating auxiliary_dependency in solver module.' 
                    end if
                    
                end if 
                
                auxiliary_dependency(:) = 0 
                auxiliary_matrix(:) =  0.0_r_kind  ;
                
                
                !! 6.4.4. Store data in auxiliary matrices from their first elements
                k = 0 ;
             
                do  i = row_id , id_max
                    if( index_set(i) .eq. 1 ) then
                        k = k + 1 ;
                        auxiliary_dependency( k ) = i ;
                        auxiliary_matrix( k ) = temp_row(i) ;
                    end if
                end do 
                
                solver_dictionary( row_id + 1 )  = k ;                  !! 8.3.3. set next value of dictionary
                
                if( is_printed ) then
                    write(100 , *) '============================ row =' ,row_id  
                    write(100 , *) row_id , permutation(row_id) , solver_dictionary(row_id:row_id+1)
                    do  i = solver_dictionary(row_id) + 1 , solver_dictionary(row_id+1) 
                        write(100 , *) auxiliary_dependency(i) , auxiliary_matrix(i)
                    end do  
                end if
                
            end if 
            
        else
        
            !! 6.5. add new data to auxiliary matrices
            k = solver_dictionary( row_id )  ;
             
            do  i = row_id , id_max
                if( index_set(i) .eq. 1 ) then
                    k = k + 1 ;
                    auxiliary_dependency( k ) = i ;
                    auxiliary_matrix( k ) = temp_row(i) ;
                end if
            end do 
            
            solver_dictionary( row_id + 1 )  = k ;                          !! 8.3.3. set next value of dictionary
            
            if( is_printed ) then
                write(100 , *) '============================ row =' ,row_id  
                write(100 , *) row_id , permutation(row_id) , solver_dictionary(row_id:row_id+1)
                do  i = solver_dictionary(row_id) + 1 , solver_dictionary(row_id+1) 
                    write(100 , *) auxiliary_dependency(i) , auxiliary_matrix(i)
                end do  
            end if
        end if
                
        
        index_set( id_min : id_max ) = 0 ;
        temp_row( id_min:id_max ) =  0.0_r_kind  ;        
    
    end subroutine  add_real_row_to_matrix
    
    
!    subroutine  add_real_row_to_matrix( row_id , istep )  
    
!        use  mod_physics , only : dictionary , dependency
!        implicit none
!        integer , intent( in  ) :: row_id , istep
!        integer :: i, j, k, k_beg , k_end , numvars , first_var, max_id, current_id, beg_id , end_id, id_max , id_min
!        real( kind = r_kind ) :: response_value, tmp , tmp2 , max_val , max_response ;
!        logical :: is_test = .false.
        
!        max_id = 0 ;
!        max_val = 0.0_r_kind ;
        
!        !if( (istep .eq. 3) .and. (row_id .eq. 4743) ) is_test = .true.
        
!        if( is_test) then
!            !! print data in new index set
!            open ( unit   = 88  , file    = "./output/FEMData1.txt"  , action  = 'write'   ,  &
!                   access = 'stream' , iostat  = i             , form    = 'formatted' )

!            if ( i /= 0) stop 'Error: cannot write output file in analysis module. '
            
!            write(88,*) ' row_id = ' , row_id
!        end if 
        
!        id_min = 1 ;
!        id_max = size( dictionary ) -1 ;
        
!        !!  1. find suitable row or column
!        do  i = row_id , size( dictionary ) -1
            
!            index_set( id_min : id_max ) =  0 ;
!            temp_row(  id_min:id_max   ) =  0.0_r_kind  ;
            
!            current_id = i ;
!            !! 2. put i-th row in the index set and matrix row
!            beg_id = dictionary( permutation( i ) ) + 1 ;
!            end_id = dictionary( permutation( i ) + 1 ) ;
        
            
!            if( is_test) then
!               write(88,*) 'i = ' , i , 'beg_id =' , beg_id , 'end_id =' , end_id , &
!                         & 'dep(beg_id) =' , dependency(beg_id) , ' P(i) =' , permutation(i) 
!            end if
            
!            !! 3. if the row cannot possibly generate x_k variable, just ignore it:
!            !!    If all of its variables are > row_id then after substitution the coeff is zero
            
!            if( dependency(beg_id) > permutation(row_id)  ) cycle ; 
            
!            !! The varialble row_id exists in row id = permutation( i )
!            id_min = dependency(beg_id) ;
!            id_max = dependency(end_id) ;
            
            
            
!            do j = beg_id , end_id
!               index_set( dependency(j) ) = 1 ;                         !! these variables are present in this row
!               temp_row(  dependency(j) ) = matrix(j) ;                 !! put matrix row in the row-vector
!            end do 
            
!            first_var = id_min ;
            
!            response_value = response( permutation( i ) ) ;
            
            
!            if( is_test) then
!                k = 0 ;
!                do j = 1 , size(temp_row)
!                    if( index_set(j) /= 0) then
!                        write(88,*) j, temp_row(j)
!                        k = k + 1 ;
!                    end if
!                end do
!                write(88,*) 'nnz = ' , k , 'end-beg = ' , end_id - beg_id + 1 , 'Resp =' , response_value
!            end if
            
!!            if( is_test ) then
!!                write(*,*) 'i = ' , i , 'P(i) = ' , permutation(i) , 'tmp_row_id=' , temp_row(row_id)
!!            end if
            
            
                                   
!            !! 3. substitute all previously singled out variables, i.e. x_1 from first row, x_2 from second row and so on.
!            !! Some variables are in the base matrix and others are stored in the auxiliary matrix, 
!            !! So the substitution is separated to two parts, i.e. variables storing in base and auxiliary matrices.
!            !! The whole variables in this row are dependency(beg_id:end_id), all variables before diagonal one, i.e. x_k, k= row_id
!            !! should be substituted. Therefore, the following loop begins from dependency(beg_id) and runs to dependency( varlocRowId - 1 )
!            !! Some are in the base matrix the others are in the auxiliary matrix.
            
!            !! substituttion loop 
!            !! There are three possibilitieis
!            !!      1. if  num_row_in_main_matrix = row_id -1, the base matrix is under considration
!            !!      2. if  num_row_in_main_matrix = row_id -2, the auxiliary matrix is chosen for the first time
!            !!      3. otherwise the auxiliary matrix is under consideration
            
!            !! the substitution loop should be over index set itself as substituting x_1 may introduce x_3 that has not been in the
!            !! main index set. The beginning value in the index set is id_min and it runs to row_id -1 in the index set.
!            !! if some index was zero then it is not calculated.
            
!            if( is_test) then
!                write(88,*) '----------------- substitution step -------------'
!            end if
            
            
!            do  j = first_var , row_id - 1
               
                
            
!                if( index_set(j) /= 1 ) cycle ;
!                if( temp_row(j) .eq. 0.0_r_kind ) cycle ;
                 
!                if( j .le. num_row_in_main_matrix ) then                !! 3.1. the j-th row is in base matrix
                    
!                    k_beg = solver_dictionary( j ) + 1 ;
!                    k_end = solver_dictionary( j + 1 ) ;
                    
!                    id_min = min( id_min , solver_dependency(k_beg) ) ;
!                    id_max = max( id_max , solver_dependency(k_end) ) ;
                
                             
!                    tmp = temp_row( j ) / upper_matrix( k_beg ) ;          
                    
                    
            
!                    do k = k_beg , k_end
!                        index_set(solver_dependency(k)) = 1 ;
!                        temp_row( solver_dependency(k)) = temp_row(solver_dependency(k)) - tmp * upper_matrix(k)
!                    end do  
                    
!                    if( is_test) then
!                        write(88,*) 'subs j = ' , j , 'upper_matrix( k_beg ) =' , upper_matrix( k_beg ) , &
!                                  &  'tmp = ' , tmp , 'temp_row_mid = ' , temp_row(row_id)
!                    end if
                    
               
!                else    !! 3.2. The j-th row is stored in auxiliary matrix
                    
!                    !! 3.2.1. In the second case mentioned above, k_beg =1 , i.e. the first item in auxiliary matrix
!                    if( j .eq. num_row_in_main_matrix + 1 ) then
!                        k_beg = 1 ;
!                    else
!                        k_beg = solver_dictionary( j ) + 1 ;
!                    end if 
                    
!                    k_end = solver_dictionary( j + 1 )  ;
                    
!                    id_min = min( id_min , auxiliary_dependency(k_beg) ) ;
!                    id_max = max( id_max , auxiliary_dependency(k_end) ) ;
                    
!                    tmp = temp_row(j) / auxiliary_matrix( k_beg ) ;
                    
!                    do k = k_beg , k_end
!                        index_set( auxiliary_dependency(k)) = 1 ;
!                        temp_row(  auxiliary_dependency(k)) = temp_row(auxiliary_dependency(k)) &
!                                                             & - tmp * auxiliary_matrix( k ) ;
!                    end do   
!                end if
                
        
!                !! 3.3. change response value 
!                response_value = response_value - tmp * transformed_response(j) ;
                
!            end do !! over j
            
!!            if( is_test ) then
!!                write(*,*) 'i= ' , i , 'P(i)=', permutation(i), 'tmp_row_id=' , temp_row(row_id) , max_val , max_id
!!            end if
            
!            if( is_test) then
!                write(88,*) 'max_val = ' , max_val , 'max_id =' , max_id , &
!                                  &  'temp_row(row_id) = ' , temp_row(row_id)
!            end if
            
!            !! choose max value
!            if( abs(temp_row(row_id)) > max_val ) then
!                max_id  = i  ;
!                max_val = abs(temp_row(row_id )) ;
!                max_response = response_value ;
!            end if 
            
!!            if( is_test ) then
!!                write(*,*) 'i= ' , i , 'max val=' , max_val , ' max_id = ' , max_id
!!            end if
                
!            if( max_val .ge. 0.10_r_kind ) exit ;
             
!        end do !! over i 
        
!        if( is_test) then
!            write(88,*) 'current_id = ' , current_id , 'max_id =' , max_id 
!        end if
        
!        !!------------------------------------------------------------!!
!        !! 4. if row of max value is different from max_id then use it
!        if( current_id /= max_id  ) then
            
!            index_set( id_min:id_max ) =  0 ;
!            temp_row(  id_min:id_max ) =  0.0_r_kind  ;
            
!            beg_id = 0 ;
!            end_id = 0 ;
            
!            do j = 1 , size( index_set )
!                if( index_set(j) /= 0 ) then
!                    beg_id = beg_id + 1 ;
!                end if 
!            end do
            
!            do j = 1 , size( temp_row )
!                if( temp_row(j) /= 0.0_r_kind ) then
!                    end_id = end_id + 1 ;
!                end if 
!            end do
            
!            beg_id = dictionary( permutation( max_id ) ) + 1 ;
!            end_id = dictionary( permutation( max_id ) + 1 ) ;
        
!            !! 4.1. The varialble row_id exists in row id = permutation( i )
!            id_min = dependency(beg_id) ;
!            id_max = dependency(end_id) ;
            
            
!            do j = beg_id , end_id
!               index_set( dependency(j) ) = 1 ;                         !! these variables are present in this row
!               temp_row(  dependency(j) ) = matrix(j) ;                 !! put matrix row in the row-vector
!            end do 
            
!            first_var = id_min ;
!            response_value = response( permutation( max_id ) ) ;
           
!            do  j = first_var , row_id - 1
               
!                if( index_set(j) /= 1 ) cycle ;
!                if( temp_row(j) .eq. 0.0_r_kind ) cycle ;
                 
!                if( j .le. num_row_in_main_matrix ) then                !! 3.1. the j-th row is in base matrix
                    
!                    k_beg = solver_dictionary( j ) + 1 ;
!                    k_end = solver_dictionary( j + 1 ) ;
                    
!                    id_min = min( id_min , solver_dependency(k_beg) ) ;
!                    id_max = max( id_max , solver_dependency(k_end) ) ;
                     
!                    tmp = temp_row( j ) / upper_matrix( k_beg ) ;          
                    
!                    do k = k_beg , k_end
!                        index_set(solver_dependency(k)) = 1 ;
!                        temp_row( solver_dependency(k)) = temp_row(solver_dependency(k)) - tmp * upper_matrix(k)
!                    end do  
                    
!                else    !! 4.2. The j-th row is stored in auxiliary matrix
                    
!                    !! 4.2.1. In the second case mentioned above, k_beg =1 , i.e. the first item in auxiliary matrix
!                    if( j .eq. num_row_in_main_matrix + 1 ) then
!                        k_beg = 1 ;
!                    else
!                        k_beg = solver_dictionary( j ) + 1 ;
!                    end if 
                    
!                    k_end = solver_dictionary( j + 1 )  ;
                    
!                    id_min = min( id_min , auxiliary_dependency(k_beg) ) ;
!                    id_max = max( id_max , auxiliary_dependency(k_end) ) ;
                    
!                    tmp = temp_row(j) / auxiliary_matrix( k_beg ) ;
                    
!                    do k = k_beg , k_end
!                        index_set( auxiliary_dependency(k)) = 1 ;
!                        temp_row(  auxiliary_dependency(k)) = temp_row(auxiliary_dependency(k)) &
!                                                             & - tmp * auxiliary_matrix( k ) ;
!                    end do   
!                end if
      
!                response_value = response_value - tmp * transformed_response(j) ;  
                
!            end do !! over j
        
!            !! 4.3. now find max column
!            max_response = response_value ;
!        end if 
        
!        if( is_test) then
!            write(88,*) 'max_val = ' , max_val , 'max_id =' , max_id , &
!                     &  'temp_row(row_id) = ' , temp_row(row_id)
            
!            close(88) ;
!        end if
        
!        if( (max_id .eq. 0) .or. ( abs(temp_row(row_id )) .eq. 0.0_r_kind) ) then
!            write(*,*) 'The variable x_{k} for k = ' , row_id , ' cannot be singled out value =', temp_row(row_id);
!            stop
!        end if
        
!        !! change permutation matrix
!        beg_id = permutation( row_id ) ;
!        permutation( row_id ) = permutation( max_id ) ;                 !! 4.1. do suitable permutaion of rows
!        permutation( max_id ) = beg_id ;
!        transformed_response(row_id) = max_response ;
        
!!        write(*,*) 'i = ' , row_id , 'P(i) = ' , permutation( row_id ) ,  &
!!                  & ' diag = ' , temp_row(row_id ) , 'max_response = ' , max_response 
        
!        !!=============================================================!!
!        !! 5. Add to upper triangular matrix
        
            
!        !! 6. add new row to upper triangular matrix if there is enough memory
!        !!  if num_row_in_main_matrix = RowId - 1 then it is checked whether there is room to add new elements.
!        !!  If so, num_row is incremnted and the difference remains the same for next iteration.
!        !!  Otherwise, the new auxiliary matrices are allocated and num_row is not changed. 
!        !!  For other iterations, RowId - numRow > 1 and the new matrices are used for storing data
        
!        if( row_id .eq. num_row_in_main_matrix + 1 ) then
            
!            !! 6.1. calculate required space
!            numvars = 0 ;
            
!            do  j = row_id , id_max
!                if( index_set(j) .eq. 1 ) then
!                    numvars = numvars + 1 ;
!                end if
!            end do 
            
            

!            !! 6.2. calculate available space in main matrix, put it to beg_id
!            ! capacity = size( solver_dependency ) - solver_dictionary( row_id ) ;
            
!            !! 6.3. if there is enough room, just add this row to main matrix
!            if(size( solver_dependency )  .ge. solver_dictionary( row_id ) + numvars ) then
                
!                k = solver_dictionary( row_id ) ;                       !! 6.3.1. get location of new variables in solver dictionary
                
                
!                do  i = row_id , id_max
!                    if( index_set(i) .eq. 1 ) then
!                        if( temp_row(i) .eq. 0.0_r_kind ) cycle ;
!                        k = k + 1 ;
!                        solver_dependency( k ) = i ;
!                        upper_matrix( k ) = temp_row(i) ;
!                    end if
!                end do 
              
                
                
!                num_row_in_main_matrix = row_id ;                       !! 6.3.2. increment num rows in main matrix
!                solver_dictionary( row_id + 1 )  = k ;                  !! 6.3.3. set next value of dictionary 
                 
                 
!            else  
!                !! 6.4. otherwise, allocated memory for new matrices
!                !! 6.4.1. Calculate required memory for next matrix to be fully filled by other rows
!                !!        It is (N - m )*( N -m + 1 ) / 2 where m = num_row_in_main_matrix
                
!                k=  size( permutation ) - num_row_in_main_matrix  ;
                
!                numvars = 1 + ( k * ( k + 1 ) ) / 2 ;                  !! find required memory
                
!                !! 6.4.2. If auxilary matrix is not allocated just allocate it.
                
!                if( .not. allocated( auxiliary_matrix ) ) then
                
!                    allocate( auxiliary_matrix( int( 1.2 * numvars ) ) , stat = k )   !! rescale memory by an offset
!                    if( k /= 0 ) stop 'Error: failed in allocating cplx_auxiliary_matrix in solver module.' 
                
!                    allocate( auxiliary_dependency( int( 1.2 * numvars ) ) , stat = k )
!                    if( k /= 0 ) stop 'Error: failed in allocating auxiliary_dependency in solver module.' 
                    
!                else  !! 6.4.3. If auxiliary matrices are allocated just check if their size is enough large, as required.
                   
!                    if( size( auxiliary_matrix ) < numvars ) then
                    
!                        deallocate( auxiliary_matrix , stat = k )
!                        deallocate( auxiliary_dependency , stat = k )
!                        if( k /= 0 ) stop ' Error in deallocating auxiliary arrays in solver module'
                        
!                        allocate( auxiliary_matrix( int( 1.2 * numvars ) ) , stat = k )
!                        if( k /= 0 ) stop 'Error: failed in allocating cplx_auxiliary_matrix in solver module.' 
                
!                        allocate( auxiliary_dependency( int( 1.2 * numvars ) ) , stat = k )
!                        if( k /= 0 ) stop 'Error: failed in allocating auxiliary_dependency in solver module.' 
!                    end if
                    
!                end if 
                
!                auxiliary_dependency(:) = 0 
!                auxiliary_matrix(:) =  0.0_r_kind  ;
                
                
!                !! 6.4.4. Store data in auxiliary matrices from their first elements
!                k = 0 ;
             
!                do  i = row_id , id_max
!                    if( index_set(i) .eq. 1 ) then
!                        if( temp_row(i) .eq. 0.0_r_kind ) cycle ;
!                        k = k + 1 ;
!                        auxiliary_dependency( k ) = i ;
!                        auxiliary_matrix( k ) = temp_row(i) ;
!                    end if
!                end do 
                
!                solver_dictionary( row_id + 1 )  = k ;                  !! 8.3.3. set next value of dictionary
                
!            end if 
            
!            index_set( id_min:id_max ) =  0 ;
!            temp_row(  id_min:id_max ) =  0.0_r_kind  ;
            
!            return ;
!        end if
        
        
        
!        !! 6.5. add new data to auxiliary matrices
!        k = solver_dictionary( row_id )  ;
             
!        do  i = row_id , id_max
!            if( index_set(i) .eq. 1 ) then
!                if( temp_row(i) .eq. 0.0_r_kind ) cycle ;
!                k = k + 1 ;
!                auxiliary_dependency( k ) = i ;
!                auxiliary_matrix( k ) = temp_row(i) ;
!            end if
!        end do 
                
!        solver_dictionary( row_id + 1 )  = k ;                          !! 8.3.3. set next value of dictionary
!        index_set( id_min:id_max ) =  0 ;
!        temp_row(  id_min:id_max ) =  0.0_r_kind  ;       
               
               
!    end subroutine  add_real_row_to_matrix
  
  

!!======================================================================
!!
!!   solve_complex_system by direct substitution method
!!   
!!
!!====================================================================== 
   
     
     
subroutine  solve_complex_system( freq )
    use  mod_physics  , only : dictionary, dependency 
    implicit none 
    real( kind = r_kind ), intent ( in ) :: freq  
    
    integer :: i , j
    real( kind = r_kind ) ::  max_error , tmp_re;  
    logical :: is_row_found 

    !! 0. normalize system
    do  i = 1 , size( dictionary ) -1
        tmp_re = 0.0_r_kind  ;
        do  j = dictionary(i)+1 , dictionary(i+1)
            tmp_re = tmp_re + cplx_matrix(j)%re ** 2 + cplx_matrix(j)%im ** 2 ;
        end do
        
        tmp_re =  sqrt( tmp_re ) ;
        
        do  j = dictionary(i)+1 , dictionary(i+1)
            cplx_matrix(j) = cplx_matrix(j) / tmp_re  ;
        end do
        cplx_response(i) = cplx_response(i) / tmp_re ;
    end do
    

    
    !! 1. set initial values
    
    !! permutation is identity matrix, store index of non-zero column at each row
    do i = 1 , size( permutation ) 
       permutation(i) = i ;
    end do
    
    
    num_row_in_main_matrix = 0 ;                                        !! Is the main matrix enough for storing transformed base matrix?
    solver_dependency(:)   = 0 ;                                        !! Initial value of dependency matrix is zero
    cplx_upper_matrix(:)   = cmplx( 0.0_r_kind , 0.0_r_kind ) ;         !! This is also zero vector
    solver_dictionary(:)   = 0 ;                                        !! The same for solver dictionary
    cplx_transformed_response(:)  = cmplx( 0.0_r_kind , 0.0_r_kind ) ;    
   
    
    !! 2. trasnform matrix to an upper triangular one
    !! 2.1. find the k-th row that contains x1 and, moreover, A(k,1) is non-zero.
    
    is_row_found = .false. ;
    
    do  i = 1 , size( dictionary ) - 1                                  !! Loop over all rows
        
        if( dependency( dictionary(i)+1 ) > 1 ) cycle                   !! 2.1.1. If x1 is not a member of this row just ignore the row.
        
        if( abs(cplx_matrix( dictionary(i)+1)) < effective_zero ) cycle !! 2.1.2. check if the matrix A(i,1) is non-zero 
        
        is_row_found = .true. ;                                         !! 2.1.3. a suitable row is found
        permutation(i) = 1 ;                                            !! 2.1.4. permute i-th row with first one
        permutation(1) = i ;                                            !! 2.1.5. permute first row with i-th one
        exit ;                                                          !! 2.1.6. exit the loop over i      
    end do
    
    !! 2.1.7. If there is no such row, the matrix is singular.
    if( .not. is_row_found ) stop 'Error: solver module: matrix is singular: x1 variable cannot be singled out.' 
    
  
    !! 2.1.8. First row may have been replaced by permutation(1)-th row, so the true first row is j = permutation(1).
    !!        Put this row in the first row of upper triangular matrix.
    
   
    solver_dictionary(2) = dictionary( permutation(1) + 1 ) - dictionary( permutation(1) ) ;  !! note that solver_dictionary(1) = 0
    solver_dependency( 1:solver_dictionary(2) ) = dependency(  dictionary( permutation(1) ) + 1  :  &
                                                            &  dictionary( permutation(1) + 1 )) ;
                                                            
    cplx_upper_matrix( 1:solver_dictionary(2) ) = cplx_matrix( dictionary( permutation(1) ) + 1  :  &
                                                            &  dictionary( permutation(1) + 1 )) ;
    
    cplx_transformed_response(1) = cplx_response( permutation(1) ) ;                 !! put response vector to transformed response without permutation
    
    
    num_row_in_main_matrix  = 1 ;                                       !! Is the main matrix enough for storing transformed base matrix? 
    
    
    !! 2.2. add remaining rows to the upper matrix
    
    cplx_temp_row(:) = cmplx( 0.0_r_kind , 0.0_r_kind ) ;
    index_set(:)     = 0 ;
    
    do  i = 2 ,   size( dictionary ) - 1                                 !! 2.2.1. Loop over all remaining rows: from row= 2 to row = n
        call  add_new_row_to_cplx_matrix( i ) ;                          !! 2.2.2. transform and add rows to the upper matrix
    end do
    
    
    cplx_solution(:) = cmplx( 0.0_r_kind , 0.0_r_kind ) ;
    
    do  i = size( permutation ) , num_row_in_main_matrix  + 2 , -1
        cplx_solution(i) = cplx_transformed_response(i) ;
        do  j = solver_dictionary(i)+ 2  , solver_dictionary(i+1)
            cplx_solution(i) = cplx_solution(i) - cplx_auxiliary_matrix(j) * cplx_solution(auxiliary_dependency(j)) ;
        end do 
        
        cplx_solution(i) =  cplx_solution(i)  /  cplx_auxiliary_matrix( solver_dictionary(i) + 1 ) ;
                    
    end do 
    
    !! 3.2. the first item in new database begins from id = 1
    if( size( permutation ) > num_row_in_main_matrix ) then
        i = num_row_in_main_matrix  + 1 ;
        cplx_solution(i) = cplx_transformed_response(i) ;
        do  j = 2  , solver_dictionary(i+1)
            cplx_solution(i) = cplx_solution(i) - cplx_auxiliary_matrix(j) * cplx_solution(auxiliary_dependency(j)) ;
        end do 
        
        cplx_solution(i) =    cplx_solution(i) / cplx_auxiliary_matrix( 1 ) ;
                    
        
    end if 
    
    
    !! 3.2. second for remaining variables
    do i = num_row_in_main_matrix , 1 , -1
    
        cplx_solution(i) = cplx_transformed_response(i) ;
        do  j = solver_dictionary(i)+ 2  , solver_dictionary(i+1)
            cplx_solution(i) = cplx_solution(i) - cplx_upper_matrix(j) * cplx_solution( solver_dependency(j) ) ;
        end do 
        
        cplx_solution(i) =  cplx_solution(i)  / cplx_upper_matrix( solver_dictionary(i) + 1 )  ;
                    
        
    end do



    
    !! 3.3 Testing solution
    cplx_temp_row = cmplx( 0.0_r_kind , 0.0_r_kind ) ;
    index_set(:)  = 0 ;
    max_error  = 0.0_r_kind ;
    
    do i = 1 , size( permutation )
        do j = dictionary(i) + 1 , dictionary(i+1)
           cplx_temp_row(i) = cplx_temp_row(i) + cplx_matrix(j) * cplx_solution( dependency(j) ) ;
        end do
        cplx_temp_row(i) = cplx_temp_row(i) - cplx_response(i) ;
        max_error  = max( max_error  , abs( cplx_temp_row(i) ) ) 
    end do
    
    
    write(*,*) 'freq = ' , freq , ' max error = ' , max_error 
    
     
end subroutine  solve_complex_system
  
  


  
!!=================================================================
!!
!!   for variable x_k (k = row_id), choose j-th rows of matrix A, such that A(k,k) /= 0 after substitution variables x_i, i = 1 , ..., k.
!!
!!================================================================= 

    subroutine  add_new_row_to_cplx_matrix( row_id )  
    
        use  mod_physics , only : dictionary , dependency
        implicit none
        integer , intent( in  ) :: row_id 
        integer :: i, j, k, k_beg , k_end , num_cols, beg_id , end_id , id_min , id_max, first_var
        logical :: found_any  ;
        complex( kind = r_kind ) :: cplx_response_value, cplx_tmp 
        
        
        !! 1. loop over all unused rows that contain x_k, the unused rows are permutation( row_id:end)
        
        found_any = .false. ;
        
        do  i = row_id ,  size( dictionary ) - 1 
        
            !! 2. put i-th row in the index set and matrix row
            
            beg_id = dictionary( permutation( i ) ) + 1 ;
            end_id = dictionary( permutation( i ) + 1 ) ;
        
            !! 3. if the row does not contain x_k variable, just ignore it
            
            if( ( dependency(beg_id) > row_id ) .or. ( dependency(end_id) < row_id ) ) cycle ;  !! simple necessary condition
            

            if( get_location( permutation( i ) , row_id ) .eq. 0 ) cycle ;                          !! sufficient condition for checking if the variable is in the row
            
            id_min = dependency(beg_id) ;
            id_max = dependency(end_id) ;
            
            
            do j = beg_id , end_id
               index_set( dependency(j) ) = 1 ;                         !! these variables are present in this row
               cplx_temp_row( dependency(j) ) = cplx_matrix(j) ;        !! put matrix row in the row-vector
            end do 
            
            first_var = id_min ;
            
            cplx_response_value = cplx_response( permutation( i ) ) ;
            
            
            !! 3. substitute all previously singled out variables, i.e. x_1 from first row, x_2 from second row and so on.
            !! Some variables are in the base matrix and others are stored in the auxiliary matrix, 
            !! So the substitution is separated to two parts, i.e. variables storing in base and auxiliary matrices.
            !! The whole variables in this row are dependency(beg_id:end_id), all variables before diagonal one, i.e. x_k, k= row_id
            !! should be substituted. Therefore, the following loop begins from dependency(beg_id) and runs to dependency( varlocRowId - 1 )
            !! Some are in the base matrix the others are in the auxiliary matrix.
            
            !! substituttion loop 
            !! There are three possibilitieis
            !!      1. if  num_row_in_main_matrix = row_id -1, the base matrix is under considration
            !!      2. if  num_row_in_main_matrix = row_id -2, the auxiliary matrix is chosen for the first time
            !!      3. otherwise the auxiliary matrix is under consideration
            
            !! the substitution loop should be over index set itself as substituting x_1 may introduce x_3 that has not been in the
            !! main index set. The beginning value in the index set is id_min and it runs to row_id -1 in the index set.
            !! if some index was zero then it is not calculated.
            
            
            do  j = first_var , row_id - 1
               
                if( index_set(j) /= 1 ) cycle ;
                 
                if( j .le. num_row_in_main_matrix ) then                !! 3.1. the j-th row is in base matrix
                
                    k_beg = solver_dictionary( j ) + 1 ;
                    k_end = solver_dictionary( j + 1 ) ;
                    
                    id_min = min( id_min , solver_dependency(k_beg) ) ;
                    id_max = max( id_max , solver_dependency(k_end) ) ;
                    
                    
                    cplx_tmp =   cplx_temp_row( j )  / cplx_upper_matrix( k_beg )  ;
                             
                    
                    
                    do k = k_beg , k_end
                        index_set( solver_dependency(k) ) = 1 ;
                        cplx_temp_row(solver_dependency(k)) = cplx_temp_row(solver_dependency(k)) &
                                                          & - cplx_tmp * cplx_upper_matrix( k ) ;
                    end do  
                    
               
                else    !! 3.2. The j-th row is stored in auxiliary matrix
                    
                    !! 3.2.1. In the second case mentioned above, k_beg =1 , i.e. the first item in auxiliary matrix
                    if( j .eq. num_row_in_main_matrix + 1 ) then
                        k_beg = 1 ;
                    else
                        k_beg = solver_dictionary( j ) + 1 ;
                    end if 
                    
                    k_end = solver_dictionary( j + 1 )  ;
                    
                    id_min = min( id_min , auxiliary_dependency(k_beg) ) ;
                    id_max = max( id_max , auxiliary_dependency(k_end) ) ;
                    
                    cplx_tmp = cplx_temp_row( j ) / cplx_auxiliary_matrix( k_beg )  ;
                    
                    
                    do k = k_beg , k_end
                        index_set( auxiliary_dependency(k) ) = 1 ;
                        cplx_temp_row(auxiliary_dependency(k)) = cplx_temp_row(auxiliary_dependency(k)) &
                                                             & - cplx_tmp * cplx_auxiliary_matrix( k ) ;
                    end do   
                end if
            
                !! 3.3. change response value 
                cplx_response_value = cplx_response_value - cplx_tmp * cplx_transformed_response(j) ;
                
            end do !! over j
            

            
            
            !! 4. check if A(k,k) /= 0, otherwise set vectors to default values
            if( abs( cplx_temp_row( row_id ) ) /= 0.0_r_kind ) then
            
                found_any = .true.
                permutation( row_id ) = i ;                                  !! 4.1. do suitable permutaion of rows
                permutation( i ) = row_id ;
                cplx_transformed_response(row_id) = cplx_response_value ;    !! 4.2. add transformed response to the vector
                exit ;                                                       !! 4.3. exit the i-loop
                
            else
            
                index_set( id_min : id_max ) = 0 ;
                cplx_temp_row( id_min:id_max ) = cmplx( 0.0_r_kind , 0.0_r_kind ) ;
            
            end if 
        
        end do ! over i
        
         
        !! 5. check if the variable x_k was singled out
        if( .not. found_any ) then
            write(*,*) 'Error: solver: the variable x_k for k = ' , row_id , ' cannot be singled out, singular matrix' ;
            stop 'Error: solver: Matrix is singular' 
        end if 
       
        
        
        
        !! 6. add new row to upper triangular matrix if there is enough memory
        !!  if num_row_in_main_matrix = RowId - 1 then it is checked whether there is room to add new elements.
        !!  If so, num_row is incremnted and the difference remains the same for next iteration.
        !!  Otherwise, the new auxiliary matrices are allocated and num_row is not changed. 
        !!  For other iterations, RowId - numRow > 1 and the new matrices are used for storing data
        
        if( row_id .eq. num_row_in_main_matrix + 1 ) then
            
            !! 6.1. calculate required space
            num_cols = 0 ;
            do  i = row_id , id_max
                if( index_set(i) .eq. 1 ) then
                    num_cols = num_cols + 1 ;
                end if
            end do 
            

            !! 6.2. calculate available space in main matrix, put it to beg_id
            beg_id = size( solver_dependency ) - solver_dictionary( row_id ) ;
            
            !! 6.3. if there is enough room, just add this row to main matrix
            if( beg_id .ge. num_cols ) then
                
                k = solver_dictionary( row_id ) ;                       !! 6.3.1. get location of new variables in solver dictionary
                
                
                do  i = row_id , id_max
                    if( index_set(i) .eq. 1 ) then
                        k = k + 1 ;
                        solver_dependency( k ) = i ;
                        cplx_upper_matrix( k ) = cplx_temp_row(i) ;
                    end if
                end do 
                
                
                num_row_in_main_matrix = row_id ;                       !! 6.3.2. increment num rows in main matrix
                solver_dictionary( row_id + 1 )  = k ;                  !! 6.3.3. set next value of dictionary
                 
            else  
                !! 6.4. otherwise, allocated memory for new matrices
                !! 6.4.1. Calculate required memory for next matrix to be fully filled by other rows
                !!        It is (N - m )*( N -m + 1 ) / 2 where m = num_row_in_main_matrix
                
                k=  size( permutation ) - num_row_in_main_matrix  ;
                
                num_cols =  ( k * ( k + 1 ) ) / 2 ;                         !! find required memory
                
                !! 6.4.2. If auxilary ariables are not allocated just allocate them.
                
                if( .not. allocated( cplx_auxiliary_matrix ) ) then
                
                    allocate( cplx_auxiliary_matrix( int( 1.2 * num_cols ) ) , stat = k )   !! rescale memory by an offset
                    if( k /= 0 ) stop 'Error: failed in allocating cplx_auxiliary_matrix in solver module.' 
                
                    allocate( auxiliary_dependency( int( 1.2 * num_cols ) ) , stat = k )
                    if( k /= 0 ) stop 'Error: failed in allocating auxiliary_dependency in solver module.' 
                    
                else  !! 6.4.3. If auxiliary matrices are allocated just check if their size is enough large, as required.
                   
                    if( size( cplx_auxiliary_matrix ) < num_cols ) then
                    
                        deallocate( cplx_auxiliary_matrix , stat = k )
                        deallocate( auxiliary_dependency , stat = k )
                        if( k /= 0 ) stop ' Error in deallocating auxiliary arrays in solver module'
                        
                        allocate( cplx_auxiliary_matrix( int( 1.2 * num_cols ) ) , stat = k )
                        if( k /= 0 ) stop 'Error: failed in allocating cplx_auxiliary_matrix in solver module.' 
                
                        allocate( auxiliary_dependency( int( 1.2 * num_cols ) ) , stat = k )
                        if( k /= 0 ) stop 'Error: failed in allocating auxiliary_dependency in solver module.' 
                    end if
                    
                end if 
                
                auxiliary_dependency(:) = 0 
                cplx_auxiliary_matrix(:) = cmplx( 0.0_r_kind , 0.0_r_kind ) ;
                
                
                !! 6.4.4. Store data in auxiliary matrices from their first elements
                k = 0 ;
             
                do  i = row_id , id_max
                    if( index_set(i) .eq. 1 ) then
                        k = k + 1 ;
                        auxiliary_dependency( k ) = i ;
                        cplx_auxiliary_matrix( k ) = cplx_temp_row(i) ;
                    end if
                end do 
                
                solver_dictionary( row_id + 1 )  = k ;                  !! 8.3.3. set next value of dictionary
                
            end if 
            
            index_set( id_min : id_max ) = 0 ;
            cplx_temp_row( id_min:id_max ) = cmplx( 0.0_r_kind , 0.0_r_kind ) ;
            return ;
        end if
        
        
        !! 6.5. add new data to auxiliary matrices
        k = solver_dictionary( row_id )  ;
             
        do  i = row_id , id_max
            if( index_set(i) .eq. 1 ) then
                k = k + 1 ;
                auxiliary_dependency( k ) = i ;
                cplx_auxiliary_matrix( k ) = cplx_temp_row(i) ;
            end if
        end do 
                
        solver_dictionary( row_id + 1 )  = k ;                          !! 8.3.3. set next value of dictionary
        index_set( id_min : id_max ) = 0 ;
        cplx_temp_row( id_min:id_max ) = cmplx( 0.0_r_kind , 0.0_r_kind ) ;        
    
    end subroutine  add_new_row_to_cplx_matrix


!!======================================================================
!!
!!  Is a given number an element of the given row?
!!  
!!  The algorithm uses a binary search in the two first steps and then 
!!  a linear search to find var_loc.
!! 
!!======================================================================
   
    function  get_location( row , col )  result( var_loc )
        use  mod_physics , only : dictionary , dependency
        implicit none
        integer , intent( in  )  :: row , col
        integer :: var_loc
        
     
        integer :: low , high , i
        
        !! 1. set low and high index of the whole set
        
        low  = dictionary(row) + 1   ;
        high = dictionary( row + 1 ) ;
        
        !! 2. do binary search twice to make the set smaller
        do i = 1 , 2
            var_loc  = ( low + high ) / 2
            if( dependency( var_loc ) .eq. col ) return
            
            if( dependency( var_loc )   <  col ) then
                low  = var_loc + 1
            else 
                high = var_loc - 1
            end if
        end do 
   
        !! 3. use linear search on the smaller set of data
        do var_loc = low  , high
           if( dependency( var_loc ) .eq. col ) return
        end do
        
        var_loc = 0 ;
   
    end function get_location
    
    
!!====================================================================!!
!!
!! initialize real solver
!!
!!
!!====================================================================!!

    subroutine  initialize_real_solver()
  
        use  mod_physics  , only : dictionary , dependency, num_max_variables_in_eqns
        use  mod_physics  , only : tension , compression, element_property, analysis_kind
        implicit none
        
        integer :: state , i , num_elastic_gauss_point = 0 ;
     
        do  i = 1 , size( element_property , 1)
            if( element_property(i,1) < 3 ) then
                num_elastic_gauss_point = num_elastic_gauss_point + element_property(i,3)**2 ;
            end if
        end do
        
        
     
        
        call  clean_up_solver()
     
        !! 1. allocated memory for main vectors
     
        allocate( matrix( size( dependency ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating matrix in solver module.' 
      
        allocate( response( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating response vector in solver module.' 
      
        allocate( solution( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solution vector in solver module.' 
     
        !! 2. Allocate memory for arrays used in direct substitution solver
     
        allocate( index_set( size( dictionary )  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating index_set  in solver module.' 
       
       
        allocate( permutation( size( dictionary ) -1   ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating permutation in solver module.' 
       
     
     
     
        allocate( temp_row( size( dictionary ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating temp_row  in solver module.'
   
        allocate( solver_dictionary( size( dictionary ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solver_dictionary  in solver module.' 
        
        !! initialize solver matrix with size 10 * size ( dependency )
        
        allocate( solver_dependency( 10 * size( dependency ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solver_dependency  in solver module.' 
        
        allocate( upper_matrix( 10 * size( dependency ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solver_dependency  in solver module.' 
        
        allocate( transformed_response( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating transformed_response  in solver module.'
        
        allocate( u_vec( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating u_vec  in solver module.'
        
        allocate( udot_vec( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating udot_vec in solver module.'
        
        
        allocate( uddot_vec( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating uddot_vec  in solver module.'
        
        allocate( temp_real_vec( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating temp_real_vec  in solver module.'
        
        allocate( static_solution( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating static_solution  in solver module.'
    
        if( analysis_kind .eq. 2) then
            state = 0 ;
            if( allocated(tension) ) deallocate(tension , stat = state )
            if( state /= 0 ) stop 'Error: failed in deallocating tension  in solver module.'
        
            if( allocated(compression ) ) deallocate(compression , stat = state )
            if( state /= 0 ) stop 'Error: failed in deallocating compression  in solver module.'
        
            allocate( tension( num_elastic_gauss_point ,9  ) , stat = state )
            if( state /= 0 ) stop 'Error: failed in allocating tension  in solver module.'
        
            allocate( compression( num_elastic_gauss_point ,9  ) , stat = state )
            if( state /= 0 ) stop 'Error: failed in allocating compression  in solver module.'
        end if 
    
    end subroutine initialize_real_solver
    
    
    
!!======================================================================
!!
!!   initialize_complex_solver
!!
!!====================================================================== 

   
   
    subroutine  initialize_complex_solver()
  
        use  mod_physics  , only : dictionary , dependency, num_max_variables_in_eqns
        implicit none
        
        integer :: state
     
        call  clean_up_solver()
     
        !! 1. allocated memory for main vectors
     
        allocate( cplx_matrix( size( dependency ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating complex matrix in solver module.' 
      
        allocate( cplx_response( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating complex response vector in solver module.' 
      
        allocate( cplx_solution( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating complex solution vector in solver module.' 
     
        !! 2. Allocate memory for arrays used in direct substitution solver
     
        allocate( index_set( size( dictionary )  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating index_set  in solver module.' 
       
       
        allocate( permutation( size( dictionary ) -1   ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating permutation in solver module.' 
        
       
     
        allocate( cplx_temp_row( size( dictionary )  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating cplx_temp_row  in solver module.'
   
        allocate( solver_dictionary( size( dictionary ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solver_dictionary  in solver module.' 
        
        !! initialize solver matrix with size 2 * size ( dependency )
        
        allocate( solver_dependency( 10 * size( dependency ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solver_dependency  in solver module.' 
        
        allocate( cplx_upper_matrix( 10 * size( dependency ) ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating solver_dependency  in solver module.' 
        
        allocate( cplx_transformed_response( size( dictionary ) -1  ) , stat = state )
        if( state /= 0 ) stop 'Error: failed in allocating cplx_transformed_response  in solver module.'
        
   
    end subroutine initialize_complex_solver 
     
     
    
 
!!=================================================================
!!
!!   clean_up_solver
!!
!!================================================================= 
     
  subroutine  clean_up_solver()
    use mod_physics, only : analysis_kind
    implicit none
    integer :: state , state1 = 0
    
 
    if( analysis_kind .eq. 10 ) then
        if( allocated(  cplx_matrix ) ) then
            deallocate( cplx_matrix , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  cplx_response ) ) then
            deallocate( cplx_response , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  cplx_solution ) ) then
            deallocate( cplx_solution , stat = state )
            state1 = state1 + state 
        end if
        
        if( allocated(  cplx_upper_matrix ) ) then
            deallocate( cplx_upper_matrix , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  cplx_auxiliary_matrix ) ) then
            deallocate( cplx_auxiliary_matrix , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  cplx_temp_row  ) ) then
            deallocate( cplx_temp_row  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  cplx_transformed_response  ) ) then
            deallocate( cplx_transformed_response  , stat = state )
            state1 = state1 + state 
        end if
    end if
    
     
    if( allocated(  solver_dictionary ) ) then
        deallocate( solver_dictionary , stat = state )
        state1 = state1 + state 
    end if
    
    
    if( allocated(  index_set ) ) then
        deallocate( index_set , stat = state )
        state1 = state1 + state 
    end if
    
    if( allocated(  solver_dependency ) ) then
        deallocate( solver_dependency , stat = state )
        state1 = state1 + state 
    end if
    
    
    if( allocated(  auxiliary_dependency ) ) then
        deallocate( auxiliary_dependency , stat = state )
        state1 = state1 + state 
    end if
    
    if( allocated(  permutation ) ) then
        deallocate( permutation , stat = state )
        state1 = state1 + state 
    end if
    
   
    if( analysis_kind .eq. 2 ) then
        if( allocated(  matrix  ) ) then
            deallocate( matrix  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  response  ) ) then
            deallocate( response  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  solution  ) ) then
            deallocate( solution  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  upper_matrix  ) ) then
            deallocate( upper_matrix  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  auxiliary_matrix  ) ) then
            deallocate( auxiliary_matrix  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  temp_row  ) ) then
            deallocate( temp_row  , stat = state )
            state1 = state1 + state 
        end if
    
        if( allocated(  transformed_response ) ) then
            deallocate( transformed_response  , stat = state )
            state1 = state1 + state 
        end if
        
        if( allocated(  u_vec ) ) then
            deallocate( u_vec  , stat = state )
            state1 = state1 + state 
        end if
        
        if( allocated(  udot_vec ) ) then
            deallocate( udot_vec  , stat = state )
            state1 = state1 + state 
        end if
        
        if( allocated(  uddot_vec ) ) then
            deallocate( uddot_vec  , stat = state )
            state1 = state1 + state 
        end if
        
        if( allocated(  temp_real_vec ) ) then
            deallocate( temp_real_vec  , stat = state )
            state1 = state1 + state 
        end if
        
        if( allocated(  static_solution ) ) then
            deallocate( static_solution  , stat = state )
            state1 = state1 + state 
        end if
        
    end if
    
    if( state1 /= 0 ) stop ' Error in deallocating arrays in solver module'
  
  end subroutine clean_up_solver   
     
     

     
 
  
 end module mod_solver
