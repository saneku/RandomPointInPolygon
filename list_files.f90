program read_files_in_folder
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int32

    !gfortran poltrian.f90 list_files.f90

    !gnuplot < i18_commands.txt;gnuplot < hand_commands.txt;gnuplot < comb_commands.txt

    implicit none

    ! Declarations
    character(len=256) :: filename
    integer :: i, ios
    integer, parameter :: max_files = 100
    integer ( kind = 4 ), parameter :: number_of_points_per_poly = 2000
    character(len=256), dimension(max_files) :: file_list
    integer :: num_files


    integer ( kind = 4 ) n,t,number_of_points_per_triangle
    integer ( kind = 4 ) dim_num,offset
    integer ( kind = 4 ) triangle_num
    integer ( kind = 4 ), allocatable :: triangles(:,:)
    real ( kind = 8 ), allocatable :: xy(:,:)
    real ( kind = 8 ) polygon_area,triangle_area
    real ( kind = 8 ) poly_whole_area,tr_area
    !real ( kind = 8 ), dimension ( number_of_points_per_poly ) :: xr,yr
    real ( kind = 8 ), allocatable :: rp(:,:)

    ! Step 1: Get the list of files using a system command
    call execute_command_line("ls -1 POLYGONS > file_list.txt", wait=.true.)
    ! Step 2: Read the list of files into an array
    open(unit=10, file='file_list.txt', status='old', action='read')
    num_files = 0
    do
        read(10, '(A)', iostat=ios) filename
        if (ios /= 0) exit
        num_files = num_files + 1
        file_list(num_files) = trim(filename)
    end do
    close(10)
    call execute_command_line("rm file_list.txt", wait=.true.)

    ! Step 3: Read the content of each file
    do i = 1, num_files
        print *, 'Reading file:', trim(file_list(i))

        call r8mat_header_read ( "POLYGONS/" // file_list(i), dim_num, n )
        allocate ( xy(1:2,1:n) )
        call r8mat_data_read ( "POLYGONS/" // file_list(i), 2, n, xy )
        
        poly_whole_area = polygon_area ( n, xy(1,:), xy(2,:) )
        print *,poly_whole_area

        triangle_num = n - 2
        
        allocate ( triangles(1:3,1:triangle_num) )
        call polygon_triangulate ( n, xy(1,:), xy(2,:), triangles )
        
        ! triangles from polygon
        call i4mat_transpose_print ( 3, triangle_num, triangles, '  Triangles' )

        ! Triangle loop
        allocate ( rp(2,number_of_points_per_poly))
        t=1     
        offset = 1  ! Starting offset   
        do while ( t <= triangle_num )
            tr_area = triangle_area ( xy(1,triangles(1,t)),xy(2,triangles(1,t)), &
                                      xy(1,triangles(2,t)),xy(2,triangles(2,t)), &
                                      xy(1,triangles(3,t)),xy(2,triangles(3,t)))

            number_of_points_per_triangle = int(number_of_points_per_poly*(tr_area/poly_whole_area))
            print *,tr_area,number_of_points_per_triangle

            call randUnifTriangle(xy(1,triangles(1,t)),xy(2,triangles(1,t)), &
                                  xy(1,triangles(2,t)),xy(2,triangles(2,t)), &
                                  xy(1,triangles(3,t)),xy(2,triangles(3,t)), &
                                    number_of_points_per_triangle,rp,offset)
            offset = offset + number_of_points_per_triangle
            t=t+1
        end do 
        print *,offset
        call my_r8mat_write (  "POINTS/" // file_list(i), 2, offset, rp )

        print *, '----------------'
        
        deallocate ( rp )
        deallocate ( triangles )
        deallocate ( xy )
    end do


    contains

    subroutine randUnifTriangle ( x1, y1, x2, y2, x3, y3,ln,lrp,loffset)
        implicit none
      
        integer ( kind = 4 ) ln !number_of_random_points
        real ( kind = 8 ), intent(inout) :: lrp(:,:)
        
        real ( kind = 8 ) x1,x2,x3
        real ( kind = 8 ) y1,y2,y3
        real ( kind = 8 ) ax,ay,bx,by
        integer ( kind = 4 ) li,loffset
        real :: u1,u2
    
        !a and b are vectors at the origin
        ax = x2-x1
        ay = y2-y1
    
        bx = x3-x1
        by = y3-y1
    
        li=0
        do while ( li <= ln )
          call random_number(u1)
          call random_number(u2)
    
          if (u1+u2 >= 1.0 ) then
            u1 = 1 - u1
            u2 = 1 - u2
          end if
          
          lrp(1,li+loffset) = u1*ax + u2*bx+x1
          lrp(2,li+loffset) = u1*ay + u2*by+y1
    
          li = li+1
        end do 
    
        !return
      end subroutine randUnifTriangle


end program read_files_in_folder
