program read_files_in_folder
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int32

    !gfortran poltrian.f90 list_final.f90

    !gnuplot < i18_commands.txt;gnuplot < hand_commands.txt;gnuplot < comb_commands.txt

    implicit none

    type :: triangle
    real ( kind = 8 ) :: x1, y1
    real ( kind = 8 ) :: x2, y2
    real ( kind = 8 ) :: x3, y3
    end type triangle


    ! Declarations
    integer :: i,t,num_polys

    integer, parameter :: max_number_of_polys = 10
    integer, parameter :: max_number_of_triangles_per_poly = 100    
    real ( kind = 8 ), dimension(max_number_of_polys) :: polys_area
    integer ( kind = 4 ), dimension(max_number_of_polys) :: tris_counter
    type(triangle), dimension(max_number_of_polys,max_number_of_triangles_per_poly) :: tris
    real ( kind = 8 ), dimension(max_number_of_polys,max_number_of_triangles_per_poly) :: tris_areas


    !integer ( kind = 4 ), parameter :: number_of_points_per_poly = 2000
    !integer ( kind = 4 ) number_of_points_per_triangle
    !integer ( kind = 4 ) triangle_num
    !real ( kind = 8 ) tr_area
    !real ( kind = 8 ), allocatable :: rp(:,:)

    

    call initPolygons(num_polys)
  
    
    do i = 1, num_polys
        do t = 1, tris_counter(i)
          print *,tris_areas(i,t)
        end do 
        print *,"\n"
    end do


!allocate ( rp(2,number_of_points_per_poly))
!offset = 1  ! Starting offset   
!number_of_points_per_triangle = int(number_of_points_per_poly*(tr_area/polys_area(i)))
!print *,tr_area,number_of_points_per_triangle

!call randUnifTriangle(xy(1,triangles(1,t)),xy(2,triangles(1,t)), &
!                      xy(1,triangles(2,t)),xy(2,triangles(2,t)), &
!                      xy(1,triangles(3,t)),xy(2,triangles(3,t)), &
!                        number_of_points_per_triangle,rp,offset)
!offset = offset + number_of_points_per_triangle
!print *,offset
!call my_r8mat_write (  "POINTS/" // poly_list(i), 2, offset, rp )
!    deallocate ( rp )




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

    subroutine getNumberOfPolygons(num_files,file_list,max_number_of_files)
        implicit none
        integer, intent(inout) :: num_files
        integer, intent(in) ::max_number_of_files
        character(len=256), intent(inout) :: file_list(:)

        character(len=256) :: filename
        integer :: ios

        ! Step 1: Get the list of files using a system command
        call execute_command_line("ls -1 POLYGONS > file_list.txt", wait=.true.)
        ! Step 2: Read the list of files into an array
        open(unit=10, file='file_list.txt', status='old', action='read')
        num_files = 0
        do
            read(10, '(A)', iostat=ios) filename
            if (ios /= 0) exit
            num_files = num_files + 1
            
            if (num_files>max_number_of_files) then
                write ( *, '(a)' ) 'Error: number of polygon files exceeds max_number_of_polys. Exiting...'
                stop 1
            end if
            file_list(num_files) = trim(filename)
        end do
        close(10)
        ! Step 3: Remove file_list.txt from the disk
        call execute_command_line("rm file_list.txt", wait=.true.)
    end subroutine getNumberOfPolygons


    subroutine initPolygons(lnum_polys)
      implicit none
      integer, intent(in):: lnum_polys
      integer :: li,lt
      real ( kind = 8 ) polygon_area,triangle_area
      integer ( kind = 4 ) dim_num!,offset
      integer ( kind = 4 ) num_nodes
      character(len=256), dimension(max_number_of_polys) :: poly_list
      
      integer ( kind = 4 ), allocatable :: triangles(:,:)
      real ( kind = 8 ), allocatable :: xy(:,:)
      
      call getNumberOfPolygons(num_polys,poly_list,max_number_of_polys)
      print *, 'Found ',num_polys, ' polygons'
      do li = 1, lnum_polys
        print *, 'Reading file:', trim(poly_list(li))

        call read_file_header ( "POLYGONS/" // poly_list(li), dim_num, num_nodes )
        allocate ( xy(2,num_nodes) )
        call read_data_from_file ( "POLYGONS/" // poly_list(li), 2, num_nodes, xy )
        
        polys_area(li) = polygon_area ( num_nodes, xy(1,:), xy(2,:) )

        tris_counter(li) = num_nodes - 2
        if (tris_counter(li)>max_number_of_triangles_per_poly) then
          write ( *, '(a)' ) 'Error: number of triangles exceeds max_number_of_triangles_per_poly. Exiting...'
          stop 1
        end if

        allocate ( triangles(3,tris_counter(li)) )
        call polygon_triangulate ( num_nodes, xy(1,:), xy(2,:), triangles )
        ! print triangles from polygon
        call i4mat_transpose_print ( 3, tris_counter(li), triangles, '  Triangles' )

        ! loop over triangles in one poly
        do lt = 1, tris_counter(li)
            ! store the triangle vertices
            tris(li,lt)%x1 = xy(1,triangles(1,lt)); tris(li,lt)%y1 = xy(2,triangles(1,lt))
            tris(li,lt)%x2 = xy(1,triangles(2,lt)); tris(li,lt)%y2 = xy(2,triangles(2,lt))
            tris(li,lt)%x3 = xy(1,triangles(3,lt)); tris(li,lt)%y3 = xy(2,triangles(3,lt))

            !tris_areas(i,t) = triangle_area ( xy(1,triangles(1,t)),xy(2,triangles(1,t)), &
            !                                  xy(1,triangles(2,t)),xy(2,triangles(2,t)), &
             !                                 xy(1,triangles(3,t)),xy(2,triangles(3,t)))

            tris_areas(li,lt) = triangle_area ( tris(li,lt)%x1,tris(li,lt)%y1, &
                                              tris(li,lt)%x2,tris(li,lt)%y2, &
                                              tris(li,lt)%x3,tris(li,lt)%y3)
        end do 
        print *, '----------------'
        
        deallocate ( triangles )
        deallocate ( xy )
    end do

    end subroutine initPolygons
end program read_files_in_folder
