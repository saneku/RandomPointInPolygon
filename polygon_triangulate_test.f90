program main
  implicit none
  !gfortran poltrian.f90 polygon_triangulate_test.f90
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: number_of_random_points = 20000

  integer ( kind = 4 ) triangles(3,n-2)

  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    8.0D+00, 7.0D+00, 6.0D+00, 5.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 1.0D+00, 0.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: y = (/ &
    0.0D+00, 10.0D+00,  0.0D+00, 10.0D+00,  0.0D+00, &
   10.0D+00,  0.0D+00, 10.0D+00,  0.0D+00, -2.0D+00 /)

   real ( kind = 8 ), dimension ( number_of_random_points ) :: xr,yr

   integer ( kind = 4 ) i,t,number_of_points_per_triangle
   real ( kind = 8 ) polygon_area,triangle_area
   real ( kind = 8 ) whole_area,tr_area

   !real ( kind = 8 ), dimension ( n-2 ) :: triangles_areas
   !integer ( kind = 4 ):: triangles(3,8) = transpose(reshape((/3,5,7,9,10,3,5,7,1,3,5,7,7,10,10,10,2,4,6,8,9,1,3,5/), (/8,3/)))
   
  !call my_polygon_triangulate ( n, x, y, triangles,triangles_areas )
  !print *, triangles_areas
  !print *, sum(triangles_areas)
   triangles=0
   xr=0
   yr=0


   call polygon_triangulate ( n, x, y, triangles )

  call i4mat_transpose_print ( 3, n - 2, triangles, '  Triangles' )

  whole_area = polygon_area ( n, x, y )
  print *,whole_area

  t=1
  do while ( t <= n-2 )
      !print *, t

      !number_of_points_per_triangle=triangles_areas(t)/sum(triangles_areas)

      !print *,x(triangles(1,t)),y(triangles(1,t))
      !print *,x(triangles(2,t)),y(triangles(2,t))
      !print *,x(triangles(3,t)),y(triangles(3,t))

      tr_area = triangle_area ( x(triangles(1,t)),y(triangles(1,t)), &
                                  x(triangles(2,t)),y(triangles(2,t)), &
                                  x(triangles(3,t)),y(triangles(3,t)))
      print *,tr_area

      number_of_points_per_triangle = (number_of_random_points*(tr_area/whole_area))
      print *,number_of_points_per_triangle

      call randUnifTriangle(x(triangles(1,t)),y(triangles(1,t)), &
                            x(triangles(2,t)),y(triangles(2,t)), &
                            x(triangles(3,t)),y(triangles(3,t)), &
                            number_of_points_per_triangle,xr,yr)

      i=0
      do while ( i < number_of_points_per_triangle )
        !print *,xr(i),yr(i)
        i=i+1
      end do 
      t=t+1

  end do 

  stop 0
end