Generate random points in polygon. coded in fortran.

Compile:
gfortran poltrian.f90 list_final.f90

Run:
./a.out
Random points for 3 polygons will be written in ./POINTS dir.

Plot:
gnuplot < i18_commands.txt;gnuplot < hand_commands.txt;gnuplot < comb_commands.txt
will plot png file for each polygon with random points inside.


based on:
https://people.math.sc.edu/Burkardt/f_src/polygon_triangulate/polygon_triangulate.html

https://blogs.sas.com/content/iml/2020/10/21/random-points-in-polygon.html#:~:text=The%20first%20step%20is%20to,the%20union%20of%20the%20triangles.

https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html
