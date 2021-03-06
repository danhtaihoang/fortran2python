import numpy as np
from subprocess import call

nx=4
ny=2
n=nx*ny

w = np.empty(n) ; h = np.empty(n)
x = np.empty(n) ; y = np.empty(n)

w[0]=0.16   
h[0]=w[0]/1.1
if n>1:
   w[1:n]=w[0]
   h[1:n]=h[0]

for i in range(nx):
   x[i]=0.05+i*(w[0]+0.06)
   for j in range(1,ny):
      x[j*nx+i]=x[i]

for j in range(ny):
   y[j*nx:j*nx+nx]=1.-(j+1)*(h[0]+0.06)

x2=x+w
y2=y+h

#===============================================================================
f = open('comment.dat', 'w')

f.write("set terminal postscript eps size 20cm,20cm enhanced color font 'Helvetica,10'\n")
f.write("set autoscale\n")

f.write("set style line 1 lt 1 pt 7 lw 2 lc rgb 'black'\n")
f.write("set style line 2 lt 2 pt 5 lw 2 lc rgb 'blue'\n")
f.write("set style line 3 lt 3 pt 13 lw 2 lc rgb 'magenta'\n")
f.write("set style line 4 lt 1 pt 9 lw 2 lc rgb 'red'\n")
f.write("set style line 5 lt 1 pt 11 lw 2 lc rgb 'orange'\n")

f.write("set style line 6 lt 1 pt 65 lw 2 lc rgb 'black'\n")
f.write("set style line 7 lt 2 pt 64 lw 2 lc rgb 'blue'\n")
f.write("set style line 8 lt 3 pt 68 lw 2 lc rgb 'magenta'\n")
f.write("set style line 9 lt 1 pt 66 lw 2 lc rgb 'red'\n")
f.write("set style line 10 lt 1 pt 67 lw 2 lc rgb 'orange'\n")

f.write("set style line 11 lt 1 pt 7 dt 3 lw 2 lc rgb 'black'\n")
f.write("set style line 12 lt 2 pt 5 dt 3 lw 2 lc rgb 'blue'\n")
f.write("set style line 13 lt 3 pt 13 dt 3 lw 2 lc rgb 'magenta'\n")
f.write("set style line 14 lt 1 pt 9 dt 3 lw 2 lc rgb 'red'\n")
f.write("set style line 15 lt 1 pt 11 dt 3 lw 2 lc rgb 'orange'\n")

f.write("set style line 16 lt 1 pt 54 dt 3 lw 2 lc rgb 'black'\n")
f.write("set style line 17 lt 2 pt 64 dt 3 lw 2 lc rgb 'blue'\n")
f.write("set style line 18 lt 3 pt 68 dt 3 lw 2 lc rgb 'magenta'\n")
f.write("set style line 19 lt 1 pt 66 dt 3 lw 2 lc rgb 'red'\n")
f.write("set style line 20 lt 1 pt 67 dt 3 lw 2 lc rgb 'orange'\n")

#--------------------------------------------------------------------------------
f.write("set output 'W-W0.eps'\n")
f.write("set multiplot\n")  
f.write("set pointsize 0.3\n")

f.write("set xlabel 'W0' offset 0,0\n")
f.write("set ylabel 'W' offset 2,0\n")
#f.write("set xrange [0.5:20.5]\n")
#f.write("set yrange [0.5:20.5]\n")
#f.write("set cbrange [-1.0:1.0]\n")
#f.write("set view map\n")

for i in range(n):
   f.write("set lmargin at screen %.5f\n" %x[i])
   f.write("set rmargin at screen %.5f\n" %x2[i])
   f.write("set bmargin at screen %.5f\n" %y[i])
   f.write("set tmargin at screen %.5f\n" %y2[i])

#------------------------------------------------------------------------------
   #f.write("set title 'figure %2d'\n" %i)
   #f.write("plot 'H0.dat' u 1:2 w lp notitle dt 2 pt 7 lt 1 lc rgb 'black' \n")

   if i==0:
      f.write("plot 'Wtest.dat' u 3:4 w p notitle ls 1,\\\n")
      f.write("x w l notitle ls 4\n")
      #f.write("2*x-1 w l notitle ls 3\n")

#------------------------------------------------------------------------------
f.write("unset multiplot\n") 
f.close()
call("./1plot_gnuplot_python.script")


