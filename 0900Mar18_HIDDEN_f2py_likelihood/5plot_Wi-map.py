import numpy as np
from subprocess import call

nx=3
ny=2
n=nx*ny

w = np.empty(n) ; h = np.empty(n)
x = np.empty(n) ; y = np.empty(n)

w[0]=0.16   
h[0]=w[0]/1.1
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
f.write("set output 'W.eps'\n")
f.write("set multiplot\n")  
f.write("set autoscale\n")
f.write("set pointsize 0.5\n")   
f.write("set xlabel 'j' offset 0,0\n")
f.write("set ylabel 'i' offset 2,0\n")
f.write("set xrange [0.5:20.5]\n")
f.write("set yrange [0.5:20.5]\n")
f.write("set cbrange [-1.0:1.0]\n")
f.write("set view map\n")

for i in range(n):
   f.write("set lmargin at screen %.5f \n" %x[i])
   f.write("set rmargin at screen %.5f \n" %x2[i])
   f.write("set bmargin at screen %.5f \n" %y[i])
   f.write("set tmargin at screen %.5f \n" %y2[i])

#------------------------------------------------------------------------------
   f.write("set title 'figure %2d'\n" %i)
   f.write("splot 'W.dat' u 1:2:3 notitle w image\n") 

#------------------------------------------------------------------------------
f.write("unset multiplot\n") 
f.close()
call("./1plot_gnuplot_python.script")

!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:4 w l t "CAD" lt 1 lc rgb "blue",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:4 w l t "MXN" lt 1 lc rgb "magenta",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:5 w l t "CHF" lt 1 lc rgb "red"'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:6 w l t "NOK" lt 1 lc rgb "pink"'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:7 w l t "SEK" lt 1 lc rgb "red",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:8 w l t "GBP" lt 1 lc rgb "red",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:9 w l t "AUD" lt 1 lc rgb "red",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:10 w l t "NZD" lt 1 lc rgb "red",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:11 w l t "JPY" lt 1 lc rgb "red",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:12 w l t "SGD" lt 1 lc rgb "yellow",\'
         !WRITE(31,'(a)')'"currency_normalized.dat" u 1:13 w l t "KRW" lt 1 lc rgb "red"'         

  !WRITE(31,'(a)')'set cbrange [-0.15:0.15]'
   !WRITE(31,'(a)')'set logscale xy 10'
   !WRITE(31,'(a)')"set format y '$%g$'"
   !WRITE(31,'(a)')"set format xy '$10^{%T}$'"
   
   !WRITE(31,'(a)',advance="no")'set xtics ("EUR"1,"JPY"2,"HKD"3,"USD"4,"CAD"5,"MXN"6,"CHF"7,"NOK"8,"SEK"9,"GBP"10,"AUD"11,'
   !WRITE(31,'(a)')'"NZD"12,"SGD"13,"KRW"14,"BRL"15)'
   !WRITE(31,'(a)',advance="no")'set ytics ("EUR"1,"JPY"2,"HKD"3,"USD"4,"CAD"5,"MXN"6,"CHF"7,"NOK"8,"SEK"9,"GBP"10,"AUD"11,'
   !WRITE(31,'(a)')'"NZD"12,"SGD"13,"KRW"14,"BRL"15)'
       
   !WRITE(31,'(a)')'set xrange [0.05:0.45]'
   !WRITE(31,'(a)')'set yrange [0.35:0.45]'
   !WRITE(31,'(a)')'set ytics 0.05' 
   !WRITE(31,'(a)')'set pointsize 1.2'
   !WRITE(31,'(a)')"set xlabel 'C-true'"
   !WRITE(31,'(a)')'set xlabel ""'  !!no tic labels 
   !WRITE(31,'(a)')'unset xtics'
   !WRITE(31,'(a)')"set ylabel 'C'"
   !WRITE(31,'(a)')'set mxtics 1'
   !WRITE(31,'(a)')'set xtics 0.1'
   !WRITE(31,'(a)')"set xtics font 'Arial,10'"

      !IF (ifig==(ny+1)) THEN
         !WRITE(31,'(a)')'unset xlabel'
         !WRITE(31,'(a)')'unset ylabel'         
         !WRITE(31,'(a)')'unset logscale xy'
         !WRITE(31,'(a)')'set logscale x 10'
        !WRITE(31,'(a)')'set title "slope"'
         !WRITE(31,'(a)')'set yrange [0.0:2.0]'
         !WRITE(31,'(a)')"set format y '$%g$'"
      !END IF

  !WRITE(31,'(a)')"set xlabel 'time'"
  !WRITE(31,'(a)')"set ylabel 'currency' offset 2,0"
  !WRITE(31,'(a)')'plot "currency_normalized.dat" u 2:3 w l notitle lt 1 lc rgb "black"' 

   !WRITE(31,'(a)')"set cbtics 0.2 offset -0.5,0" 

   !CALL system('gnuplot44 -persist comment.dat')

!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


