set terminal pdf enhanced font "Times,24" size 12in, 8in
set output "out_x.pdf"
set title "" font "Times,32"
set xlabel "nstep" font "Times,32"
set ylabel "Disp" font "Times,32"
set key top right
plot "outtest.dat" using($1):($2) with lines lw 10   



set terminal pdf enhanced font "Times,24" size 12in, 8in
set output "out_y.pdf"
set title "1e6 " font "Times,32"
set xlabel "nstep" font "Times,32"
set ylabel "Disp" font "Times,32"
set key top right
plot "outtest.dat" using($1):($3) with lines lw 10   


#plot "ban.csv" using($1<50e-3 ? $1/1e-3:1/0):(-$3) with lines lw 1 lt -1, \
#       "ban.csv" using($1/1e-3):(-$2) with lines lw 1 lt 1
