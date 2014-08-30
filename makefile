all:
	gcc newton.c -o newton_fractal -O3 -fomit-frame-pointer -lm -w
	gcc polynom_root_finder.c -o root_finder -O3 -lm -fomit-frame-pointer -w

clean:
	-rm -rf core* a.out newton_fractal feh* *.dat *.ppm root_finder

icc:
	icc newton.c -o newton_fractal -O3 -fomit-frame-pointer -lm -ipo -ip -parallel
	icc polynom_root_finder.c -o root_finder -O3 -fomit-frame-pointer -lm -ipo -ip -parallel

help:
	echo "Use make for compiling the newthon-raphson fractal"
	echo "make icc compiles it with the intel compiler, the program gets WAY faster"
	echo "make clean... well... clean the directory."
