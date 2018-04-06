sldg.exe: main.f90 LU.o element_mod.o globals2d.o
	gfortran main.f90 LU.o element_mod.o  globals2d.o -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -O3 -o sldg.exe
	rm *.o
LU.o: LU.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none LU.f90
element_mod.o:element_mod.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none element_mod.f90
globals2d.o:globals2d.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none globals2d.f90
clean:
	rm sldg.exe