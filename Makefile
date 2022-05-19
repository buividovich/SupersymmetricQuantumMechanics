SRC = ./linalg.cpp ./susyqm.cpp

HDR = $(SRC:.cpp=.hpp)
HDR += ./timing.hpp ./ansi_io.hpp ./D4d.hpp ./C4v.hpp
HDR += ./arpackpp/*.h  

CC = g++ --std=c++14 -O2 -I./ -L./ -fopenmp -fmax-errors=1 

CC += -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

LIB = -lm -larpack -lopenblas -lgfortran -lm

#Some platform-specific libraries
LIB +=  -lboost_program_options-mt
CC +=  -I/mingw64/include/OpenBLAS/

otoc: otoc.cpp $(SRC) $(HDR)
	$(CC) ./$< $(SRC) $(LIB) -o ./$@
	
otoc_light: otoc_light.cpp $(SRC) $(HDR)
	$(CC) ./$< $(SRC) $(LIB) -o ./$@
	
clean:
	rm -f -v ./otoc