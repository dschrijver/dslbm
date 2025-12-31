COMPILER = mpicc
OPT_FLAGS = -std=c11 -O3 
DEBUG_FLAGS = -Wall -Wextra
LIB_FLAGS = -I$(shell pwd)/hdf5/include -L$(shell pwd)/hdf5/lib -l:libhdf5.a -lm -lz -lsz

SRC = $(wildcard main.c src/*.c)

all: cleandata clean dslbm
	mpirun -n 6 dslbm

dslbm:
	$(COMPILER) $(SRC) -o $@ $(OPT_FLAGS) $(DEBUG_FLAGS) $(LIB_FLAGS)

clean:
	rm -f dslbm

cleandata:
	rm -f *.h5

install_hdf5:
	mkdir -p hdf5
	tar -xvzf archives/hdf5* -C hdf5 --strip-components=2
	cd hdf5;\
	CC=mpicc ./configure --enable-parallel --enable-shared;\
	make;\
	make install;\
	cd ../;\
	mv hdf5/hdf5/include hdf5/;\
	mv hdf5/hdf5/lib hdf5/