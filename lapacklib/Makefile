MAIN_DIR = ..

include $(MAIN_DIR)/make.inc

SRC = $(wildcard *.f)
OBJ = $(patsubst %.f,%.o,$(SRC))

lib: 	$(OBJ)
	$(ARCH) $(ARFLAG) $(LAPACK_LIB) $(OBJ)
	$(RANLIB) $(LAPACK_LIB)

%.o:    %.f
	$(FC) $(FC_FLAGS) -c $< 

clean:
	rm -f $(MAIN_DIR)/liblapack.a *.o 
