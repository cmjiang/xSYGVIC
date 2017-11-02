MAIN_DIR = ..

include $(MAIN_DIR)/make.inc

SRC = $(wildcard *.f90)
OBJ = $(patsubst %.f90,%.o,$(SRC))

lib: 	$(OBJ)
	$(ARCH) $(ARFLAG) $(MAIN_DIR)/$(GSEP_LIB) $(OBJ)
	$(RANLIB) $(MAIN_DIR)/$(GSEP_LIB)

%.o:    %.f90
	$(FC) $(FC_FLAGS) -c $< 

clean:
	rm -f ../libgsep.a *.o 
