MAIN_DIR = ..
include $(MAIN_DIR)/make.inc

OBJ = dgvres.o dgvbor.o dormm.o dprtm10.o

LIB = $(MAIN_DIR)/$(GSEP_LIB) $(LAPACK_LIB)

all: tstdsygvic

tstdsygvic: $(OBJ)

%: %.o
	$(FC) $< $(OBJ) -o $@ $(LIB)

%.o: %.f90
	$(FC) $(FC_FLAGS) -c $<

clean:
	rm -f tstdsygvic *.o *.txt
