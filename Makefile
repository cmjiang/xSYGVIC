include make.inc

all: example_

ifdef LAPACK_DIR
lib = libgsep
else
lib = libgsep liblapack
endif

example_: $(lib)
	(cd example; $(MAKE))

liblapack:
	(cd lapacklib; $(MAKE))

libgsep:
	(cd src; $(MAKE))
clean:
	rm -f *.txt
	(cd src; $(MAKE) clean)
	(cd example; $(MAKE) clean)
	(cd lapacklib; $(MAKE) clean)
