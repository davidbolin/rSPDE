toInclude = ${R_LIBRARY_DIR}/INLA/include/

sources = $(wildcard *.c)

all : $(sources:.c=.o) so

%.o : %.c
	@echo gnerating $@
	gcc -c -O3 -Wall -Wextra -fpic -lblas -I$(toInclude) -I/opt/homebrew/Cellar/gsl/2.7.1/include $^ -o $@

so : $(sources:.c=.o)
	@echo generate .so
	gcc -shared -lblas -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl *.o -o rspde_cgeneric_models.so


clean :
	@##-rm -f *.o

.PHONY: all clean