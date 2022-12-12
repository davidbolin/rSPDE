toInclude = "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/INLA/include/"

obj = cgeneric_test.o cgeneric_test_aux.o

all : cgeneric_test.so

flags = -O2 -Wall -Wextra -fpic

%.o: %.c
	clang $(flags) -lstdc++ -Iinclude -I$(toInclude)  -c $^ -o $@

%.o: %.cpp
	clang++ $(flags)  -I$(toInclude) -I"/opt/homebrew/opt/eigen/include/eigen3" -c $^ -o $@

cgeneric_test.so: $(obj)
	clang++ -shared *.o -o cgeneric_test.so -lblas -llapack