toInclude = "/Library/Frameworks/R.framework/Versions/4.2/Resources/library/INLA/include/"
  
obj = cgeneric_test_aux.o cgeneric_test.o 

all : cgeneric_test.so
  
flags = -O2 -Wall -Wextra -fpic

cgeneric_test.o: cgeneric_test.c 
	gcc $(flags)  -Iinclude -I$(toInclude)  -c $^ -o $@ -lstdc++ 

cgeneric_test_aux.o: cgeneric_test_aux.cpp
	g++ $(flags)  -I$(toInclude)  -c $^ -o $@ 

cgeneric_test.so: $(obj)
	g++ -shared $(flags)  -o cgeneric_test.so $(obj) -lblas -llapack -lstdc++ 