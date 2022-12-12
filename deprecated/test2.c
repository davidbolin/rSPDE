#include "stdio.h"
#include <math.h>
#include <stdlib.h>
#include "assert.h"
#include "test.h"

#ifdef __cplusplus
#error C++ compiler, but should be C
#endif

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define DCalloc(n_, type_) (type_ **)calloc((n_), sizeof(type_ *))
#include <stdio.h>
#define SQR(x) ((x)*(x))

void test_function(double *entries, int *i, int *j, int size);
 
int main () {

    double *entries;
    int *i, *j;

    int k;

    entries = Calloc(5, double);
    i = Calloc(5, int);
    j = Calloc(5, int);

    i[0] = 4;
    i[1] = 3;
    i[2] = 2;
    i[3] = 1;
    i[4] = 0;

    j[0] = 4;
    j[1] = 3;
    j[2] = 2;
    j[3] = 1;
    j[4] = 0;

    entries[0] = 33;
    entries[1] = 27;
    entries[2] = 9;
    entries[3] = 3;
    entries[4] = 1;

    test_function(entries, i, j, 7);

    for(k = 0; k < 5; k++){
        printf("Entries[%d] = %f\n", k+1, entries[k]);
    }

   return 0;
}

// clang++ -I/opt/homebrew/opt/eigen/include/eigen3 -stdlib=libc++ -O -c test.cpp -o test.o
// clang -lstdc++ test.o test2.c -o test2