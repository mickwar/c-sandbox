#include <stdio.h>

// from the MCMC book
long int kiss (unsigned long *i, unsigned long *j,
    unsigned long *k){
    *j = *j ^ (*j << 17);
    *k = (*k ^ (*k << 18)) & 0X7FFFFFFF ;
    return  ((*i = 69069 * (*i) + 23606797) + 
        (*j = (*j >> 15)) + (*k = (*k >> 13))) ;
    }

int main(){
    unsigned long starti = 1;
    printf("%f\n", kiss(*starti,1,1));
    }
