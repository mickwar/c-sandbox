#include <stdio.h>

/* Suppose you need N booleans */

int main(){
    /* long long int is 64 bits */
    unsigned long long int test = 9223372036854775807;
    unsigned long long int bit = -1;
    bit = bit / 2 + 1;
    for (int i = 0; i < sizeof(test) * 8; i++){
        if (test >= bit){
            test = test - bit;
            printf("1 ");
        } else {
            printf("0 ");
            }
        
        printf("%llu\n", bit);
        bit = bit / 2;
        }
    printf("\n");
    return 0;
    }
