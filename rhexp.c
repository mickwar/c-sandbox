#include <stdio.h>          // printf
#include <stdlib.h>         // atoi, atof
#include <math.h>           // pi
#include <time.h>           // time (for seeding)
#include <unistd.h>         // getpid

double runif(double a, double b){
    return ((double)rand()/(double)RAND_MAX)*(b-a)+a;
    }

double sum(double* x, int n){
    double out = 0;
    for (int i=0; i < n; i++)
        out += x[i];
    return (out);
    }

double rand_draw(double *lambda, double *p, int nparam){
    // select exponential from which to draw
    double rand = runif(0.0, 1.0);
    int this = nparam;
    for (int i=nparam-1; i>=0; i--){
        if ( p[i] > rand )
            this = i;
        }
    // get draw from this exponential
    rand = runif(0.0, 1.0);
    return -1.0/lambda[this]*log(rand);
    }

int main(int argc, char* argv[]){
    // check arguments
    // -- right number
    if ( argc % 2 != 0 ){
        printf("usage: rhyperexp n lambda_1 ... lambda_n p_1 ... p_n\n");
        return 1;
        }
    int n = atoi(argv[1]);
    int nparam = (argc-1) / 2;
    double lambda[nparam];
    double p[nparam];
    for (int i=0; i<nparam; i++){
        lambda[i] = atof(argv[i+2]);
        p[i] = atof(argv[nparam+i+2]);
        }
    // p's must sum to 1
    if (sum(p, nparam) != 1.0){
        printf("p must sum to 1.\n");
        return 1;
        }
    // make p array cumulative
    for (int i=1; i<nparam; i++)
        p[i] += p[i-1];

    // get the n draws
    for (int i=0; i<n; i++)
        printf("%f\n", rand_draw(lambda, p, nparam));
    return 0;
    }
