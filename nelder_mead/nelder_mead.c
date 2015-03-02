#include <stdio.h>          // printf
#include <stdlib.h>         // atoi, atof
#include <math.h>           // pi
#include <time.h>           // time (for seeding)
#include <unistd.h>         // getpid


const double ALPHA = 1.0;
const double GAMMA = 2.0;
const double RHO   = -0.5;
const double SIGMA = 0.5;

double himmelblau(double* x){
    double a = x[0];
    double b = x[1];
    return pow(pow(a, 2) + b - 11.0, 2) + pow(a + pow(b, 2) - 7.0, 2);
    }

double runif(double a, double b){
    return ((double)rand()/(double)RAND_MAX)*(b-a)+a;
    }

// get the order of the sorted index
void order(double* x, int n, int out[n]){
//  int out[n];
    int temp = 0;
    for (int i=0; i < n; i++)
        out[i] = i;

    for (int i=0; i < n-1; i++){
        for (int j=i+1; j < n; j++){
            if (x[out[j]] < x[out[i]]){
                temp = out[i];
                out[i] = out[j];
                out[j] = temp;
                }
            }
        }

//  for (int i=0; i < n; i++)
//      printf("%d ", out[i]);
//  printf("\n");

//  return out;
    }


int main(){
    srand(time(0)*getpid());

    // start with N+1 vertices, N = 2 for himmelblau function
    int N = 2;
    double x[N+1][N];   // x,y-coordindates
    double f_val[N+1];  // functional values at x, y
    for (int i=0; i < N+1; i++){
//      for (int j=0; j < N; j++){
//          x[i][j] = runif(-10.0, 10.0);
//          x[i][j] = runif(1.5, 3.5);
//          }
        x[i][0] = runif(-2.5, -3.5);
        x[i][1] = runif(2.5, 3.5);
        }

    // some initialize
    int sorted[N+1];
    int iter = 0;
    double x0[N]; // centroid
    double xr[N]; // reflection
    double xe[N]; // expanded
    double xc[N]; // contracted
    double fr;
    double fe;
    double fc;

    int yes_print = 0;

    while (iter < 100){
        iter++;
        if (yes_print)
            printf("Iteration: %d\n", iter);

        // (0) recompute functional values
        for (int i=0; i < N+1; i++)
            f_val[i] = himmelblau(x[i]);

        // (1) Order
        if (yes_print)
            printf("Order\n");
        order(f_val, N+1, sorted);
        // x[sorted[0]][.] is best point
        // x[sorted[N]][.] is worst point

        // (2) Calculate centroid (except for point N+1)
        if (yes_print)
            printf("Centroid\n");
        for (int i=0; i < N; i++){
            for (int j=0; j < N; j++){
                x0[j] = x0[j] + x[sorted[i]][j];
                }
            }
        for (int j=0; j < N; j++)
            x0[j] = 1.0 * x0[j] / N;

        // (3) Reflection
        if (yes_print)
            printf("Reflection\n");
        for (int j=0; j < N; j++)
            xr[j] = x0[j] + ALPHA*(x0[j] - x[sorted[N]][j]);
        fr = himmelblau(xr);
        if (f_val[sorted[0]] <= fr && fr < f_val[sorted[N-1]]){
            for (int j=0; j < N; j++)
                x[sorted[N]][j] = xr[j];
        } else {
            // (4) Expansion
            if (yes_print)
                printf("Expansion\n");
            if (fr < f_val[sorted[0]]){
                for (int j=0; j < N; j++)
                    xe[j] = x0[j] + GAMMA*(x0[j] - x[sorted[N]][j]);
                fe = himmelblau(xe);
                if (fe < fr){
                    for (int j=0; j < N; j++)
                        x[sorted[N]][j] = xe[j];
                } else {
                    for (int j=0; j < N; j++)
                        x[sorted[N]][j] = xr[j];
                    }
            } else {
                // (5) Contraction
                if (yes_print)
                    printf("Contraction\n");
                for (int j=0; j < N; j++)
                    xc[j] = x0[j] + RHO*(x0[j] - x[sorted[N]][j]);
                fc = himmelblau(xc);
                if (fc < f_val[sorted[N]]){
                    for (int j=0; j < N; j++)
                        x[sorted[N]][j] = xc[j];
                } else {
                    // (6) Reduction
                    if (yes_print)
                        printf("Reduction\n");
                    for (int i=1; i<N+1; i++){
                        for (int j=0; j<N; j++){
                            x[sorted[i]][j] = x[sorted[0]][j] + SIGMA * (x[sorted[i]][j] - 
                                x[sorted[0]][j]);
                            }
                        }
                    }
                }
            }
        }

//  for (int i=0; i < N+1; i++)
//      printf("%f ", f_val[i]);
//  printf("\n");


    for (int i=0; i<N+1; i++){
        for (int j=0; j<N; j++){
            printf("%f ", x[i][j]);
            }
        printf("\n");
        }
        

    

    return 0; 
    }
