#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


void cholesky(int n, double x[n][n], double out[n][n]){
    for (int j=0; j<n; j++){

        // diagonal component
        out[j][j] = x[j][j];
        for (int k=0; k<j; k++){
            out[j][j] -= pow(out[j][k], 2);
            }
        out[j][j] = sqrt(out[j][j]);
        
        // off diagonals
        for (int i=j+1; i<n; i++){
            out[i][j] = 0.0;
            for (int k=0; k<j; k++){
                out[i][j] -= out[i][k] * out[j][k];
                }
            out[i][j] = 1.0 * (out[i][j] + x[i][j]) / out[j][j];
            }

        }
    }
