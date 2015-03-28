#include <math.h>

// Gaussian elimination
void solve(int n, double x[n][n], double out[n][n]){
    double temp = 0;
    // init out
    for (int k=0; k<n; k++){
        out[k][k] = 1;
        }
    // down the row (first time)
    for (int i=0; i<n; i++){
        temp = x[i][i];
        // divide row i by element i,i
        for (int j=0; j<n; j++){
                x[i][j] = x[i][j] / temp;
                out[i][j] = out[i][j] / temp;
            }
        // go through each all other i-1 rows and subtract
        for (int k=0; k<n; k++){
            if (k == i){
                if (k == n-1){
                    break;
                    }
                k++;
                }
            temp = x[k][i];
            for (int j=0; j<n; j++){
                x[k][j] = x[k][j] - x[i][j]*temp;
                out[k][j] = out[k][j] - out[i][j]*temp ;
                }
            }
        } 
    }
