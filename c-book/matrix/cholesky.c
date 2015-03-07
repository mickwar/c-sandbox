#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


double runif(double a, double b){
    return ((double)rand()/(double)RAND_MAX)*(b-a)+a;
    }

void mat_print(int n, int m, double x[n][m]){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            printf("%f ", x[i][j]);
            }
        printf("\n");
        }
    }

void transpose(int n, int m, double in[n][m], double out[m][n]){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            out[j][i] = in[i][j];
            }
        }
    }

void mat_mult(int a_n, int a_m, int b_n, int b_m,
    double a[a_n][a_m], double b[b_n][b_m], double out[a_n][b_m]){
    for (int i=0; i<a_n; i++){
        for (int j=0; j<b_m; j++){
            out[i][j] = 0;           
            for (int k=0; k<a_m; k++)
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }

// Computes L in X = LL' 
// For ill-conditioned matrices, may need to use LDL' decomposition,
// which eliminates the need to take square roots
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


int main(void){
    srand(time(0)*getpid());
    int n = 15;
    int p = 7;
    double x[n][p];
    double y[p][n];
    double z[p][p];
    double L[p][p];
    for (int i=0; i<n; i++){
        for (int j=0; j<p; j++){
            x[i][j] = runif(-2, 2);
            }
        }

    transpose(n, p, x, y);

    mat_mult(p, n, n, p, y, x, z);
    cholesky(p, z, L);

    mat_print(p, p, z);
    mat_print(p, p, L);

    

    }
