#include <stdio.h>

// print a matrix
void mat_print(int n, int m, double x[n][m]){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            printf("%f ", x[i][j]);
            }
        printf("\n");
        }
    }

// transpose the n by m matrix x so it is m by n
void transpose(int n, int m, double in[n][m], double out[m][n]){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            out[j][i] = in[i][j];
            }
        }
    }

// matrix multiplication
// a is an a_n by a_m matrix and b is a b_n by b_m matrix.
// hence we must have a_m == b_n (could probably then remove
// one of those variables
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

int main () {
    int n = 4;
    int m = 7;
    double matA[n][m];
    double matB[m][n]; // for the transpose
    double matC[n][n]; // for the multipication

    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            matA[i][j] = (i+1.0)*(j+1.0);
            matB[j][i] = 0;
            }
        }

    transpose(n, m, matA, matB);
    mat_mult(n, m, m, n, matA, matB, matC);

    mat_print(n, m, matA);
    printf("\n");
    mat_print(m, n, matB);
    printf("\n");
    mat_print(n, n, matC);

    return 0;
    }
