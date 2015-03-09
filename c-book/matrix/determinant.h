#include <math.h>
#include <string.h>

// Get a submatrix of x by removing row 'remove_row' and
// removing column 'remove_col'
void submatrix(int remove_row, int remove_col, int n, double x[n][n],
    double out[n-1][n-1]){
    int p = 0;
    int q = 0;
    for (int i=0; i<n-1; i++){
        p = i;
        if (i >= remove_row)
            p++;
        for (int j=0; j<n-1; j++){
            q = j;
            if (j >= remove_col)
                q++;
            out[i][j] = x[p][q];
            }
        }
    }

// Recursive algorithm for calculating determinant
double determinant(int n, double x[n][n], char* matrix_type){
    double sum = 0.0;

    // for diagonal or triangular matrices, take the product of the diagonal elements
    if (strcmp(matrix_type, "triangle") == 0){
        sum = 1.0;
        for (int i=0; i<n; i++){
            sum *= x[i][i];
            }
        return sum;
        }

    // otherwise compute recursively
    if (n == 1){
        return x[0][0];
    } else {
        for (int i=0; i<n; i++){
            double xsub[n-1][n-1];
            submatrix(0, i, n, x, xsub);
            sum += pow(-1.0, i) * x[0][i]*determinant(n-1, xsub, matrix_type);
            }
        return sum;
        }
    }
