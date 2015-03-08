#include <math.h>

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
double det(int n, double x[n][n]){
    double sum = 0.0;
    if (n == 1){
        return x[0][0];
    } else {
        for (int i=0; i<n; i++){
            double xsub[n-1][n-1];
            submatrix(0, i, n, x, xsub);
            sum += pow(-1.0, i) * x[0][i]*det(n-1, xsub);
            }
        return sum;
        }
    }
