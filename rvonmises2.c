#include <stdio.h>          // printf
#include <stdlib.h>         // atoi, atof
#include <math.h>           // pi
#include <time.h>           // time (for seeding)
#include <unistd.h>         // getpid

#define PI 3.14159265358979323846

int argcheck(double dvonarg, char* type){
    if (type == "size"){
        if (dvonarg <= 0)
            return 1;
        }
    if (type == "pi"){
        if (dvonarg <= -1*PI || dvonarg > PI)
            return 1;
        }
    if (type == "inf"){
        if (dvonarg <= 0)
            return 1;
        }
    if (type == "zero"){
        if (dvonarg == 0)
            return 1;
        }
    return 0;
    }

double dvonmises(double phi, double psi, double mu, double nu,
    double kappa1, double kappa2, double lambda){
    // phi/psi bound check
    int argument_check[2] = {0, 0};
    argument_check[0] = argcheck(phi, "pi");
    argument_check[1] = argcheck(psi, "pi");
    for (int i=0; i<2; i++){
        if (argument_check[i]){
            return -RAND_MAX;
            }
        }
    // log-density evaluation (not normalized)
    return kappa1*cos(phi-mu) + kappa2*cos(psi-nu)+
        lambda*sin(phi-mu)*sin(psi-nu);
    }

// Get the sign of a number
double sign(double x){
    if (x > 0)
        return +1.0;
    return -1.0;
    }

double runif(double a, double b){
    return ((double)rand()/(double)RAND_MAX)*(b-a)+a;
    }

// Return maximum of two values (for bimodal case)
double max(double x, double y){
    if (x >= y)
        return x;
    return y;
    }

// get the mode
double get_mode(double mu, double nu, double kappa1, double kappa2, double lambda){
    double mode_val = 0;
    if (kappa1*kappa2 > lambda*lambda){
        // unimodal
        mode_val = dvonmises(mu,nu,mu,nu,kappa1,kappa2,lambda);
    } else {
        // bimodal fabs: floating point absolute value
        double phi0 = sign(lambda)*acos(kappa1/fabs(lambda)*sqrt((lambda*
            lambda+kappa2*kappa2)/(lambda*lambda+kappa1*kappa1)));
        double psi0 = acos(kappa2/fabs(lambda)*sqrt((lambda*
            lambda+kappa1*kappa1)/(lambda*lambda+kappa2*kappa2)));
        double locs[2][2] = {{mu+phi0, nu+psi0}, {mu-phi0, nu-psi0}};
        // make sure the modes to check are in (-pi, pi]
        for (int i=0; i<2; i++){
            for (int j=0; j<2; j++){
                if (locs[i][j] > PI)
                    locs[i][j] -= 2*PI;
                if (locs[i][j] < -PI)
                    locs[i][j] += 2*PI;
                }
            } 
        mode_val = max(
            dvonmises(locs[0][0],locs[0][1],mu,nu,kappa1,kappa2,lambda),
            dvonmises(locs[1][0],locs[1][1],mu,nu,kappa1,kappa2,lambda));
        } 
    return mode_val;
    }

// doubles over as an R function
void rejection_sampler(int* n, double* mu, double* nu, double* kappa1, double* kappa2,
    double* lambda, double* output){
    srand(time(0)*getpid());
    // get height at mode (since uniform envelope, this is also like alpha)
    double mode_val = get_mode(*mu, *nu, *kappa1, *kappa2, *lambda);
    int counter = 0;
    double draw1 = 0.0;
    double draw2 = 0.0;
    double fun_value = 0.0;
    double prob_check = 0.0;
    do {
        draw1 = runif(-PI, PI);
        draw2 = runif(-PI, PI);
        fun_value = dvonmises(draw1, draw2, *mu, *nu, *kappa1, *kappa2, *lambda);
        prob_check = log(runif(0, 1)); 
        // accept with probability [0, f(x)/e(x)], but in log
        if (prob_check <= (fun_value - mode_val)){
            output[counter++] = draw1;
            output[counter++] = draw2;
            }
        } while (counter < 2 * *n);
    }

int main(int argc, char* argv[]){
    // Check for right number of arguments 
    if (argc != 7){
        printf("usage: rvonmises n mu nu kappa1 kappa2 lambda\n");
        return 1;
        }
    // Assign command line arguments to variabes
    int n = atoi(argv[1]);
    double mu = atof(argv[2]);
    double nu = atof(argv[3]);
    double kappa1 = atof(argv[4]);
    double kappa2 = atof(argv[5]);
    double lambda = atof(argv[6]);
    // Check variables in right bounds
    int argument_check[6] = {0, 0, 0, 0, 0, 0};
    argument_check[0] = argcheck(n, "size");
    argument_check[1] = argcheck(mu, "pi");
    argument_check[2] = argcheck(nu, "pi");
    argument_check[3] = argcheck(kappa1, "inf");
    argument_check[4] = argcheck(kappa2, "inf");
    argument_check[5] = argcheck(lambda, "zero");
    for (int i=0; i<6; i++){
        if (argument_check[i]){
            printf("invalid bounds on a parameter.\n");
            return 1;
            }
        }
    // Set seed for random number.
    // srand(time(0));

    double* output = malloc(2.0*n*sizeof(double));
    rejection_sampler(&n, &mu, &nu, &kappa1, &kappa2, &lambda, output);
    // print to bash
    for (int i=0; i<2*n; i+=2){
        printf("%f %f\n", output[i], output[i+1]);
        }

    free (output);
    return 0;
    }

