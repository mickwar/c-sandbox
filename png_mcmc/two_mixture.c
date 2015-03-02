#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>             // getpid(), (process id)
//%#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

#define PI 3.14159265358979323846

const int NCOL = 4; // number of columns in read-in text file

// get random uniform draws on the interval [a,b]
double runif(double a, double b){
    return ((double)rand()/(double)RAND_MAX)*(b-a)+a;
    }

/*  n independent draws from a normal(mu, sig_i)
    array x is modified to contain random draws,
    int n is length of the array
    mu and sig are arrays specifying what the mean
    and standard deviations should be for each draw
    */
void rnorm(double* x, int n, double* mu, double* sig){
    double u1 = 0.0;
    double u2 = 0.0;
    int counter = 0;
    if (n % 2 == 0){
        do {
            u1 = runif(0.0, 1.0);
            u2 = runif(0.0, 1.0);
            // the counter is added AFTER the line is run
            x[counter++] = mu[counter] + sig[counter] *
                sqrt(-2*log(u1)) * cos(2*PI*u2);
            x[counter++] = mu[counter] + sig[counter] *
                sqrt(-2*log(u1)) * sin(2*PI*u2);
            } while (counter < n);
    } else {
        do {
            u1 = runif(0.0, 1.0);
            u2 = runif(0.0, 1.0);
            x[counter++] = mu[counter] + sig[counter] *
                sqrt(-2*log(u1)) * cos(2*PI*u2);
            x[counter++] = mu[counter] + sig[counter] *
                sqrt(-2*log(u1)) * sin(2*PI*u2);
            } while (counter < n-1);
        u1 = runif(0.0, 1.0);
        u2 = runif(0.0, 1.0);
        x[counter++] = mu[counter] + sig[counter] *
            sqrt(-2*log(u1)) * cos(2*PI*u2);
        }
    }

/*  the autotuning function, compute the acceptance rate (a) for
    each parameter in a window of iterations, then adjust parameter's
    candidate sigma according to: sig_new = sig_old * autotune(a) 
    */
double autotune(double x){
    double m = 0.25;    // target acceptance rate, 0 < m < 1
    double k = 2.50;    // maximum multiplicative change when x > m,
                        // when x < m, the smallest decrease is x=x/k
                        // k > 1
    if (x >= 0.25)
        return (1.0+(cosh(x-m)-1.0)*(k-1.0)/(cosh(m-1.0)-1.0));
    return (1.0-(cosh(x-m)-1.0)*(1.0-1.0/k)/(cosh(m)-1.0));
    }

double calc_mean(int* x, int iter, int window){
    int out = 0;
    for (int i=iter-window; i<iter; i++)
        out += x[i];
    return (1.0*out/window);
    }

// logit(pi) = alpha + beta * msu, solve for pi
double calc_pi(double x, double alpha, double beta){
    return (1.0/(1.0+exp(-(alpha+beta*x))));
    }

//  posterior
//  n is associated with length of msu and mep
//  m is the length of params
double calc_posterior(double* msu, float* mep, int n, double* params,
    int m, int ind_eta, int ind_alpha, int ind_beta, int* ind_alp_i,
    int* ind_bet_j, int* ind_theta, int* igroup, int* jgroup, int* i2j,
    double* hyperparam){
    double out = 0;
    // likelihood part
    for (int i=0; i<n; i++){
        if (mep[i] > 0){
            out += log(1-calc_pi(msu[i], params[ind_alpha], params[ind_beta])) +
                log(gsl_ran_beta_pdf(mep[i], params[ind_eta] * 
                params[ind_theta[igroup[i]-1]], params[ind_eta] *
                (1-params[ind_theta[igroup[i]-1]])));
        } else {
            out += log(calc_pi(msu[i], params[ind_alpha], params[ind_beta]));
            }
        }
    // priors
    // gsl_ran_gaussian_pdf(double x, double sigma), not mean so do x-mu
    out += log(gsl_ran_gaussian_pdf(params[ind_alpha], 1.0));   // alpha, N(0,1)
    out += log(gsl_ran_gaussian_pdf(params[ind_beta], 1.0));    // beta, N(0,1)
    out += (hyperparam[0]-1)*log(params[ind_eta]) -             // eta, Gamma
        params[ind_eta]/hyperparam[1];
    for (int i=0; i<igroup[n-1]; i++)                           // alp_i, Gamma
        out += (hyperparam[2]-1)*log(params[ind_alp_i[i]]) - 
            params[ind_alp_i[i]]/hyperparam[3];
    for (int j=0; j<jgroup[n-1]; j++)                           // bet_j, Gamma
        out += (hyperparam[4]-1)*log(params[ind_bet_j[j]]) - 
            params[ind_bet_j[j]]/hyperparam[5];
    for (int i=0; i<igroup[n-1]; i++)                           // theta, Beta
        out += log(gsl_ran_beta_pdf(params[ind_theta[i]], params[ind_alp_i[i]],
            params[ind_bet_j[i2j[i]]]));
    return (out);
    }

int main(int argc, char* argv[]){
    // read in file
    char* file = argv[1];
    int n = -1;
    FILE* myFile = fopen(file, "r");
    if (myFile == NULL)
        perror ("error opening file");
    while (!feof(myFile)){
        double dummy[NCOL];
        for (int j=0; j<NCOL; j++) // the number of columns must be known
            fscanf(myFile, "%lf", &dummy[j]);
        n++;
        }
    fclose(myFile);
    fopen(file, "r");
    double dat[n][NCOL];
    for (int i=0; i<n; i++)
        for (int j=0; j<NCOL; j++)
            fscanf(myFile, "%lf", &dat[i][j]);
    fclose(myFile);

    // variables to hold those from myFile
    int igroup[n]; // igroup[n-1] = n_i = length(i.group) from R code
    int jgroup[n]; // jgroup[n-1] = n_j = length(j.group) 
    double msu[n];
    float mep[n]; // use float since mep is only two decimals

    // with png data, we assume the individual number is sorted in
    // in column 1, then type number is column 2, msu is 3, mep is 4
    for (int i=0; i<n; i++){
        igroup[i] = (int)dat[i][0];
        jgroup[i] = (int)dat[i][1];
        msu[i] = (double)dat[i][2];
        mep[i] = (float)dat[i][3];
        }

    //
    // variable initialization for mcmc
    //
    // m = number of parameters
    // the n-1 are because c indexes from 0,...,n-1
    int ind_eta = 0;
    int ind_alpha = 1;
    int ind_beta = 2;
    int ind_alp_i[igroup[n-1]];
    for (int i=0; i<igroup[n-1]; i++)
        ind_alp_i[i] = ind_beta + i + 1;
    int ind_bet_j[jgroup[n-1]];
    for (int j=0; j<jgroup[n-1]; j++)
        ind_bet_j[j] = ind_alp_i[igroup[n-1]-1] + j + 1;
    int ind_theta[igroup[n-1]];
    for (int i=0; i<igroup[n-1]; i++)
        ind_theta[i] = ind_bet_j[jgroup[n-1]-1] + i + 1;
    int m = ind_theta[igroup[n-1]-1] + 1;
    int i2j[igroup[n-1]];
    for (int i=0; i<n; i++)
        i2j[igroup[i]-1] = jgroup[i]-1;

    int nburn = 100;
    int nmcmc = 150;
    int window = 50; // window to compute accept rate
    double params[m]; // matrix of parameter draws
    double cand[m]; // array to hold the candidate value
    double cand_vec[m]; // array to hold all cands for the iteration
    int accept[m][window]; // matrix of logicals for acceptances
    for (int i=0; i<nburn+nmcmc+1; i++)
        for (int j=0; j<m; j++)
            accept[j][i] = 0;
    double sigs[m]; // array for the candidate sigmas
    // make lower/upper different types like a float (if possible)
    double lower[m]; // array for lower bound on parameter i
    double upper[m]; // array for upper bound
    for (int i=0; i<m; i++){
        sigs[i] = 1.0; // sigs begin at 1
        lower[i] = 0.0;
        upper[i] = 1.0;
        // arbitrarily high
        if (i == 1 || i == 2)
            lower[i] = -10000000.0; 
        if (i < 3 + igroup[n-1] + jgroup[n-1])
            upper[i] = 10000000.0;
        }
    double hyperparam[6];
    hyperparam[0] = 5.0; // a_eta
    hyperparam[1] = 1.0; // b_eta
    hyperparam[2] = 5.0; // a_alpha
    hyperparam[3] = 1.0; // b_alpha
    hyperparam[4] = 5.0; // a_beta
    hyperparam[5] = 1.0; // b_beta
    // give starting values
    for (int i=0; i<m; i++){
        params[i] = 0.5;
        cand[i] = params[i];
        }
    double post = calc_posterior(msu, mep, n, params, m, ind_eta,
        ind_alpha, ind_beta, ind_alp_i, ind_bet_j, ind_theta, igroup,
        jgroup, i2j, hyperparam);
    double cand_post = calc_posterior(msu, mep, n, cand, m, ind_eta,
        ind_alpha, ind_beta, ind_alp_i, ind_bet_j, ind_theta, igroup,
        jgroup, i2j, hyperparam);

    srand(time(0)*getpid());
    for (int iter=1; iter<nburn+nmcmc+1; iter++){
//      printf("Burn-in: %d, NMCMC: %d, Iteration: %d\n", nburn, nmcmc,
//          iter);
        // the candidates for all parameters
        rnorm(cand_vec, m, params, sigs); 
        // j is the loop through each of the m parameters
        for (int j=0; j<m; j++){
            if (cand_vec[j] >= lower[j] && cand_vec[j] <= upper[j]){
                cand[j] = cand_vec[j];
                cand_post = calc_posterior(msu, mep, n, cand, m,
                    ind_eta, ind_alpha, ind_beta, ind_alp_i, ind_bet_j,
                    ind_theta, igroup, jgroup, i2j, hyperparam);
                if (log(runif(0.0, 1.0)) < cand_post - post){
                    post = cand_post;
                    params[j] = cand[j];
                    accept[j][iter] = 1;
                } else {
                    cand[j] = params[j];
                    }
                }
            if (floor(1.0*iter/window) == 1.0*iter/window &&
                iter <= nburn){
                sigs[j] *= autotune(calc_mean(&accept[j][0],
                    iter, window));
//              printf("%f\n", (calc_mean(&accept[j][0],
//                  iter, window)));
//              printf("%f\n", sigs[j]);

                }
            }
        }
//  for (int j=0; j<m; j++)
//      printf("%f ", sigs[j]);
//  printf("\n");
//  for (int i=1; i<nburn+nmcmc+1; i++){
//      for (int j=0; j<m; j++)
//          printf("%f ", accept[j][i]);
//      printf("\n");
//      }
//  for (int i=1; i<nburn+nmcmc+1; i++){
//      for (int j=0; j<m; j++)
//          printf("%f ", params[i][j]);
//      printf("\n");
//      }


    return (0);
    }
