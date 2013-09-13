#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <string>
#include "time.h"
#include <armadillo>


using namespace std;

void solve_tridiagonal_matrix_eq(int n, double *x, double *v);
int make_solution_file(int n, double *x, double *v);
int relative_error(int n, double *x, double *v);
int time_usage(int n, double *x, double *v);

int main(int argc, char* argv[])
{
    int n;
    int i;
    int solved = 0;

    n = atoi(argv[argc-1]);

    double *x = new double[n+2];
    double *v = new double[n+2];

    for(i=1; i < argc; i++){
        if((string(argv[i]).find("-") == 0 && string(argv[i]).find("t")!=string::npos) ||
                argc == 2){
            time_usage(n,x,v);
            solved = 1;
        }
        if((string(argv[i]).find("-") == 0 && string(argv[i]).find("s")!=string::npos) ||
                argc == 2){
            if(solved == 0){
                    solve_tridiagonal_matrix_eq(n, x, v);
            }
            make_solution_file(n, x, v);
            solved = 1;
        }
        if(string(argv[i]).find("-") == 0 && string(argv[i]).find("e")!=string::npos){
            if(solved == 0){
                solve_tridiagonal_matrix_eq(n,x,v);
            }
            relative_error(n, x, v);

        }
    }
    return 0;
}

void solve_tridiagonal_matrix_eq(int n, double *x, double *v)
{
    int i;

    double *d = new double[n];
    double *f = new double[n];
    double h;

    // Set step length:
    h = 1.0/(n+1);

    // Initialize the array d:
    for(i=0; i < n; i++){
        d[i] = 2;
    }

    // Construct the array x, which go from the value 0 to 1 with n+2 steps of distance h:
    for(i=0; i < n+2; i++){
        x[i] = h*i;
    }

    // Dirichlet boundary conditions:
    v[0] = 0;
    v[n+1] = 0;


    // Compute the array f, which contains the source term for our eqation:
    for(i=0; i < n; i++){
        f[i] = h*h*100*exp(-10*x[i+1]);
    }


    // Forward substitution:
    for(i=1; i < n; i++){
        d[i] -= 1/d[i-1];
        f[i] += f[i-1]/d[i-1];
    }

    // Backward substitution:
    // First step
    v[n] = f[n-1]/d[n-1];
    // Loop:
    for(i=n-1; i>0; i--){
        v[i] = (f[i-1] + v[i+1])/d[i-1];
    }
}

int make_solution_file(int n, double *x, double *v)
{
    int i;

    char filename[30];

    // Write to file:
    fstream myfile;
    sprintf(filename, "solution_n%d.txt", n);
    myfile.open(filename, ios::out);
    for(i=0; i <= n+1; i++){
        myfile << x[i] << "   " << v[i] << endl;
    }
    myfile.close();
    return 0;
}

int relative_error(int n, double *x, double *v)
{
    int i;

    double *u = new double[n+2];
    double *epsilon = new double[n];
    double epsilon_max;
    double h;
    double h_log10;

    // Compute h and h_log10:
    h = 1.0/(n+1);
    h_log10 = log10(h);

    // Calculate analytical solution u:
    for(i=0; i<n+2; i++){
        u[i] = 1 - (1 - exp(-10))*x[i] - exp(-10*x[i]);
    }

    // Calculate relative error from i=1 to i=n and find max value for relative error:
    epsilon_max = -10;
    for(i=1; i<n+1; i++){
 //       cout << i << " " << v[i] << " " << u[i] << " " << v[i]-u[i] << " " << (v[i] - u[i])/u[i] << endl;
        epsilon[i-1] = log10(fabs((v[i]-u[i])/u[i]));
        if(epsilon[i-1] > epsilon_max){
            epsilon_max = epsilon[i-1];
        }
    }

    // Write to file:
    ofstream outfile;
    outfile.open("relative_error.txt", ios::app);
    outfile << endl << n <<	"   " << h_log10 << "   " << epsilon_max;
    outfile.close();
    return 0;
}

double getUnixTime()
{
    struct timespec tv;

    if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;

    return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}

void time_unit(double duration, double *time_value, char *time_unit)
{
    double duration_log = log10(duration);
    double value;

    if(duration_log>0){
        time_value = &duration;
        time_unit = "seconds";
    }
    else if(duration_log>-3){
        value = duration*1e3;
        time_value = &value;
        time_unit = "milliseconds";
    }
    else if(duration_log>-6){
        value = duration*1e6;
        time_value = &value;
        time_unit = "microseconds";
    }
}


int time_usage(int n, double *x, double *v)
{
    using namespace arma;

    double start_tridiag, finish_tridiag, start_LU, finish_LU;
    double duration_tridiag, duration_LU;
    double *time_tridiag, *time_LU;
    char *unit_tridiag, *unit_LU;

    start_LU = getUnixTime();

    int i, j;
    double h;

    mat A = zeros<mat>(n,n);
    vec f(n);
    vec x_LU(n);
    mat L;
    mat U;
    vec v_LU;
    vec z;

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i == j){
                A(i,j) = 2;
            }
            else if(fabs(i-j) == 1){
                A(i,j) = -1;
            }
            else{
                A(i,j) = 0;
            }
        }
    }

    h = 1.0/(n+1);

    for(i=0; i < n; i++){
        x_LU[i] = h*(i+1);
    }

    for(i=0; i < n; i++){
        f[i] = h*h*100*exp(-10*x_LU[i]);
    }

    lu(L,U,A);
    solve(z, L,f);
    solve(v_LU, U,z);

    finish_LU = getUnixTime();

    duration_LU = finish_LU - start_LU;

    start_tridiag = getUnixTime();
    solve_tridiagonal_matrix_eq(n, x, v);
    finish_tridiag = getUnixTime();

    duration_tridiag = finish_tridiag - start_tridiag;

    //time_unit(duration_tridiag, time_tridiag, unit_tridiag);

    cout << "n = " << n << " :" << endl;
    cout << "Time for tridiagonl matrix method: " << duration_tridiag << " s" << endl;
    cout << "Time for LU decomposition method : " << duration_LU << " s" << endl;

    return 0;
}


