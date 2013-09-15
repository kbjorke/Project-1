/* Program for Project 1, FYS3150.
 *
 * Made by: Kristian Bjørke
 * */
#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <string>
#include "time.h"
#include <armadillo>


using namespace std;
using namespace arma;


// Functions declarations, actual functions follow int main() in this order:
void solve_tridiagonal_matrix_eq(int n, double *x, double *v);
int make_solution_file(int n, double *x, double *v);
int relative_error(int n, double *x, double *v);
vec solve_by_LU_decomp(int n, double *x);
int make_solution_file(int n, double *x, vec v_LU);
int relative_error(int n, double *x, vec v_LU);
int time_usage(int n, double *x, double *v, vec v_LU);
double getUnixTime();

int main(int argc, char* argv[])
{
    /* Main program.
     * Takes input from command line:
     *
     * n : The size of the n x n matrix to describe the problem is
     *     given as the last commandline arguments
     *
     * Flags:
     * 			-s  --  Solves the problem and writes to file
     *				 	"solution_n<n_value>.txt
     * 					( -sLU solves with LU decomposition)
     *          -e  -- 	Calculates relative error (epsilon) and appends
     * 					the values of n, log_10(h) and log_10(epsilon)
     *					to the file "relative_error.txt"
     *  					( -eLU calculates with LU decompostion)
     * 			-t  --  Runs and times both the LU decomposition and the
     * 					optimized method and writes the results to the
     * 					screen.
     *
     * Use: $~ ./Project1 <flag[s]> <n-value>
     * */
    int n;
    int i;
    // Variable to avoid calling solve function two times in flag handler:
    int solved = 0;

    n = atoi(argv[argc-1]);

    // Initialize x and v array for optimized method:
    double *x = new double[n+2];
    double *v = new double[n+2];

    // Initialize armadillo vector v_LU for the LU method:
    vec v_LU = zeros<vec>(n+2,1);


    // Handles flags and calls desired functions:
    for(i=1; i < argc; i++){
        if((string(argv[i]).find("-") == 0 && string(argv[i]).find("t")!=string::npos)){
            // Time flag (-t):
            time_usage(n,x,v, v_LU);
            solved = 1;
        }
        if((string(argv[i]).find("-") == 0 && string(argv[i]).find("s")!=string::npos) ||
                argc == 2){
            // Solve flag (-s) or no flag:
            if(string(argv[i]).find("LU")!=string::npos){
                // Solve by LU decomposition (-sLU)
                v_LU.subvec(1,n) = solve_by_LU_decomp(n, x);
                make_solution_file(n, x, v_LU);
            }
            else{
                // Solve by optimized method:
                if(solved == 0){
                    solve_tridiagonal_matrix_eq(n, x, v);
                }
            make_solution_file(n, x, v);
            solved = 1;
            }
        }
        if(string(argv[i]).find("-") == 0 && string(argv[i]).find("e")!=string::npos){
            // Error-calulation flag (-e):
            if(string(argv[i]).find("LU")!=string::npos){
                // Error-calculation for LU decomposition solution (-eLU):
                v_LU.subvec(1,n) = solve_by_LU_decomp(n, x);
                relative_error(n, x, v_LU);
            }
            else{
            // Error-calculation by optimized method:
            if(solved == 0){
                solve_tridiagonal_matrix_eq(n,x,v);
            }
            relative_error(n, x, v);
            }
        }
    }
    return 0;
}

void solve_tridiagonal_matrix_eq(int n, double *x, double *v)
{
    /* Function for solving the tridiagonal matrix equation by
     * the optimal method we introduced in the rapport.
     * Takes the arrays x and v by reference and the size variable
     * n as input.
     *
     * Important difference from algorithm in rapport:
     * In this function the variable f-tilde is called f.
     * */
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
    /* Function to write solution data to file.
     * Used for solution by the optimized method and
     * takes the size varible n and the arrays x and v
     * by reference as input.
     *
     * Prints to file: "solution_n<n_value>.txt"
     * */
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
    /* Function to calculate relative error for solutions made
     * by the optimized method. Taks size varible n and the
     * arrays x and v by reference as input.
     *
     * Appends n-value, log_10(h) and log_10(epsilon) to file
     * 			"relative_error.txt"
     * */
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
        epsilon[i-1] = log10(abs((v[i]-u[i])/u[i]));
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



vec solve_by_LU_decomp(int n, double *x)
{
    /* Function for solving the tridiagonal matrix equation by
     * LU decomposition and armadillo types and functions.
     * Takes the array x reference and the size variable
     * n as input.
     * Returns solution as armadillo vector, does not include
     * boundary points v_0 = v_n+1 = 0.
     *
     * Important difference from algorithm in rapport:
     * In this function the variable f-tilde is called f.
     * */
    using namespace arma;

    int i, j;
    double h;

    mat A = zeros<mat>(n,n);
    vec f(n);
    mat L;
    mat U;
    vec z;

    // Initialize tridiagonal matrix A:
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i == j){
                A(i,j) = 2;
            }
            else if(fabs(i-j) == 1){
                A(i,j) = -1;
            }
        }
    }

    // Step length:
    h = 1.0/(n+1);

    // Initialze vector x:
    for(i=0; i < n+2; i++){
        x[i] = h*i;
    }

    // Calculate source term:
    for(i=0; i < n; i++){
        f[i] = h*h*100*exp(-10*x[i+1]);
    }

    // Solves by preforming LU decomposition of A, then solving
    // the two matrix equations L*z = f , U*v = z:
    lu(L,U,A);
    z = solve(L,f);
    return solve(U,z);
}


int make_solution_file(int n, double *x, vec v_LU)
{
    /* Function to write solution data to file.
     * Used for solution by the LU decompositon method and
     * takes the size varible n and, array x by reference
     * and the armadillo vector v_LU as input.
     *
     * Prints to file: "solution_n<n_value>.txt"
     * */
    int i;

    char filename[30];
    // Write to file:
    fstream myfile;
    sprintf(filename, "solution_n%d.txt", n);
    myfile.open(filename, ios::out);
    for(i=0; i <= n+1; i++){
        myfile << x[i] << "   " << v_LU[i] << endl;
    }
    myfile.close();
    return 0;
}


int relative_error(int n, double *x, vec v_LU)
{
    /* Function to calculate relative error for solutions made
     * by the LU decompostion method. Taks size varible n, the
     * array x by reference and the armadillo vector v_LU as input.
     *
     * Appends n-value, log_10(h) and log_10(epsilon) to file
     * 			"relative_error.txt"
     * */
    int i;

    vec u(n);
    double h;
    double h_log10;

    // Compute h and h_log10:
    h = 1.0/(n+1);
    h_log10 = log10(h);

    // Calculate analytical solution u:
    for(i=0; i<n; i++){
        u[i] = 1 - (1 - exp(-10))*x[i+1] - exp(-10*x[i+1]);
    }

    // Write to file:
    ofstream outfile;
    outfile.open("relative_error.txt", ios::app);
    outfile << endl << n <<	"   " << h_log10 << "   "
            << max(log10(abs((v_LU.subvec(1, n)-u)/u)));
    outfile.close();
    return 0;
}





int time_usage(int n, double *x, double *v, vec v_LU)
{
    /* Function to test time usage of both the optimized method and the
     * LU decomposition method for the size varible n.
     * Takes n, the arrays x and v by refferece and a armadillo vector
     * v_LU. The array v and the vector v_LU should be solved solutions.
     *
     * Prints the result of the timing to the screen.
     * */
    double start_tridiag, finish_tridiag, start_LU, finish_LU;
    double duration_tridiag, duration_LU;

    // Test time of LU decompostions method:
    start_LU = getUnixTime();
    v_LU = solve_by_LU_decomp(n,x);
    finish_LU = getUnixTime();

    duration_LU = finish_LU - start_LU;

    // Test time of optimized method:
    start_tridiag = getUnixTime();
    solve_tridiagonal_matrix_eq(n, x, v);
    finish_tridiag = getUnixTime();

    duration_tridiag = finish_tridiag - start_tridiag;

    // Prints results to screen:
    cout << "n = " << n << " :" << endl;
    cout << "Time for tridiagonal matrix method: " << duration_tridiag << " s" << endl;
    cout << "Time for LU decomposition method  : " << duration_LU << " s" << endl;

    return 0;
}


double getUnixTime()
{
    /* Function to get accurate time, prescission to about some microseconds.
     *
     * Important: to get the function clock_gettime() to work link to library
     * -lrt when linking.
     * */
    struct timespec tv;

    // Get time:
    if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;

    // Returns time in seconds:
    return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}
