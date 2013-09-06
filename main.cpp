#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cstdlib>

using namespace std;

int solve_tridiagonal_matrix_eq(int n);

int main(int argc, char* argv[])
{
    int n;
    n = atoi(argv[1]);
    solve_tridiagonal_matrix_eq(n);
    return 0;
}

int solve_tridiagonal_matrix_eq(int n)
{
    int i;
    char filename[20];
    int n_filename;

    double u[n];
    double f[n];
    double v[n];
    double x[n];
    double h;

    // Set step length:
    h = 1.0/(n+1);

    // Initialize the arrays a, c and d:
    for(i=0; i <= n-1; i++){
        d[i] = 2;
    }

    // Construct the array x, which go from the value 0+h to 1-h with n steps:
    for(i=0; i <= n-1; i++){
        x[i] = h*(1+i);
    }

    // Compute the array f, which contains the source term for our eqation:
    for(i=0; i <= n-1; i++){
        f[i] = h*h*100*exp(-10*x[i]);
    }


    // Forward substitution:
    for(i=1; i <= n-1; i++){
        d[i] -= 1/d[i-1];
        f[i] += f[i-1]/d[i-1];
    }

    // Backward substitution:
    // First step:
    v[n-1] = f[n-1]/d[n-1];
    // Loop:
    for(i=n-2; i>=0; i--){
        v[i] = (f[i] + v[i+1])/d[i];
    }

    // Write to file:
    fstream myfile;
    n_filename = sprintf(filename, "solution_n%d.txt", n);
    myfile.open(filename, ios::out);
    for(i=0; i <= n-1; i++){
        myfile << x[i] << "   " << v[i] << endl;
    }
    myfile.close();
    return 0;
}
