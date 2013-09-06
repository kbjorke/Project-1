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

    double d[n];
    double f[n];
    double v[n+2];
    double x[n+2];
    double h;

    // Set step length:
    h = 1.0/(n+1);

    // Initialize the arrays a, c and d:
    for(i=0; i <= n-1; i++){
        d[i] = 2;
    }

    // Construct the array x, which go from the value 0 to 1 with n+2 steps of distance h:
    for(i=0; i <= n+1; i++){
        x[i] = h*i;
    }

    // Dirichlet boundary conditions:
    v[0] = 0;
    v[n+1] = 0;


    // Compute the array f, which contains the source term for our eqation:
    for(i=0; i <= n-1; i++){
        f[i] = h*h*100*exp(-10*x[i+1]);
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
        v[i+1] = (f[i] + v[i+2])/d[i];
    }

    // Write to file:
    fstream myfile;
    n_filename = sprintf(filename, "solution_n%d.txt", n);
    myfile.open(filename, ios::out);
    for(i=0; i <= n+1; i++){
        myfile << x[i] << "   " << v[i] << endl;
    }
    myfile.close();
    return 0;
}
