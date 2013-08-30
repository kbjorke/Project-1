#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int solve_tridiagonal_matrix_eq(int n);

int main()
{ 
    solve_tridiagonal_matrix_eq(100);
    return 0;
}

int solve_tridiagonal_matrix_eq(int n)
{
    int i;

    double u[n];
    double f[n];
    double d[n];
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
    u[n-1] = f[n-1]/d[n-1];
    // Loop:
    for(i=n-2; i>=0; i--){
        u[i] = (f[i] + u[i+1])/d[i];
    }

    // Write to file:
    fstream myfile;
    myfile.open("solution.txt", ios::out);
    for(i=0; i <= n-1; i++){
        myfile << x[i] << "   " << u[i] << endl;
    }
    myfile.close();
    return 0;
}
