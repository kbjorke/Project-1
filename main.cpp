#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cstdlib>

using namespace std;

int solve_tridiagonal_matrix_eq(int n);
int relative_error(int n);

int main(int argc, char* argv[])
{
    int n;
    int r;

    n = atoi(argv[1]);
    r = atoi(argv[2]);
    if(r = 1){
        solve_tridiagonal_matrix_eq(n);
   }
    if(r = 2){
        relative_error(n);
    }
    return 0;
}

int solve_tridiagonal_matrix_eq(int n)
{
    int i;
    char filename[30];
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
    v[n] = f[n-1]/d[n-1];
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

int relative_error(int n)
{
    int i;

    char filename[50];
    int n_filename;

    double v[n+2];
    double x[n+2];
    double u[n+2];
    double epsilon[n];
    double epsilon_max;
    double h;
    double h_log10;

    // Compute h and h_log10:
    h = 1.0/(n+1);
    h_log10 = log10(h);

    // Open and read x and v data from file:
    fstream infile;
    n_filename = sprintf(filename, "../Project1/data/solution_n%d.txt", n);
    infile.open(filename, ios::in);
    i = 0;
    while(!infile.eof())
    {
        infile >> x[i];
        infile >> v[i];

        i++;
    }
    infile.close();

    // Calcualte analytical solution u:
    for(i=0; i <= n+1; i++){
        u[i] = 1 - (1 - exp(-10))*x[i] - exp(-10*x[i]);
    }

    // Calculate relative error from i=1 to i=n and find max value for relative error:
    epsilon_max = -7;
    for(i=1; i <= n; i++){
        epsilon[i] = log10(fabs((v[i] - u[i])/u[i]));
        if(epsilon[i] > epsilon_max){
            epsilon_max = epsilon[i];
        }
    }

    // Write to file:
    ofstream outfile;
    outfile.open("relative_error.txt", ios::app);
    outfile << endl << n <<	"   " << h_log10 << "   " << epsilon_max;
    outfile.close();
    return 0;
}
