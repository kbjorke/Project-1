/* Program to preform matrix multiplication of two random n x n
 * sized  matrices. Write to file the time used for the multiplication.
 *
 * Usage: $~ ./matrix_multiplication 1000 
 *
 * Important to add library -lrt when linking, to get the function
 * clock_gettime() to work.
 *
 * Made by: Kristian Bjørke
 * */
#include <iostream>
#include <armadillo>
#include "time.h"
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace arma;

double getUnixTime();

int main(int argc, char* argv[])
{
	/* Main program, takes the n value for the matrix sizes
	 * from commandline.
	 */
	int i;
	int j;
	int k;
	int n;

	char filename[50];

	double start, finish, duration;
	int hours, minutes;
	double seconds;

	n = atoi(argv[1]);

	mat B = randu<mat>(n,n);
	mat C = randu<mat>(n,n);
	mat A = zeros<mat>(n,n);
	
	// Loop for matrx multiplication, which is timed:
	start = getUnixTime();
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			for(k=0; k<n; k++){
				A(i,j) += B(i,k)*C(k,j);
			}
		}
	}
	finish = getUnixTime();

	duration = finish - start;

	// Converts duration [s] to hours, minuttes and second.
	hours = (int) duration/(60.0*60.0);
	duration = duration - hours*60*60;
	minutes = (int) duration/60.0;
	duration = duration - minutes*60;
	seconds = duration;

	// Write to file:
	sprintf(filename, "time_matrix_multiplication_n%d.txt", n);
	fstream file;
	file.open(filename, ios::out);
	file << "Time used on matrix multiplication for n = " << n << " :" << endl;
	file << "      " << hours << " h  " << minutes << " m  " << seconds << 
		" s" << endl;
	file.close();

	return 0;
}

double getUnixTime()
{
	/* Function to get current time, with precission of microsecond */
    struct timespec tv;

    if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;

    return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}
