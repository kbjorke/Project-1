#include <iostream>
#include <armadillo>
#include "time.h"
#include <fstream>

using namespace std;
using namespace arma;

double getUnixTime();

int main()
{
	int i;
	int j;
	int k;
	int n;

	char filename[50];

	double start, finish, duration;
	int hours, minutes;
	double seconds;

	n = 1e3;

	mat B = randu<mat>(n,n);
	mat C = randu<mat>(n,n);
	mat A = zeros<mat>(n,n);
	
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
    struct timespec tv;

    if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;

    return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}
