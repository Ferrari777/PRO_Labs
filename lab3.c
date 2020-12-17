#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <omp.h>

const long binaryNumber32 = 4294967296;
const long R = 1048576;
const long Rmod = 1048575;
const long Rsft = 20;
const long sleepConst = 800;
const long binaryNumber9 = 512;
const long mask = 1023;
const long binaryNumber = 4611686018427387904;

long BinPow(long a, long m, long n) {

	long res = 1;
	long i = 1;
	a = a % m;

	while (n > 0) {
		if (n & 1) {
			res = (res * a) % m;
		}
		a = (a * a) % m;
		n >>= 1;
	}

	return res;
}

long FindInverseNumber(long a, long m) {

	long ac, mc;
	long xi, yi;
	long xip, yip;
	long xc, yc;
	long q;
	long mod;
	long exch;

	exch = a;
	a = m;
	m = exch;

	xi = 1;
	yi = 0;
	xip = 0;
	yip = 1;
	ac = a;
	mc = m;
	while (1 == 1) {

		q = ac / mc;
		mod = ac % mc;

		xc = xi - q * xip;
		yc = yi - q * yip;

		xi = xip;
		yi = yip;
		xip = xc;
		yip = yc;

		ac = mc;
		mc = mod;

		if (mc == 0) {
			break;
		}

	}

	return yi;

}

long MultiplyNumbers(long a, long b, long m) {

	long aCurr;
	long bCurr;
	long tCurr;
	long resCurr;
	long mmin;
	long Rmin;

	aCurr = (a * R) % m;
	bCurr = (b * R) % m;
	tCurr = aCurr * bCurr;

	Rmin = -FindInverseNumber(R, m);
	mmin = FindInverseNumber(m, R);

	resCurr = (tCurr + m * ((tCurr * mmin) & Rmod)) >> Rsft;
	if (resCurr >= m) {
		resCurr = resCurr - m;
	}
	resCurr = (resCurr * Rmin) % m;

	return resCurr;
}

long BinPowMontgomeri(long a, long m, long n) {

	long res = 1;
	long i = 1;
	a = a % m;


	while (n > 0) {
		if (n & 1) {
			res = MultiplyNumbers(res, a, m);
		}
		a = MultiplyNumbers(a, a, m);
		n >>= 1;
	}


	if (res < 0) {
		res = m + res;
	}
	return res;
}


long LeftPow(long a, long m, long n) {

	long currBinaryNumber = binaryNumber;
	long res, ind;
	if ((n & currBinaryNumber) == 0) {
		res = 1;
	} else {
		res = a % m;
	}


	currBinaryNumber >>= 1;
	while (currBinaryNumber) {

		ind = n & currBinaryNumber;
		if (ind == 0) {
			res = (res * res) % m;
		} else {
			res = (res * res * a) % m;
		}
		currBinaryNumber >>= 1;

	}

	return res;
}

int main() {


	long a, m, n;
	long res;
    long a1, a2, nHigh, nLow;
    long res9, res32, resHigh, resLow;
	struct timeval start, end;
	double timeSnap;


	printf("Enter the value of a: ");
	scanf("%li", &a);
	printf("Enter the value of n: ");
	scanf("%li", &n);
	printf("Enter the value of m: ");
	scanf("%li", &m);

	//////////////////////////////////////////////////
    //       Realization of binary method           //
    //////////////////////////////////////////////////


	gettimeofday(&start, NULL);
	res = BinPow(a, m, n);
	gettimeofday(&end, NULL);
	timeSnap = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	
    printf("\nBinary method result: %li\n", res);
	printf("Binary method time: %f\n\n", timeSnap);

    nHigh = n >> 32;
	nLow = (n << 32) >> 32;
	a1 = a;
    a2 = a1;

	omp_set_dynamic(0);
	omp_set_num_threads(2);
	gettimeofday(&start, NULL);

	#pragma omp parallel sections
	{
    	#pragma omp section
    	{ 
		res32 = BinPow(a1, m, binaryNumber32);
		resHigh = BinPow(res32, m, nHigh);
		}
		#pragma omp section
		{
		resLow =  BinPow(a2, m, nLow);
		}
	}
	res = (resHigh * resLow) % m;
	gettimeofday(&end, NULL);
	timeSnap = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	printf("\nBinary method parallel result: %li\n", res);
	printf("Binary method parallel time: %.10f\n\n", timeSnap);

    
    //////////////////////////////////////////////////
    //       Realization of Montgomery method       //
    //////////////////////////////////////////////////

	gettimeofday(&start, NULL);
	res = BinPowMontgomeri(a, m, n);
	gettimeofday(&end, NULL);
	timeSnap = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	printf("\nMontgomery method result: %li\n", res);
	printf("Montgomery method time: %f\n\n", timeSnap);

    nHigh = n >> 10;
	nLow = n & mask;
	a1 = a;
    a2 = a1;


	omp_set_dynamic(0);
	omp_set_num_threads(2);
	gettimeofday(&start, NULL);

	#pragma omp parallel sections
	{
    	#pragma omp section
    	{ 
		res9 = BinPow(a1, m, binaryNumber9);
		resHigh = BinPowMontgomeri(res9, m, nHigh);
		}
		#pragma omp section
		{
		resLow =  BinPowMontgomeri(a2, m, nLow);
		}
	}
	res = (resHigh * resLow) % m;
	gettimeofday(&end, NULL);
	timeSnap = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	printf("\nMontgomery parallel method result: %li\n", res);
	printf("Montgomery parallel method time: %f\n\n", timeSnap);

    
    ////////////////////////////////////////////////////////
    //    Realization of "from right to left" method      //
    ////////////////////////////////////////////////////////

	gettimeofday(&start, NULL);
	res = LeftPow(a, m, n);
	gettimeofday(&end, NULL);
	timeSnap = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	printf("\n\"From right to left\" method result: %li\n", res);
	printf("\"From right to left\" method time: %f\n\n", timeSnap);

    nHigh = n >> 32;
	a1 = a;
    a2 = a1;

	omp_set_dynamic(0);
	omp_set_num_threads(2);
	gettimeofday(&start, NULL);

	#pragma omp parallel sections
	{
    	#pragma omp section
    	{ 
		res32 = LeftPow(a1, m, binaryNumber32);
		resHigh = LeftPow(res32, m, nHigh);
		}
		#pragma omp section
		{
		resLow =  LeftPow(a2, m, n);
		}
	}
	res = (resHigh * resLow) % m;

	gettimeofday(&end, NULL);
	timeSnap = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;


	printf("\n\"From right to left\" parallel method result: %li\n", res);
	printf("\"From right to left\" parallel method time: %f\n\n", timeSnap);

    return 0;
}