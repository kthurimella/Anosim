/*******************************************************/
/* To compile                                          */
/* % g++ similar.cpp -o anosim                         */
/*                                                     */
/* To run                                              */
/* anosim -d <sample_filename>                         */
/*        -g <group_label_filename>                    */
/*       [-p no-of-permutations]                       */
/*                                                     */
/* example 1:                                          */
/* % ./anosim -d test.txt -g gtest.txt -p 3            */
/*                                                     */
/* example 2:                                          */
/* % ./anosim -g gd.txt -d d1.txt -p 2000000           */
/*                                                     */
/* example 3:                                          */
/* % ./anosim -g gd.txt -d d2.txt                      */
/*                                                     */
/*******************************************************/

#include <fstream>
#include <iostream> 
#include <sstream> 
#include <string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

using namespace std; 

//structure to hold distance, row, column values
//this is needed so that when sorted on dist values
//row, column values move with them

struct dgPair {
	double dist;
	int rowNo;
	int columnNo;
};

//standard quicksort routine with a random pivot
void quickSort(dgPair arr[], int left, int right);

//function to compute R value
double anosim_stat( dgPair dataArray[], int data_size, 
		int gArray[], int gArray_size);

int main(int argc, char *argv[]){
	int n;           //number of samples
	int NO_OF_PERMS = 1000000;

	if ( argc < 5 ) // argc should be at least 5 (0 through 4)
		// for correct execution
	{               // argv[0] is the program name
		cout << "usage: \n" <<
				"anosim -d <sample_filename> \n"
				"       -g <group_label_filename>  \n"
				"       [-p no-of-permutations]  \n";
		cout << "example: \n"
				<< "./anosim -d test.txt -g gtest.txt -p 3\n";
		exit(1);
	}

	//read command line arguments
	ifstream f, g;
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			case 'd':
				f.open (argv[++i]);
				break;
			case 'g':
				g.open (argv[++i]);
				break;
			case 'p':
				NO_OF_PERMS = atoi(argv[++i]);
			default: /* do nothing */ break;
			}
		}
	}

	//get the top line--column descriptor and compute n from it
	if (f.is_open())
	{
		string descriptor;
		getline(f, descriptor); //read first line into descriptor variable
		stringstream descStream(descriptor); //convert it to a stream
		while(descStream) { //count how many columns
			// assumption that last number is n
			descStream >> n;
		}
	}
	else
	{
		cout<<"Could not open the sample file\n";
		exit(1);
	}

	//populate the data array from the file
	double d, data[n][n];
	for (int i=0;i<n;i++)
		for (int j=0;j<(n+1);j++)
			if (!(f>>d)) //read d from the file stream f
			{
				cout << "insufficient data \n";
				exit(1);
			}
			else if (j != 0) data[i][j-1] = d; //ignore the first column

	//populate the group label array from the file
	int groupLabels[n];
	if (g.is_open())
	{
		int x;
		for (int i=0;i<n;i++)
			if ( g>>x ) groupLabels[i] = x; //read x from the file stream g
			else {
				cout << "insufficient data \n";
				exit (1);
			}
	}
	else {
		cout<< "Could not open the group label file\n";
		exit(1);
	}

	//--------populate the upper-triangular matrix)-------------
	int utSize = n*(n-1)/2; // total = n*n, diagonal = n, upper half = (total - diagonal)/2
	dgPair dg[utSize]; //upper triangular matrix in 1 dimensional form
	int k=0;

	for (int i=0;i < (n-1);i++)
		for (int j=(i+1);j < n;j++)
		{
			dg[k].dist = data[i][j];
			dg[k].rowNo = i;
			dg[k++].columnNo = j;
		}

	//sort the dg array by dist values (first data member)
	quickSort(dg,0,utSize-1);

	//fix ties
	double avg, sum;
	int l, r;

	for (int i=0; i < utSize-1;i++)
	{
		int runLength;
		avg = sum = 0;
		l = r = i;

		while ( i < (utSize -1) && dg[i].dist == dg[i+1].dist){
			sum += i+1; //raw rank values
			r++;
			i++;
		}
		sum += i+1;

		runLength = r-l+1;
		avg = sum/runLength;
		for (int i=0; i<runLength;i++)
			dg[r-i].dist = avg;
	}
	//address the last value
	if (l == r) dg[utSize-1].dist = utSize;

	// compute R while permuting group label array
	double R, R_true, P;
	int noOfPermutes = 1,
			noOfBiggerR = 0;

	R_true = anosim_stat(dg, utSize, groupLabels, n);

	cout << "R_true = " << R_true << endl;

	// for next_permutation documentation, refer to
	// http://www.cplusplus.com/reference/algorithm/next_permutation/

	while ( noOfPermutes < NO_OF_PERMS &&
			next_permutation (groupLabels,groupLabels+n) ){
		R = anosim_stat(dg, utSize, groupLabels, n);
		if (R > R_true)  noOfBiggerR++;
		noOfPermutes++;
	}

	//compute P vale
	P = (double) noOfBiggerR /  noOfPermutes;

	cout << "no of permutations = " << noOfPermutes << endl;
	cout << "p-value = " << P << endl;
}

double anosim_stat( 
		dgPair dataArray[],
		int data_size,
		int gArray[],
		int gArray_size) {
	double R, R_w = 0, R_b = 0, R_wSum = 0, R_bSum = 0;
	int wCount = 0, bCount = 0;

	//scan dataArray and compute R_w and R_b
	for (int i=0;i< data_size;i++)
		if (gArray[dataArray[i].rowNo] ==
				gArray[dataArray[i].columnNo]) {
			//within the same group
			R_wSum += dataArray[i].dist;
			wCount++;
		}
		else {
			//between two different groups
			R_bSum += dataArray[i].dist;
			bCount++;
		}

	R_w = R_wSum/wCount;
	R_b = R_bSum/bCount;
	R = (R_b - R_w) / (data_size/2);

	return R;
}

void quickSort(dgPair arr[], int left, int right) {
	//quick sort adapted from
	//http://www.algolist.net/Algorithms/Sorting/Quicksort
	//key differences
	// -- data types
	// -- randomized pivot selection

	int i = left, j = right;
	dgPair tmp, pivot;

	srand ( time(NULL) );
	pivot = arr[left + rand() % (right-left+1)];

	/* partition */
	while (i <= j) {
		while (arr[i].dist < pivot.dist) i++;
		while (arr[j].dist > pivot.dist) j--;
		if (i <= j) {
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	};
	/* recursion */
	if (left < j)
		quickSort(arr, left, j);
	if (i < right)
		quickSort(arr, i, right);
}
