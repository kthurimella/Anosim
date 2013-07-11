This code was compiled on 64-bit Linux Ubuntu OS, using the G++ 
compiler that comes shipped with it.

To compile the source code:
% g++ similar.cpp -o anosim

To compute the similarity index use it as

anosim -d <sample_filename> 
       -g <group_label_filename>  
       [-p no-of-permutations]  

example: 
% ./anosim -d test.txt -g gtest.txt -p 3

It expects a minimum of two switches pointing to the files that contain
sample data and the group labels.

The sample data is assumed contain integers 1..n in the first row as well
as the first column, where n is the number of samples.  These descriptor
values are ignored while populating the data matrix.

You can optionally specify the number of permutations to be used in
computing the p-value with -p switch which is set to the number of
permutations variable NO_OF_PERMS.  The code permutes the group labels
MIN(factorial(n), NO_OF_PERMS).  The default value for NO_OF_PERMS is
set to 1 million. (Up to 2 million permutations can be used
without a noticeable delay.)

TO MAKE THE CODE ROBUST: 
-----------------------
1) Test against random input values
2) If speed becomes an issue, 
   a) replace recursion in quicksort with
      iteration
   b) since the code that is executed multiple times is
      anosim_stat function, replace the function call and compute
      R in-place.
3) Add exception handling
4) Check data for integrity, e.g. see if the input matrix has
   only real values and it is a symmetric matrix.
5) Make it object-oriented, turn it into a library, make it suitable for
   distribution. 
