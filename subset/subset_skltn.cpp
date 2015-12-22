#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>

using namespace std;

vector<double> S; // this is the input vector
double global_min = 1000000000000.0;
double main_sum;

//index is the index of the new element to be added into S
void GenerateSubset(int index, vector<double> &S, vector<double> subset) {

    if (index == S.size) { //then this subset is done
        update if it is a better solution than what we have so far
            return
    }
    
    // this is only a pseudocode!
    // you will need to make the appropriate modifications
    subset = subset U S(index)          //add the index element in subset
    GenerateSubset(index+1, S, subset)  //generate subsets with the index element in them
    remove the added element from the subset
    GenerateSubset(index+1, S, subset)  //generate subsets without the index element
}



int main(int argc, char* argv[]) {
  double start_time, end_time, time_diff;
  int i, nweights, nthreads;
  vector<double> subset; // this is where you will generate the subsets
  subset.clear();
 
  // Get the number of weights and number of threads from the command line
  if (argc == 3) {
    nweights = strtol(argv[1], NULL, 10);
    nthreads = strtol(argv[2], NULL, 10);
  }
  else {
    printf("\nWhen running this program, please include number of weights and number of threads on command line.\n");
    return 0;
  }
  
  printf("\nnumber of weigts: %d\n", nweights);
  printf("number of threads: %d\n", nthreads);
    
  main_sum = 0;
  //srand(time(NULL));
  srand(0);
  for (i = 0 ; i < nweights; i++) {
    S.push_back(((double)rand())/RAND_MAX);
    main_sum += S[i];
  }
  
  printf("main set : ");
  for (i = 0 ; i < S.size() ; i++)
    printf("%lf ", S[i]);
  printf("\n");
  printf("main sum = %lf\n", main_sum);
  
  start_time = omp_get_wtime();

  GenerateSubset(0, S, subset);
  
  end_time = omp_get_wtime();
  
  printf("\n minimum diff = %.14lf\n", global_min);
  
  printf("time needed = %f\n", end_time - start_time);
}
