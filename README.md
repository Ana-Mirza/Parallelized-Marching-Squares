Name: Mirza Ana-Maria

Group: 331CA


### Homework 1 - Marching Squares -
The scope of this homework was parallelizing the sequencial code of the marching squares algorithm.

### Solution
In order to speed up the processing of images in this algorithm, I started with parallelizing the scaling process bringing images to a dimension of 2048 x 2048, which was the costliest operation. Next, I parallelized the two steps of the marching squares algorithm by dividing the matrices' lines between the threads, so that each thread had an individual number of lines to work with.

### Scalability
The parallel version of this algorithm uses efficiently resources by allocating the number of threads requested at the beginning and dividing the workload between the available threads.
