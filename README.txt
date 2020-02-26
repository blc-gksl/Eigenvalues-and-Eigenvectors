The code take three commend line arguments =>
1-) A matrix file  
2-) Tolerance value to check if we are close enough to find solution
3-) Another file to create and write the solution in it

So in order to use source.cpp file you need to have a matrix file and it will give you another file that contain eigenvalues(largest and smallest ones) and their solutions

The code file name is source with .cpp extension(source.cpp)


--> PROCESS STEPS:
1 - First read command line arguments and decide whether they are desireable input
2 - Take the arguments and assign them to some usable variables
3 - Creating a matrix by using our matrix file info and check if its singular or not
4 - Creating a vector to use in power iteration method
5 - Finding largest eigenvalue and its eigenvector
6 - To find smallest eigenvalue and its vector creating a matrix
7 - By Gaussian elimination method making last created matrix into inverse matrix of main matrix, to find smallest eigenvalue
8 - Form another vector and value and make them the smallest eigenvalue and eigenvector
9 - Write the solution into a file and delete the eigenvectors and matrices that are created dynamically
