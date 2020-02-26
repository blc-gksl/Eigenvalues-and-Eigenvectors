#include<iostream>
#include<fstream>// using fstream library to read from and write on files
#include<string> //using string library to read lines and find n value of matrix in reading file
#include<cmath> //using cmath library to use abs function
#define Emaschine 0.0000001 //defining Emaschine to detect singularity
#include<iomanip> // to use fixed & setprecision to get the value in the form that exactly how we want
// do not forget to place sth like this (A.txt 0.000006 Output.txt)  to commend line --> Execute->parameters

using namespace std;

int findn(char *fileA,int n){
	int b=0;
	string line; 
  ifstream projectfileA(fileA);    //Making projectfileA to read fileA file
  if (projectfileA.is_open()){       //opening fileA file
    while ( getline (projectfileA,line) )  //counting lines to find n
    {
      b++;
    }
    projectfileA.close();  //closing fileA file
  }
  else { cout << "Unable to open file A (first try)"<<endl;} //to warn what if program fail to read txt
  n=b;
  //cout<<n<<endl;	// -->can be activated to check if n is controlled or not
  return n;
}

double** buildmatriceA(int n,char *fileA){
	//	 form a matrice by dynamic allocation
	int i,j;
	double **matriceX=new double*[n];// forming A matrice(matrice A)
	for(i=0;i<n;i++)
	matriceX[i]=new double[n];
    // reading fileA and fill formed matrices 'A' (matrisA[n][n])
	ifstream projectfileA(fileA);
	if(projectfileA.is_open()){

	for(i= 0; i < n; i++){	
    for(j = 0; j < n; j++){
    	projectfileA>>matriceX[i][j]; 	//get matrice's elements' value through reading files
    }}
	   projectfileA.close();	
	}
	else { cout << "Unable to open file A (second try)"<<endl;}
	
	return matriceX;
}

int singularitycheck(int n,char *fileA){
    //	 form two matrices by dynamic allocation
	int i,j;
	double **matriceA=new double*[n];// forming A matrice(matrice A)
	for(i=0;i<n;i++)
	matriceA[i]=new double[n];
    // reading fileA and fill formed matrices 'A' (matrisA[n][n])
	ifstream projectfileA(fileA);
	if(projectfileA.is_open()){

	for(i= 0; i < n; i++){	
    for(j = 0; j < n; j++){
    	projectfileA>>matriceA[i][j]; 	//get matrice's elements' value through reading files
	} }
	   projectfileA.close();	
	}
	else { cout << "Unable to open file A (third try)"<<endl;}
   
   //  making needed row exchanges while also doing partial pivoting
   int l,m;
   double temp,temp3; 
   for(i=0;i<n;i++){ 		
   	   for(l=0;l<n-1;l++){ // the for's  used to exchanges row properly
   	   for(j=i;j<n-1;j++){
   	   if(abs(matriceA[j+1][i])>abs(matriceA[j][i])){
   	   	    for(m=0;m<n;m++){
   	   	  	   temp=matriceA[j+1][m];
   	   	  	   matriceA[j+1][m]=matriceA[j][m];
   	   	  	   matriceA[j][m]=temp;	   
       }}}}
     // partial pivoting
   	    for(j=i+1;j<n;j++){
   			temp3=(matriceA[j][i])/(matriceA[i][i]);
   			for(m=i;m<n;m++){
   			matriceA[j][m]=matriceA[j][m]-temp3*matriceA[i][m];
   		}}}   	

   // checking if 'A' matrice (matriceA[n][n]) is singular or not by flagk variable
	int flagk=0;
	for(i=0;i<n;i++){
		if(abs(matriceA[n-1][i])>Emaschine) flagk=1;
	}
	
   int singularityflag=flagk;
   // deleting dynamically formed matriceA
   for(i=0;i<n;i++){     //first delete one of two dimensions then  other dimension
   	delete[] matriceA[i];
   }
   delete[] matriceA;
	return singularityflag;
}

class Matrix{
	private:
	int n;
	double **matriceA;
	double *eigenvector;
	double tolerance;
    public:
    // take the private(here could be protected too) value by Matrix constructor
    Matrix (int x,double **y,double *z,double t) : n(x), matriceA(y),eigenvector(z),tolerance(t) {};
    // by findEvector function find and return eigenvector
    double *findEvector(){
    double tempv[n];  //to update and find right value of eigenvector ,a temporary vector is used
    int i,j;
    double temp,d=0; //to find eigenvalue and update it d and a temporary variable temp is used
    //with dowhile loop make needed operations -->
    do{
    // first for loop multiply main matrice and current selected vector and form a temporary vector 
	// to get closer and closer to eigenvector and eigenvalue
	for(i=0;i<n;i++)
    {       tempv[i]=0;
            for(j=0;j<n;j++)
                tempv[i]+=matriceA[i][j]*eigenvector[j];        }
                
    for(i=0;i<n;i++)  eigenvector[i]=tempv[i]; // update eigenvector
    
    temp=d; //update current eigenvalue
    d=0;    //first reset d to use it to normalize current eigenvector and find new eigenvalue
    // this for loop is used to find biggest element in the vector
    for(i=0;i<n;i++)
    {       if(abs(eigenvector[i])>abs(d))
                d=eigenvector[i];                          }
                
    // and this for is used to normalization
	for(i=0;i<n;i++)  eigenvector[i]/=d;     
    }while(abs(d-temp)>tolerance); //check if we are close enough to the eigenvale

    return eigenvector; 
	}
	
	// the function below (findEvalue) is the same as the above one but above function return the eigen vector 
	//while the below one the eigen value
	double findEvalue(){
    double tempv[n];
    int i,j;
    double temp,d=0;
    do{
	for(i=0;i<n;i++)
    {       tempv[i]=0;
            for(j=0;j<n;j++)
                tempv[i]+=matriceA[i][j]*eigenvector[j];        }
                
    for(i=0;i<n;i++)  eigenvector[i]=tempv[i];
    
    temp=d;
    d=0;
    
    for(i=0;i<n;i++)
    {       if(abs(eigenvector[i])>abs(d))
                d=eigenvector[i];                          }
                
    
	for(i=0;i<n;i++)  eigenvector[i]/=d;     
    }while(abs(d-temp)>tolerance);	
    
	return d;
	}
};

class Matrix2{
	private:
	int n;
	double **matriceX;
    public:
    // take the private(here could be protected too) value by Matrix constructor
    Matrix2 (int x,double **y) : n(x), matriceX(y) {};
    // form matriceX  as an identity matrix by formidentitymatrix function
    double **formidentitymatrix(){
    	int i,j;
    	for(i=0;i<n;i++){
    		for(j=0;j<n;j++){
    			if(i==j){
    				matriceX[i][j]=1;
				}
				else matriceX[i][j]=0;
			}}
    	return matriceX;
	}
};

class Matrix3{
	private:
	int n;
	double **matriceA;
	double **matriceX;
    public:
    // take the private(here could be protected too) value by Matrix constructor
    Matrix3 (int x,double **y,double **z) : n(x), matriceA(y),matriceX(z) {};
    // form matriceX  as inverse matrix of matriceA by forminversematrix function
    // by gauss elimination  and with augmented matrix method find the inverse matrix
    double **forminversematrix(){
    	
   int i,j,l,m;
   double temp,temp2,temp3,temp4,temp5,temp6; 
   for(i=0;i<n;i++){ 		
   	   for(l=0;l<n-1;l++){ // the for's used to exchanges rows properly
   	   for(j=i;j<n-1;j++){
   	   if(abs(matriceA[j+1][i])>abs(matriceA[j][i])){
   	   	    for(m=0;m<n;m++){
   	   	  	   temp=matriceA[j+1][m];
   	   	  	   matriceA[j+1][m]=matriceA[j][m];
   	   	  	   matriceA[j][m]=temp;
   	   	  	   
			   temp4=matriceX[j+1][m];
   	   	  	   matriceX[j+1][m]=matriceX[j][m];
   	   	  	   matriceX[j][m]=temp4;	   
       }}}}
         
     //  partial pivoting
   	    for(j=i+1;j<n;j++){
   			temp3=(matriceA[j][i])/(matriceA[i][i]);
   			for(m=i;m<n;m++){
   			matriceA[j][m]=matriceA[j][m]-temp3*matriceA[i][m];}
   	        
   	       for(m=0;m<n;m++){
   			matriceX[j][m]=matriceX[j][m]-temp3*matriceX[i][m]; 
		
		}}}
	// full pivotng
    	for(i=n-1;i>0;i--){ 
    	for(j=i-1;j>(-1);j--){
   			temp5=(matriceA[j][i])/(matriceA[i][i]);
   			for(m=n-1;m>(-1);m--){
   			matriceA[j][m]=matriceA[j][m]-temp5*matriceA[i][m];
   	    }
   	        for(m=n-1;m>(-1);m--){
   			matriceX[j][m]=matriceX[j][m]-temp5*matriceX[i][m];
		
		}}}
		//
    	for(i=0;i<n;i++){ 
    	temp6 = matriceA[i][i];
    	for(j=0;j<n;j++){
    		matriceX[i][j] /= temp6;
    		matriceA[i][j] /= temp6;  //can be checked that matriceA is an identity matrix now
   		}   
		}
    	///
    	return matriceX;
	}
};

int main(int argc, char *argv[]){
   //First check if there are needed number of inputs
   if(argc<4){
   cout<<"ERROR type 1: Not enough command line arguments"<<endl; 
   return 0;
   }
   else if(argc>4){
   cout<<"ERROR type 2: Too many command line arguments"<<endl;
   return 0;
   }
   else {
   	
	char *fileA = argv[1];
	double tolerance = atof(argv[2]);
	char *outputfile =argv[3];
	// cout<<fileA<<endl<<tolerance<<endl<<outputfile<<endl;  //-->can be activated to check if names assigned or not
	 
	// reading fileA file and decide matrice's form (finding n's value)
	int n=0;
	n=findn(fileA,n);   //find the value of n with findn function
	
	int i;
	double **matriceA=new double*[n];// forming A matrice(matrice A)
	for(i=0;i<n;i++)
	matriceA[i]=new double[n];
	matriceA=buildmatriceA(n,fileA);   //find the matriceA with buildmatriceA function
	
	//defining a parameter to check the matrice if it is singular or not
	int singularityflag=0;
	singularityflag = singularitycheck(n,fileA); // singularity check function
	if(singularityflag==0){
	cout<<"ERROR type 3: The matrix is singular"<<endl; // if the matrix is singular give error message and quit
    return 0;
	}
	else{
	double maxeigenvalue; //form maximum eigenvalue
	// making eigenvector first form
	double *eigenvector =new double[n] ;
	eigenvector[0]= 1;                     //first element is 1
	for(i=1;i<n;i++)   eigenvector[i]= 0;  //other elements are 0
	
	Matrix MYmatrix(n,matriceA,eigenvector,tolerance);  // forming the class MYmatrix
	maxeigenvalue = MYmatrix.findEvalue();  // assign maximum eigenvalue
	eigenvector = MYmatrix.findEvector();	// assign maximum eigenvalue's eigenvector
		
	double **matriceX=new double*[n];// forming A matrice(matrice X) in the form of identity matrix 
	for(i=0;i<n;i++)
	matriceX[i]=new double[n];	
	
	//assign matriceX elements like identity matrix
	Matrix2 MYmatrix2(n,matriceX);
	matriceX= MYmatrix2.formidentitymatrix();
	
	//form inverse matrix
	Matrix3 MYmatrix3(n,matriceA,matriceX);
	matriceX= MYmatrix3.forminversematrix();	
		
	double mineigenvalue; //form minimum eigenvalue
	// making eigenvector2's first form
	double *eigenvector2 =new double[n] ;
	eigenvector2[0]= 1;                     //first element is 1
	for(i=1;i<n;i++)   eigenvector2[i]= 0;  //other elements are 0
	
	Matrix MYmatrix4(n,matriceX,eigenvector2,tolerance);  // forming the class MYmatrix4
	mineigenvalue = MYmatrix4.findEvalue();  // assign maximum eigenvalue
	//reverse maximum eigenvalue for A^-1 matrix to get minimum eigenvalue of A matrix
	mineigenvalue =1/(mineigenvalue);        
	eigenvector2 = MYmatrix4.findEvector();	// assign minimum eigenvalue's eigenvector
			
	/* writing eigenvalue and eigenvector into a txt file that is Output.txt
	   used fixed and setprecision to get desired form of the answer        */
   ofstream projectfilex(outputfile);
   if (projectfilex.is_open()){
   	// first largest eigenvalue and its vector
   	projectfilex<<"Eigenvalue#1: "<< fixed << setprecision(2) << maxeigenvalue<<endl;
    for(i = 0; i < n; i++){
	    projectfilex<< fixed << setprecision(2) << eigenvector[i]<<endl;
	}
	// then smallest eigenvalue and its vector
	projectfilex<<"Eigenvalue#2: "<< fixed << setprecision(2) << mineigenvalue<<endl;
    for(i = 0; i < n; i++){
	    projectfilex<< fixed << setprecision(2) << eigenvector2[i]<<endl;
	}
  }
   else { cout << "Unable to open file outputfile"<<endl;}
   
   delete[] eigenvector2;  // deleting dynamically formed eigenvector2
   // deleting dynamically formed matriceX
    for(i=0;i<n;i++){     //first delete one of two dimensions then  delete other dimension
   	delete[] matriceX[i];
    }
    delete[] matriceX;	
	delete[] eigenvector;  // deleting dynamically formed eigenvector
	}
    // deleting dynamically formed matriceA
    for(i=0;i<n;i++){     //first delete one of two dimensions then  delete other dimension
   	delete[] matriceA[i];
    }
    delete[] matriceA;
}
    return 0;
}
