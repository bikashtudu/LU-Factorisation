#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
double **mat,**P,**L,**U,**matTemp, *B,*Btemp, *x, *v;
double start=0,stop=0,res=0,ptotal=0;
void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}
void printMatrix(double **mat,int size)
{
	printf("\n \n");
	int i,j;
    for (i=0;i<size;i++) {
        for(j=0;j<size;j++) {
           	 if(mat[i][j]>=0)
           	 printf(" ");
 			 printf(" %f ",mat[i][j]);
			        }
        printf("\n");
    }
	printf("\n");
	}
	
void lu_pivot(int size)
{
	start=omp_get_wtime();	
    int k;
    for(k = 0; k < size-1; k++){
    	int maxi=k;
    	for(int i=k+1;i<size;i++)
    		if(abs(U[i][k]) > abs(U[maxi][k]))
    			maxi=i;
    	if(maxi!=k){
   
    	for(int i=k;i<size;i++)
			swap(&U[k][i],&U[maxi][i]);
	
		for(int i=0;i<=k-1;i++)
		 	swap(&L[k][i],&L[maxi][i]);
	
		for(int i=0;i<size;i++)
		 	swap(&P[k][i],&P[maxi][i]);
		 	}
		 if(U[k][k]==0)
		 {
		 	printf("\n INFINITE MANY SOLUTION\n");
		 	exit(0);
		 	}
	
      	for(int row = k+1 ; row < size; row++){
       		L[row][k] = U[row][k]/U[k][k];
       	
       		for(int col = k; col < size; col++)
         		U[row][col] = U[row][col] - L[row][k]*U[k][col];
        }
    }
    stop=omp_get_wtime();
    res=stop-start;
    printf("\nExecution Time for Partial Pivot LU Factorization = %f\n",res);
    ptotal+=res;
}

void lu_npivot(int size)
{
	start=omp_get_wtime();	
    int k;
    for(k = 0; k < size-1; k++){
    	 if(U[k][k]==0)
		 {
		 	printf("\nError...\n");
		 	exit(0);
		 	}
	
      	for(int row = k+1 ; row < size; row++){
       		L[row][k] = U[row][k]/U[k][k];
       	
       		for(int col = k; col < size; col++)
         		U[row][col] = U[row][col] - L[row][k]*U[k][col];
        }
    }
    stop=omp_get_wtime();
    res=stop-start;
    printf("\nExecution Time for No Pivot LU Factorization = %f\n",res);
    ptotal+=res;
}

void generate(int size)
{
	mat = malloc(size *(sizeof *mat));
	matTemp = malloc(size *(sizeof *mat));
	P = malloc(size * (sizeof *P));
	L = malloc(size *(sizeof *L));
	U = malloc(size *(sizeof *U));
 		int i =0, j=0;
    for(i=0;i<size;i++){
    	mat[i] = malloc((sizeof *mat[i]) * size);
    	matTemp[i] = malloc((sizeof *mat[i]) * size);
    	P[i] = malloc((sizeof *P[i]) * size);
    	L[i] = malloc((sizeof *L[i]) * size);
    	U[i] = malloc((sizeof *U[i]) * size);
        }
        
   	for(i=0;i<size;i++){
    	for( j=0; j< size; j++) {
            mat[i][j] = i+1;
            mat[j][i] = i+1;
            U[i][j]=mat[i][j];
            U[j][i]=mat[j][i];
            if(i==j)
            {
            	P[i][j]=1;
            	L[i][j]=1;
            	}
            else
            {
            	P[i][j]=0;
            	L[i][j]=0;
            	}
            matTemp[i][j]=0;
        }
        }
	x = malloc(size*(sizeof(double)));
	v = malloc(size*(sizeof(double)));
	B= malloc(size*(sizeof(double)));
	Btemp= malloc(size*(sizeof(double)));
	for(int i=0;i<size;i++)
		{
			B[i]=2*i+7;
			Btemp[i] = 0;
			}
  }
  
void printVector(double *t,int size)
{
	printf("\n\n");
	if(t==B){
	
	for(int i=0;i<size;i++)
		printf("B[%d]: %f\n",i,B[i]);
		}
	else if( t== v){
	
	for(int i=0;i<size;i++)
		printf("v[%d]: %f\n",i,v[i]);
		}
	else if(t== x)
	{
	for(int i=0;i<size;i++)
		printf("x[%d]: %f\n",i,x[i]);
		}
	else
	{
	for(int i=0;i<size;i++)
		printf("PB[%d]: %f\n",i,Btemp[i]);
		}
}
		
void SolveV(int size)
{
	start=omp_get_wtime();	
	for(int i=0;i<size;i++)
	{
		v[i]=B[i];
		double sub=0;
		for(int j=0;j<i;j++)
			sub+=L[i][j]*v[j];
		v[i]-=sub;
    }
    stop=omp_get_wtime();
    res=stop-start;
    printf("\nExecution Time for Evaluating 'V' = %f\n",res);
    ptotal+=res;
    }
    
    
    
void SolveX(int size)
{
	start=omp_get_wtime();	
 	for(int i=size-1;i>=0;i--)
    {
    	x[i]=v[i];
    	double sum=0;
		for(int j=i+1;j<size;j++)
			sum+=U[i][j]*x[j];
		x[i]-=sum;
		x[i]/=U[i][i];
    }
    stop=omp_get_wtime();
    res=stop-start;
    printf("\nExecution Time for Evaluating 'X' = %f\n",res);
    ptotal+=res;
}
int main()
{
 	int size,choice;
 	printf("Enter The Square Matrix Size: ");
 	scanf("%d",&size);
 	generate(size);
    //printf("\nDone:\n");
	int i,j;
	printf("\n A");
   	printMatrix(mat,size);
   	printVector(B,size);
   	//printMatrix(L,size);
   	//printMatrix(U,size);
   	printf("\n Partial Pivoted ....1");
   	printf("\n No Pivoted..........2\n");
    scanf("%d",&choice);
    if(choice==1){
		lu_pivot(size);
	//printMatrix(mat,size);
		printf("\n LU PARTIAL PIVOTED\n");
		printf("\n L");
		printMatrix(L,size);
		printf("\n U");
		printMatrix(U,size);
		printf("\n P");
		printMatrix(P,size);
		for(int i=0;i<size;i++)
			for(int k=0;k<size;k++)
				Btemp[i]+=P[i][k]*B[k];
		for(int i=0;i<size;i++)
			for(int j=0;j<size;j++)
				for(int k=0;k<size;k++)
					matTemp[i][j]+=P[i][k]*mat[k][j];
		printf("\n PA");
   		printMatrix(matTemp,size);
		printVector(Btemp,size);
		B = Btemp;
		SolveV(size);
		printVector(v,size);
		SolveX(size);
		printVector(x,size);
		}
	else if(choice ==2){
		printf("\n LU \n");
		lu_npivot(size);
		printf("\n L");
		printMatrix(L,size);
		printf("\n U");
		printMatrix(U,size);
		SolveV(size);
		printVector(v,size);
		SolveX(size);
		printVector(x,size);
		}
		printf("\n\nTotal Time of Execution = %f",ptotal);
}
