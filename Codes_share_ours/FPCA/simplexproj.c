#include "mex.h"
// 使用MEX必须包含的头文件
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#define datatype double /* type of the elements in y */


static void simplexproj_Condat(datatype* y, datatype* x, const unsigned int length, const double a); 

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){

    if (nrhs != 1)
        mexErrMsgTxt("Wrong number of input arguments.\n");

    if (nlhs > 1)   
        mexErrMsgTxt("Too many output argumnents.\n");
    
    #define y_IN  prhs[0]
    #define x_OUT plhs[0]
    
    int M = mxGetM(y_IN);
    int N = mxGetN(y_IN);

    x_OUT = mxCreateDoubleMatrix(M, N, mxREAL);
    
    double *y = mxGetPr(y_IN);
    double *x = mxGetPr(x_OUT);
    
    simplexproj_Condat( y, x, M, 1);
}



static void simplexproj_Condat(datatype* y, datatype* x, 
const unsigned int length, const double a) {
	datatype*	aux = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	datatype*  aux0=aux;
	int		auxlength=1; 
	int		auxlengthold=-1;	
	double	tau=(*aux=*y)-a;
	int 	i=1;
	for (; i<length; i++) 
		if (y[i]>tau) {
			if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))
			<=y[i]-a) {
				tau=y[i]-a;
				auxlengthold=auxlength-1;
			}
			auxlength++;
		} 
	if (auxlengthold>=0) {
		auxlength-=++auxlengthold;
		aux+=auxlengthold;
		while (--auxlengthold>=0) 
			if (aux0[auxlengthold]>tau) 
				tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
	}
	do {
		auxlengthold=auxlength-1;
		for (i=auxlength=0; i<=auxlengthold; i++)
			if (aux[i]>tau) 
				aux[auxlength++]=aux[i];	
			else 
				tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
	} while (auxlength<=auxlengthold);
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
	if (x==y) free(aux0);
} 
