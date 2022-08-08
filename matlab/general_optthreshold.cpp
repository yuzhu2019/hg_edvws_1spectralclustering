#include <stdio.h>
#include <math.h>
#include <random>
#include <string>
#include <iostream>
#include "mex.h"

void swap_double(double * a, int p, int q)
{
    double temp = a[p];
    a[p] = a[q];
    a[q] = temp;
}

double sum(double * data, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += data[i]; }
    return res;
}

void MYsort_descend(double * a, double * b, int len)
{
    // sort a, permutate b according to a
    if(len <= 1) return;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, len-1);
    int pos = distribution(generator), p = 0;
    double pivot = a[pos];
    swap_double(a, pos, len-1);
    swap_double(b, pos, len-1);
    for(int i = 0; i < len-1; i++)
    {
        if(a[i] > pivot)
        {
            swap_double(a, i, p);
            swap_double(b, i, p);
            p++;
        }
    }
    swap_double(a, p, len-1);
    swap_double(b, p, len-1);
    MYsort_descend(a, b, p);
    MYsort_descend(a+p+1, b+p+1, len-p-1);
}

void edvw2para_e_part(double * subw, double * para, int len, int mode, double delta)
{
    // para has length len
    // subw has length len+1
    subw[0] = 0;    // the first entry
    subw[len] = 0;  // the last entry
    double total, ss;
    total = sum(para, len);
    for(int i = 1; i < len; i++)
    {
        ss = sum(para, i);
        switch(mode)
        {
            case 0: // c
                subw[i] = ss * (total - ss);
                break;
            case 1: // s
                subw[i] = std::min(ss, total - ss);
                break;
            case 2: // l
                subw[i] = std::min({ss, total - ss, delta * total});
                break;
            case 3: // ct
                subw[i] = std::min(ss * (total - ss), delta * total * total);
                break;
            default:
                throw std::invalid_argument( "specify a mode: {0: c, 1: s, 2: l, 3: ct}" );
        }
    }
}

double eval_Q1(double * x, int ** incidence_list, double ** parameter_list, int R, int * esize_list, int mode, double delta)
{
   double Q1 = 0;
   int esize;
   for(int i = 0; i < R; i++) // for each hyperedge
   {
       esize = esize_list[i];
       double * tempx = new double [esize];
       double * para = new double [esize];
       double * subw = new double [esize+1];
       for(int j = 0; j < esize; j++)
       {
           tempx[j] = x[incidence_list[i][j]];
           para[j] = parameter_list[i][j];
       }
       MYsort_descend(tempx, para, esize);
       edvw2para_e_part(subw, para, esize, mode, delta);
       for(int j = 0; j < esize; j++)
           Q1 += tempx[j] * (subw[j+1] - subw[j]);
       delete[] tempx;
       delete[] para;
       delete[] subw;
   }
   return Q1;
}

void optthreshold(double * labels, double * NCut,
                 double * x, int ** incidence_list, double ** parameter_list, double * mu,
                 int N, int R, int * esize_list, int mode, double delta)
{
   double * index = new double [N];
   double * xsort = new double [N];
   for(int i = 0; i < N; i++)
   {
       index[i] = i;
       xsort[i] = x[i];
   }
   MYsort_descend(xsort, index, N);
   int optpos = 0;
   double tempvol = 0, tempval, Q1;
   double * xtemp = new double [N](); // initialize all zeros
   *NCut = INFINITY;
   double sum_mu = sum(mu, N);
   for(int i = 0; i < N-1; i++)
   {
       xtemp[(int)index[i]] = 1;
       Q1 = eval_Q1(xtemp, incidence_list, parameter_list, R, esize_list, mode, delta);
       tempvol += mu[(int)index[i]];
       tempval = Q1/fmin(tempvol, sum_mu - tempvol);
       if(tempval < *NCut)
       {
           *NCut = tempval;
           optpos = i;
       }
   }
   for(int i = 0; i < N; i++)
   {
       if(i <= optpos)
           labels[(int)index[i]] = 1;
       else
           labels[(int)index[i]] = 0;
   }
   delete[] index;
   delete[] xsort;
   delete[] xtemp;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 2 || nrhs != 8)
    {
        mexWarnMsgTxt("Check Parameters");
        return;
    }
    
    // INPUT
    
    double * mu = mxGetPr(prhs[2]);
    int N = *(mxGetPr(prhs[3]));
    int R = *(mxGetPr(prhs[4]));
    int mode = *(mxGetPr(prhs[5]));
    double delta = *(mxGetPr(prhs[6]));
    double * x = mxGetPr(prhs[7]);
            
    
    // read incidence_list
    const mxArray * incidence_list_org = prhs[0];
    mxArray * incidence_Element;
    int ** incidence_list = new int * [R];
    int * esize_list = new int [R];
    double * templist;
    for(int j = 0; j < R; j++)
    {
       incidence_Element = mxGetCell(incidence_list_org, j);
       esize_list[j] = (int)mxGetN(incidence_Element);
       incidence_list[j] = new int [esize_list[j]];
       templist = mxGetPr(incidence_Element);
       for(int k = 0; k < esize_list[j]; k++)
           incidence_list[j][k] = (int)templist[k]-1; // notice!
    }
    
    // read parameter_list
    const mxArray * parameter_list_org = prhs[1];
    mxArray * parameter_Element;
    double ** parameter_list = new double * [R];
    for(int j = 0; j < R; j++)
    {
       parameter_Element = mxGetCell(parameter_list_org, j);
       parameter_list[j] = new double [esize_list[j]];
       templist = mxGetPr(parameter_Element);
       for(int k = 0; k < esize_list[j]; k++)
           parameter_list[j][k] = templist[k];
    }

    // OUTPUT
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * labels = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * NCut = mxGetPr(plhs[1]);
    
    optthreshold(labels, NCut,
                 x, incidence_list, parameter_list, mu,
                 N, R, esize_list, mode, delta);
    
}

