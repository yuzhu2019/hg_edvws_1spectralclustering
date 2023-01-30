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

void copy(double * a, double * b, int len)
{
    // a <- b
    for(int i = 0; i < len; i++) { a[i] = b[i]; }
}

double sum(double * data, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += data[i]; }
    return res;
}

double square_sum(double * data, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += data[i] * data[i]; }
    return res;
}

double inner(double * a, double * b, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += a[i] * b[i]; }
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

void remove_median(double * x, double * mu, int len)
{
    // return a non-constant vector x satisfying 0 ∈ argmin_c || x - c1 ||_{l1, mu}
    int i, pos=0;
    double * xsort = new double [len];
    double * mu_index = new double [len];
    copy(xsort, x, len);
    copy(mu_index, mu, len);
    MYsort_descend(xsort, mu_index, len);
    for(i = 1; i < len; i++)
        mu_index[i] += mu_index[i-1];
    for(i = 0; i < len; i++)
    {
        if(mu_index[i] >= mu_index[len-1]/2)
        {
            pos = i;
            break;
        }
    }
    for(i = 0; i < len; i++)
        x[i] -= xsort[pos];
    delete[] xsort;
    delete[] mu_index;
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

double eval_R1(double Q1, double * x, double * mu, int N)
{
   double temp = 0;
   for(int i = 0; i < N; i++)
       temp += abs(x[i]) * mu[i];
   double R1 = Q1/temp;
   return R1;
}

void derivative_mu_norm(double * g, double * x, double * mu, int N)
{
   double epsilon = 1e-10;
   double mu_x_pos = 0;
   double mu_x_neg = 0;
   double mu_x_zeros = 0;
   for(int i = 0; i < N; i++)
   {
       if(x[i] > epsilon)
           mu_x_pos += mu[i];
       else if(x[i] < -epsilon)
           mu_x_neg += mu[i];
       else
           mu_x_zeros += mu[i];
   }
   double muzero = (mu_x_neg - mu_x_pos)/mu_x_zeros;
   for(int i = 0; i < N; i++)
   {
       if(x[i] > epsilon)
           g[i] = mu[i];
       else if(x[i] < -epsilon)
           g[i] = -mu[i];
       else
           g[i] = muzero * mu[i];
   }
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

// len(y) = len(ybar) = n, len(z) = len(edge_weights) = len(iy) = len(jy) = M, len(g) = N
void inner_problem_single_step_pdhg(double * y, double * z, double * ybar, double * sigma, double * tau, double * g, double * edge_weights, int * iy, int * jy, int N, int n, int M)
{
    int u, v;
    double tempval;
    double * y_old = new double [n];
    // update z
    for(int i = 0; i < M; i++)  // iterate over edges u -> v
    {
        u = iy[i];
        v = jy[i];
        tempval = z[i] + *sigma * edge_weights[i] * (ybar[u] - ybar[v]);
        z[i] = fmax(fmin(tempval, 1), 0);
    }
    // update y
    copy(y_old, y, n);
    for(int i = 0; i < M; i++)
    {
        u = iy[i];
        v = jy[i];
        tempval = *tau * edge_weights[i] * z[i];
        y[u] -= tempval;
        y[v] += tempval;
    }
    for(int i = 0; i < N; i++)
    {
        y[i] += *tau * g[i];
        y[i] /= 1 + *tau;
    }
    for(int i = N; i < n; i++)
        y[i] /= 1 + *tau;
    // update sigma and tau
    double theta = 1 / sqrt(1 + *tau);
    *tau *= theta;
    *sigma /= theta;
    // update ybar
    for(int i = 0; i < n; i++)
        ybar[i] = y[i] + theta * (y[i] - y_old[i]);
    delete[] y_old;
}

double eval_primal_obj(double * y, double * g, double * edge_weights, int * iy, int * jy, int N, int n, int M)
{
    // compute the objective of the primal problem f(y) + g(By)
    double primal_obj = 0;
    double tempval;
    int u, v;
    // compute f(y) = 0.5 * ||y - g||^2
    for(int i = 0; i < N; i++)
    {
        tempval = y[i] - g[i];
        primal_obj += tempval * tempval;
    }
    for(int i = N; i < n; i++)
        primal_obj += y[i] * y[i];
    primal_obj /= 2;
    // compute g(By)
    for(int i = 0; i < M; i++)  // iterate over edges u -> v
    {
        u = iy[i];
        v = jy[i];
        tempval = edge_weights[i] * (y[u] - y[v]);
        primal_obj += fmax(tempval, 0);
    }
    return primal_obj;
}

double eval_dual_obj(double * z, double * g, double * edge_weights, int * iy, int * jy, int N, int n, int M)
{
    // compute the objective of the dual problem (0 <= z <= 1)
    double * nBTz = new double [n](); // -B^T * z
    double tempval;
    int u, v;
    for(int i = 0; i < M; i++)
    {
        u = iy[i];
        v = jy[i];
        tempval = edge_weights[i] * z[i];
        nBTz[u] -= tempval;
        nBTz[v] += tempval;
    }
    double dual_obj = -square_sum(nBTz, n)/2 - inner(nBTz, g, N);
    delete[] nBTz;
    return dual_obj;
}

void main_func(double * labels, double * NCut, double * x_final,
       int ** incidence_list, double ** parameter_list, double * mu,
       int N, int R, int * esize_list, int mode, double delta,
       double dec_outloop, double err_inloop, double * warmstart,
       double B, double * edge_weights, int * iy, int * jy, int n, int M,
       int maxiter_inner_bound, int maxiter_inner, int maxiter_outer)
{
    double Q1, tempeta, eta, gap, primal_obj, dual_obj, sigma, tau;
    
    double * x = new double [N];
    double * g = new double [N];
    double * y = new double [n];
    double * z = new double [M];
    double * ybar = new double [n];
    
    // initialization of the outer-loop problem
    copy(x, warmstart, N);
    remove_median(x, mu, N);
    Q1 = eval_Q1(x, incidence_list, parameter_list, R, esize_list, mode, delta);
    tempeta = eval_R1(Q1, x, mu, N);
    
    // outer-loop (IPM)
    for(int iter_outer = 0; iter_outer < maxiter_outer; iter_outer++)
    {
        std::cout << " iter_outer = " << iter_outer << std::endl;
       
        // update g
        derivative_mu_norm(g, x, mu, N);
        eta = tempeta;
        for(int i = 0; i < N; i++)
            g[i] *= eta;
        
        // initialization of the inner-loop problem
        sigma = 1/B;
        tau = 1/B;
        copy(y, x, N);
        for(int i = N; i < n; i++)
            y[i] = 0;
        copy(ybar, y, n);
        for(int i = 0; i < M; i++)
            z[i] = 0;
        
        // inner-loop (PDHG)
        for(int iter_inner = 0; iter_inner < maxiter_inner_bound; iter_inner++)
        {
            inner_problem_single_step_pdhg(y, z, ybar, &sigma, &tau, g, edge_weights, iy, jy, N, n, M);
            
            if((iter_inner+1)%maxiter_inner == 0)
            {
                // compute the dual gap
                primal_obj = eval_primal_obj(y, g, edge_weights, iy, jy, N, n, M);
                dual_obj = eval_dual_obj(z, g, edge_weights, iy, jy, N, n, M);
                gap = abs(primal_obj - dual_obj);
                
                copy(x, y, N);
                remove_median(x, mu, N);
                Q1 = eval_Q1(x, incidence_list, parameter_list, R, esize_list, mode, delta);
                tempeta = eval_R1(Q1, x, mu, N);
                
                std::cout << " iter = " << iter_inner << " gap = " << gap << " eta = " << tempeta << std::endl;
                
                if(gap < err_inloop || iter_inner == maxiter_inner_bound - 1)
                {
                   if(tempeta < eta)
                   {
                       std::cout << "update x_final" << std::endl;
                       copy(x_final, x, N);
                   }
                   break; // end the inner-loop
                }
            }
        }
        if(tempeta > (1-dec_outloop)*eta)
            break; // end the outer-loop
    }

    optthreshold(labels, NCut, x_final, incidence_list, parameter_list, mu, N, R, esize_list, mode, delta);
    
    delete[] x;
    delete[] g;
    delete[] y;
    delete[] z;
    delete[] ybar;
}

// N: number of nodes in the hypergraph
// R: number of hyperedges in the hypergraph
// n: number of nodes in the reduced graph (augmented vertex set)
// M: number of (directed) edges in the reduced graph
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 3 || nrhs != 19)
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
    double dec_outloop = *(mxGetPr(prhs[7]));
    double err_inloop = *(mxGetPr(prhs[8]));
    double * warmstart = mxGetPr(prhs[9]);
       
    double B = *(mxGetPr(prhs[10]));
    double * edge_weights = mxGetPr(prhs[11]);
    int n = *(mxGetPr(prhs[14]));
    int M = *(mxGetPr(prhs[15]));
    int maxiter_inner_bound = *(mxGetPr(prhs[16]));
    int maxiter_inner = *(mxGetPr(prhs[17]));
    int maxiter_outer = *(mxGetPr(prhs[18]));
    
    double * iy_org = mxGetPr(prhs[12]);
    double * jy_org = mxGetPr(prhs[13]);
    int * iy = new int [M];
    int * jy = new int [M];
    for(int k = 0; k < M; k++)
    {
        iy[k] = (int) iy_org[k];
        jy[k] = (int) jy_org[k];
    }
    
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
    plhs[2] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * x_final = mxGetPr(plhs[2]);
    
    main_func(labels, NCut, x_final,
              incidence_list, parameter_list, mu,
              N, R, esize_list, mode, delta,
              dec_outloop, err_inloop, warmstart,
              B, edge_weights, iy, jy, n, M,
              maxiter_inner_bound, maxiter_inner, maxiter_outer);
   
}
