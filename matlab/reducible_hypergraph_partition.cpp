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

// len(W_triu) = len(W_tril) = len(ix) = len(jx) = m
double eval_Q1_digraph(double * y, double * W_triu, double * W_tril, int * ix, int * jx, int m)
{
    double Q1 = 0;
    int u, v;
    for(int cnt = 0; cnt < m; cnt++)
    {
        u = ix[cnt];
        v = jx[cnt];
        Q1 += W_triu[cnt] * fmax(y[u] - y[v], 0) + W_tril[cnt] * fmax(y[v] - y[u], 0);
    }
    return Q1;
}

// len(fA) = n, len(alpha) = m
void f_A(double * fA, double * alpha, double * W_triu, double * W_tril, int * ix, int * jx, int n, int m)
{
    for(int i = 0; i < n; i++)
        fA[i] = 0;
    
    int u, v;
    double dummy;
    for(int cnt = 0; cnt < m; cnt++)
    {
        u = ix[cnt];
        v = jx[cnt];
        dummy = (W_triu[cnt] + W_tril[cnt]) * alpha[cnt] - W_tril[cnt];
        fA[u] += dummy;
        fA[v] -= dummy;
    }
}

// len(g) = N
void inner_problem_single_step(double * alpha, double * beta, double * t, double L, double * g, double * W_triu, double * W_tril, int * ix, int * jx, int N, int n, int m)
{
    double * alpha_old = new double [m];
    double * fA = new double [n];
    double * D = new double [n];
    
    copy(alpha_old, alpha, m);
    
    f_A(fA, beta, W_triu, W_tril, ix, jx, n, m);
    
    copy(D, fA, n);
    for(int i = 0; i < N; i++)
        D[i] -= g[i];
    
    int u, v;
    double alphacur;
    for(int cnt = 0; cnt < m; cnt++)
    {
        u = ix[cnt];
        v = jx[cnt];
        alphacur = beta[cnt] - 2 * (W_triu[cnt] + W_tril[cnt]) * (D[u] - D[v]) / L;
        alpha[cnt] = fmax(fmin(alphacur, 1), 0);
    }
    
    double t_old = *t;
    *t = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
    
    for(int cnt = 0; cnt < m; cnt++)
        beta[cnt] = alpha[cnt] + (t_old - 1) / *t * (alpha[cnt] - alpha_old[cnt]);
    
    delete[] alpha_old;
    delete[] fA;
    delete[] D;
}


void main_func(double * labels, double * NCut, double * x_final,
       int ** incidence_list, double ** parameter_list, double * mu,
       int N, int R, int * esize_list, int mode, double delta,
       double dec_outloop, double err_inloop, double * warmstart,
       double L, double * W_triu, double * W_tril, int * ix, int * jx, int n, int m, int maxiter_inner_bound, int maxiter_inner, int maxiter_outer)
{
    double Q1, tempeta, eta, normy, gap;
    double tempgap = 0;
    double t;
    int u, v;

    double * x = new double [N];
    double * g = new double [N];
    double * y = new double [n]();
    double * fA = new double [n];
    double * alpha = new double [m];
    double * beta = new double [m];
    
    copy(x, warmstart, N);
    remove_median(x, mu, N);
    Q1 = eval_Q1(x, incidence_list, parameter_list, R, esize_list, mode, delta);
    tempeta = eval_R1(Q1, x, mu, N);
    
    copy(y, x, N);
    for(int cnt = 0; cnt < m; cnt++)
    {
        u = ix[cnt];
        v = jx[cnt];
        if(y[u] > y[v])
            alpha[cnt] = 1;
        else
            alpha[cnt] = 0;
    }

    for(int iter_outer = 0; iter_outer < maxiter_outer; iter_outer++)
    {
        std::cout << " iter_outer = " << iter_outer << std::endl;
       
        derivative_mu_norm(g, x, mu, N);
        eta = tempeta;
        for(int i = 0; i < N; i++)
            g[i] *= eta;
        
        copy(beta, alpha, m);
        t = 1;
        for(int iter_inner = 0; iter_inner < maxiter_inner_bound; iter_inner++)
        {
            inner_problem_single_step(alpha, beta, &t, L, g, W_triu, W_tril, ix, jx, N, n, m);
           
            if((iter_inner+1)%maxiter_inner == 0)
            {
                f_A(fA, alpha, W_triu, W_tril, ix, jx, n, m);
                copy(y, fA, n);
                for(int i = 0; i < N; i++)
                    y[i] -= g[i];
                normy = sqrt(square_sum(y, n));
                for(int i = 0; i < n; i++)
                    y[i] /= -normy;
                gap = tempgap;
                tempgap = abs(eval_Q1_digraph(y, W_triu, W_tril, ix, jx, m) - inner(y, g, N) + normy);
                
                copy(x, y, N);
                remove_median(x, mu, N);
                Q1 = eval_Q1(x, incidence_list, parameter_list, R, esize_list, mode, delta);
                tempeta = eval_R1(Q1, x, mu, N);
                
                std::cout << " iter = " << iter_inner << " gap = " << tempgap << " eta = " << tempeta << std::endl;
                
                if(tempgap < err_inloop || iter_inner == maxiter_inner_bound - 1)
                {
                   if(tempeta < eta)
                   {
                       std::cout << "update x_final" << std::endl;
                       copy(x_final, x, N);
                   }
                   break;
                }
           }
            
        }
        if(tempeta > (1-dec_outloop)*eta)
            break;
        
    }

    optthreshold(labels, NCut, x_final, incidence_list, parameter_list, mu, N, R, esize_list, mode, delta);
    
    delete[] x;
    delete[] g;
    delete[] y;
    delete[] fA;
    delete[] alpha;
    delete[] beta;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 3 || nrhs != 20)
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
        
    double L = *(mxGetPr(prhs[10]));
    double * W_triu = mxGetPr(prhs[11]);
    double * W_tril = mxGetPr(prhs[12]);
    int n = *(mxGetPr(prhs[15]));
    int m = *(mxGetPr(prhs[16]));
    int maxiter_inner_bound = *(mxGetPr(prhs[17]));
    int maxiter_inner = *(mxGetPr(prhs[18]));
    int maxiter_outer = *(mxGetPr(prhs[19]));
    
    double * ix_org = mxGetPr(prhs[13]);
    double * jx_org = mxGetPr(prhs[14]);
    int * ix = new int [m];
    int * jx = new int [m];
    for(int k = 0; k < m; k++)
    {
        ix[k] = (int) ix_org[k];
        jx[k] = (int) jx_org[k];
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
              L, W_triu, W_tril, ix, jx, n, m, maxiter_inner_bound, maxiter_inner, maxiter_outer);
   
}

