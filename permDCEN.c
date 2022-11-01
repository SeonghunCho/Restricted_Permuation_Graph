
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<math.h>
#include<Rmath.h>

void calc_cdf_value(double *x, int *M, double *cdf_x, double *value)
{
	int m = *M;
	int i, ind = m - 1;
	double xval = *x;
	if (xval < cdf_x[0])
	{
		*value = 0.0;
		return;
	}
	if (xval >= cdf_x[ 2*m - 1])
	{
		*value = 1.0;
		return;
	}
	for (i = 0; i < (m - 1); i++)
	{
		if ( (xval >= cdf_x[i]) & (xval < cdf_x[i + 1]))
		{
			ind = i;
			break;
		}
	}

	if (cdf_x[ind] == cdf_x[m+ind])
	{
		*value = cdf_x[3 * m + ind];
		return;
	}
	if ( (cdf_x[ind]< cdf_x[m+ind] )& (ind > 0))
	{
		*value = cdf_x[3 * m + ind - 1] + (xval - cdf_x[ind]) / (cdf_x[m + ind] - cdf_x[ind]) * cdf_x[2 * m + ind];
		return;
	}
	if ( (cdf_x[ind]< cdf_x[m+ind]) & (ind == 0))
	{
		*value = (xval - cdf_x[ind]) / (cdf_x[m + ind] - cdf_x[ind]) * cdf_x[2 * m + ind];
		return;
	}
}


void case2_pair(int *i, int *j, int *N, double *dat, int *M, double *cdf_dat, double *res)
{
	int n = *N, m= *M;
	double Li, Lj, Ri, Rj;
	double FLi, FLj, FRi, FRj;
	double val = 0;

	Li = dat[*i];
	Lj = dat[*j];
	Ri = dat[n + *i];
	Rj = dat[n + *j];

	calc_cdf_value(&Li, &m, cdf_dat, &FLi);
	calc_cdf_value(&Lj, &m, cdf_dat, &FLj);

	calc_cdf_value(&Ri, &m, cdf_dat, &FRi);
	calc_cdf_value(&Rj, &m, cdf_dat, &FRj);
	
	if ( (Li <= Ri) & (Ri < Lj) & (Lj <= Rj))
		val = -1;

	if ( (Li < Lj) & (Lj < Ri) & (Ri < Rj) & (FRj != FLj) & (FRi != FLi))
		val = pow((FRi - FLj), 2.0) / ((FRj - FLj) * (FRi - FLi)) - 1;

	if ( (Li < Lj) & (Lj < Rj) & (Rj < Ri) & (FRi != FLi))
		val = (2 * FRi - FRj - FLj) / (FRi - FLi) - 1;

	if ((Lj <= Rj) & (Rj < Li) & (Li <= Ri))
		val = 1;

	if ((Lj < Li) & (Li < Rj) & (Rj < Ri) & (FRj != FLj) & (FRi != FLi))
		val = 1 - pow((FRj - FLi),2.0) / ((FRj - FLj) * (FRi - FLi));

	if ((Lj < Li) & (Li < Ri) & (Ri < Rj) & (FRj != FLj))
		val = (FRi + FLi - 2 * FLj) / (FRj - FLj) - 1;

	*res = val;

	return;

}

void tau_est(int *N, double *dat_x, double *dat_y, int *Mx, int *My, double *cdf_x, double *cdf_y, double *tau)
{
	int n = *N, mx = *Mx, my = *My;
	int i, j;
	double sum_ab, sum_aa, sum_bb;

	size_t n2 = n * n * sizeof(double);

	double *a_mat, *b_mat;
	a_mat = (double *)malloc(n2);
	b_mat = (double *)malloc(n2);

	memset(a_mat, 0, n2);
	memset(b_mat, 0, n2);

	sum_ab = 0;
	sum_aa = 0;
	sum_bb = 0;
	for (i = 0; i < (n - 1); i++)
	{
		for (j = i; j < n; j++)
		{
			case2_pair(&i, &j, &n, dat_x, &mx, cdf_x, &a_mat[i + j * n]);
			a_mat[j + i * n] = a_mat[i + j * n];

			case2_pair(&i, &j, &n, dat_y, &my, cdf_y, &b_mat[i + j * n]);
			b_mat[j + i * n] = b_mat[i + j * n];
			sum_ab += a_mat[i + j * n] * b_mat[i + j * n];
			sum_aa += a_mat[i + j * n] * a_mat[i + j * n];
			sum_bb += b_mat[i + j * n] * b_mat[i + j * n];
		}

	}

	*tau = sum_ab / sqrt(sum_aa*sum_bb);

	free(a_mat);
	free(b_mat);
	return;
}


void is_possible_rank(int *N, double *rank, double *A, int *flag)
{
	int i, n = *N, r = 0;
	double sum = 0, crit = n * (n + 1)*0.5, sum_rank = 0;

	for (i = 0; i < n; i++)
	{
		r = (int)rank[i] -1;
		if (A[i + r * n] == 1.0)
			sum += 1.0;

		sum_rank += rank[i];
	}

	if ((sum == (double)n) & (sum_rank == crit))
	{
		*flag = 1;
	}
	else
		*flag = 0;
}




void init_rank(int *N, double *A, double *sample, int *flag)
{
	int i, j, k, n = *N;
	int count = 0, oth_r;
	int *pos_index, *initial_order;
	double *Atemp;
	double *initial_rank, *row_sum, cur_deg = 0;

	size_t sn = n * sizeof(double);
	size_t sn2 = n * n * sizeof(double);

	Atemp = (double *)malloc(sn2);
	initial_rank = (double *)malloc(sn);
	initial_order = (int *)malloc(sn);
	row_sum = (double *)malloc(sn);
	pos_index = (int *)malloc(sn);

	memcpy(Atemp, A, sn2);

	memset(initial_rank, 0, sn);
	memset(row_sum, 0, sn);
	memset(pos_index, 0, sn);

	for (j = 0; j < n; j++)
	{
		// Rank 1 - n (checking possibility)
		count = 0;
		for (i = 0; i < n; i++)
		{
			if (Atemp[i + j * n] == 1.0)
			{
				pos_index[count] = i;
				row_sum[i] = 0;
				for (k = j + 1; k < n; k++)
				{
					if (Atemp[i + k * n] == 1.0)
						row_sum[i] += 1.0;
				}
				count += 1;
			}
		}

		if (count == 0)
		{
			Rprintf("\n There is a problem in the given adjacency matrix.\n");
			*flag = 0;
			return;
		}

		initial_order[j] = pos_index[0];
		cur_deg = row_sum[pos_index[0]];
		for (i = 1; i < count; i++)
		{
			oth_r = pos_index[i];
			if (cur_deg >= row_sum[oth_r])
			{
				initial_order[j] = oth_r;
				cur_deg = row_sum[oth_r];
			}
		}

		//Rprintf("\n initial order[%d] = %d ", j, initial_order[j]);
		for (k = j; k < n; k++)
		{
			Atemp[initial_order[j] + n * k] = 0;
		}

	}

	for (j = 0; j < n; j++)
	{
		initial_rank[initial_order[j]] = (double) (j + 1);
	}

	memcpy(sample, initial_rank, sn);

	*flag = 1;

	free(initial_rank);
	free(Atemp);
	free(initial_order);
	free(row_sum);
	free(pos_index);

	return;
}

double dmod(double x, double y)
{
	return(x - ((int)(x / y)) * y);
}

void sampler(int *N, double *A, int *B, int *burnin, int *thining, double *sample)
{
	int i,n = *N, b = *B, th = *thining, br = *burnin, BB = b*th + br;
	int swap, c, d, flag, ival, irem;
	int *base_index, *temp_index;
	double *rank_sample, *rank_update, *rank_temp, temp=0;
	double u = 0;

	size_t sd_nb = n * b * sizeof(double);
	size_t si_n = n * sizeof(int);
	size_t sd_n = n * sizeof(double);


	base_index = (int *)malloc(si_n);
	temp_index = (int *)malloc(si_n);

	rank_sample = (double *)malloc(sd_nb);
	rank_update = (double *)malloc(sd_n);
	rank_temp = (double *)malloc(sd_n);

	init_rank(&n, A, rank_update, &flag);
	if (flag == 0)
	{
		Rprintf("\n There is a problem in the given adjacency matrix.\n");
		return;
	}

	
	for (i = 0; i < n; i++)
		base_index[i] = i;

	for (i = 0; i < BB; i++)
	{
		memcpy(temp_index, base_index, si_n);
		memcpy(rank_temp, rank_update, sd_n);

		u = runif(0.0, 1.0);
		c = (int)(n * u);
		
		swap = temp_index[n - 1];
		temp_index[n - 1] = temp_index[c];
		temp_index[c] = swap;
		u = runif(0.0, 1.0);
		d = temp_index[(int)((n - 1)*u)];

		//Rprintf("%.4f \t c = %d, d = %d \n", u, c, d);


		temp = rank_temp[c];
		rank_temp[c] = rank_temp[d];
		rank_temp[d] = temp;

		is_possible_rank(&n, rank_temp, A, &flag);
		if (flag)
		{
			memcpy(rank_update, rank_temp, sd_n);
		}

		//memcpy(&sample[0], rank_temp, sd_n);
		//return;


		ival = (int)(i - br) / th;
		irem = (int)(i - br) % th;

		if ((ival >= 0) & (irem == 0))
		{
			memcpy(&rank_sample[ival*n], rank_update, sd_n);
			if( dmod( (double) (ival+1)*10, (double) b) == 0.0)
			{
				Rprintf("%d%% ", (ival+1) * 100 / b);
			}
		}
	}
	Rprintf("\n");
	memcpy(sample, rank_sample, sd_nb);
	
	free(base_index);
	free(temp_index);
	free(rank_sample);
	free(rank_update);
	free(rank_temp);

	return;
}


// sign(X_l - X_k) 
double rank_sign(double Xk, double Xl)
{
	double temp = 0.0;
	if (Xk > Xl)
		temp = 1.0;
	if (Xk < Xl)
		temp = -1.0;
	return(temp);
}


void k_tau(int *N, double *X, double *Y, double *out)
{
	int i, j, n = *N;
	//int n2 = n*(n-1)/2;
	double temp = 0;
	double tau;
	double a2, b2, L, M, C;

	a2 = 0;
	b2 = 0;
	for (i = 0; i<(n - 1); i++)
	{
		for (j = i + 1; j<n; j++)
		{
			L = rank_sign(X[i], X[j]);
			M = rank_sign(Y[i], Y[j]);
			temp = temp + L * M;
			a2 = a2 + L * L;
			b2 = b2 + M * M;
		}
	}

	// tau
	C = sqrt(a2*b2);
	if (C == 0)
		tau = 0;
	else
		tau = temp / C;

	*out = tau;

}



void k_tau_tie(int *N, double *X, double *Y, double *Tx, double *Ty, double *out)
{
	int i, j, n = *N;
	//int n2 = n*(n-1)/2;
	double temp = 0;
	double tau;
	double a2, b2, L, M, C;

	a2 = 0;
	b2 = 0;
	for (i = 0; i<(n - 1); i++)
	{
		for (j = i + 1; j<n; j++)
		{
			L = rank_sign(X[i], X[j])*Tx[i+j*n];
			M = rank_sign(Y[i], Y[j])*Ty[i+j*n];
			temp = temp + L * M;
			a2 = a2 + L * L;
			b2 = b2 + M * M;
		}
	}

	// tau
	C = sqrt(a2*b2);
	if (C == 0)
		tau = 0;
	else
		tau = temp / C;

	*out = tau;

}

void calc_tau_vec(int *N, int *B, double *SX, double *SY, double *tau_vec)
{
	int i = 0, n = *N, b = *B;

	for (i = 0; i < b; i++)
	{
		k_tau(&n, &SX[i*n], &SY[i*n], &tau_vec[i]);
	}
	return;
}

void calc_tau_vec_tie(int *N, int *B, double *SX, double *SY, double *Tx, double *Ty,double *tau_vec)
{
	int i = 0, n = *N, b = *B;

	for (i = 0; i < b; i++)
	{
		k_tau_tie(&n, &SX[i*n], &SY[i*n], Tx, Ty, &tau_vec[i]);
	}
	return;
}
