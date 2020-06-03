#ifndef LDLTSOLVER_H
#define LDLTSOLVER_H

#include <cholmod.h>

class LDLTTrip
{
private:
	cholmod_triplet *m_trip;
	cholmod_sparse *m_sparse;
	int *prow, *pcol;
	double *pval;
	cholmod_common com;

public:
	LDLTTrip(size_t nrow, size_t ncol, size_t nmax)
	{
		cholmod_start(&com);
		com.itype = CHOLMOD_INT;
		com.dtype = CHOLMOD_DOUBLE;
		m_sparse = NULL;
		m_trip = cholmod_allocate_triplet(nrow, ncol, nmax, 1, CHOLMOD_REAL, &com);
		prow = (int *)m_trip->i;
		pcol = (int *)m_trip->j;
		pval = (double *)m_trip->x;
		m_trip->nnz = 0;
	}

	~LDLTTrip()
	{
		cholmod_free_triplet(&m_trip, &com);
		if (m_sparse)
			cholmod_free_sparse(&m_sparse, &com);
		cholmod_finish(&com);
	}

	void clear()
	{
		m_trip->nnz = 0;
	}

	void addTrip(int i, int j, double v)
	{
		prow[m_trip->nnz] = i;
		pcol[m_trip->nnz] = j;
		pval[m_trip->nnz] = v;
		m_trip->nnz++;
	}

	cholmod_sparse * build()
	{
		if(m_sparse)
			cholmod_free_sparse(&m_sparse, &com);
		m_sparse = cholmod_triplet_to_sparse(m_trip, m_trip->nnz, &com);
		FILE *fp = fopen("sparse.txt", "w");
		cholmod_write_sparse(fp, m_sparse, NULL, "sparse.txt", &com);
		fclose(fp);
		return m_sparse;
	}
};

class LDLTSolver
{
private:
	cholmod_common m_cholmod;
	cholmod_factor* m_cholmodFactor;
	cholmod_dense *m_x, *m_Y, *m_E;
	double m_shiftOffset[2];
	
public:
	LDLTSolver()
	{
		m_x = m_Y = m_E = NULL;
		m_shiftOffset[0] = m_shiftOffset[1] = 0;
		cholmod_start(&m_cholmod);
		init();
	}

	void compute(cholmod_sparse *A)
	{
		analyzePattern(A);
		factorize(A);
	}

	double * solve(double *b, int m, int n = 1)
	{
		cholmod_dense b_cd;
		b_cd.nrow = m;
		b_cd.ncol = n;
		b_cd.nzmax = m * n;
		b_cd.d = m;
		b_cd.x = (void *)b;
		b_cd.xtype = CHOLMOD_REAL;
		b_cd.dtype = CHOLMOD_DOUBLE;

		cholmod_solve2(CHOLMOD_LDLt, m_cholmodFactor, &b_cd, NULL, &m_x, NULL, &m_Y, &m_E, &m_cholmod);
		return (double *)(m_x->x);
	}

	~LDLTSolver()
	{
		if (m_cholmodFactor)
			cholmod_free_factor(&m_cholmodFactor, &m_cholmod);
		if (m_x) cholmod_free_dense(&m_x, &m_cholmod);
		if (m_Y) cholmod_free_dense(&m_Y, &m_cholmod);
		if (m_E) cholmod_free_dense(&m_E, &m_cholmod);
		cholmod_finish(&m_cholmod);
	}

private:
	void init()
	{
		m_cholmod.final_asis = 1;
		m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
	}

	void analyzePattern(cholmod_sparse *A)
	{
		if (m_cholmodFactor)
		{
			cholmod_free_factor(&m_cholmodFactor, &m_cholmod);
			m_cholmodFactor = 0;
		}
		m_cholmodFactor = cholmod_analyze(A, &m_cholmod);
	}

	void factorize(cholmod_sparse *A)
	{
		cholmod_factorize_p(A, m_shiftOffset, 0, 0, m_cholmodFactor, &m_cholmod);
	}
};

#endif // !LDLTSOLVER_H
