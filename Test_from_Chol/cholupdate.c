//------------------------------------------------------------------------------
// SPEX_Update/Test/cholupdate.c: test Cholesky rank-1 update functionality
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/*
 * performance test
 */

#define FREE_WORKSPACE                           \
{                                                \
    glp_delete_prob(LP);                         \
    glp_free_env();                              \
    if (result_file != NULL) fclose(result_file);\
    SPEX_matrix_free(&Prob_A, option);           \
    SPEX_matrix_free(&Prob_b, option);           \
    SPEX_matrix_free(&Prob_c, option);           \
    SPEX_matrix_free(&A1, option);               \
    SPEX_matrix_free(&A2, option);               \
    SPEX_matrix_free(&A3, option);               \
    SPEX_matrix_free(&Atmp, option);             \
    SPEX_matrix_free(&L1, option);               \
    SPEX_matrix_free(&L2, option);               \
    SPEX_matrix_free(&L3, option);               \
    SPEX_matrix_free(&Ltmp, option);             \
    SPEX_matrix_free(&rhos1, option);            \
    SPEX_matrix_free(&rhos2, option);            \
    SPEX_matrix_free(&rhos3, option);            \
    SPEX_matrix_free(&rhostmp, option);          \
    SPEX_matrix_free(&rhos_update, option);      \
    SPEX_matrix_free(&L_update, option);         \
    SPEX_matrix_free(&A0, option);               \
    SPEX_vector_free(&w, option);                \
    SPEX_Chol_analysis_free(&analysis);          \
    SPEX_Chol_analysis_free(&best_analysis);     \
    SPEX_FREE(option);                           \
    mpz_clear(tmpz);                             \
    SPEX_FREE(P_update_inv);                     \
    SPEX_FREE(P_update);                         \
    SPEX_FREE(basis);                            \
    SPEX_FREE(used);                             \
    SPEX_finalize() ;                            \
}
#define PRINT_TO_FILE

#include "test.h"
#include "SPEX_Chol.h"
#include <assert.h>

int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    //char *prob_name = "lp_25fv47";
    char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
    //char *prob_name = "aa5";
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    SPEX_info info;
    //------------------------------------------------------------------ 
    // Initialize for SPEX library 
    //------------------------------------------------------------------ 
 
    SPEX_initialize () ; 
 
    //------------------------------------------------------------------ 
    // Initialize and allocate workspace 
    //------------------------------------------------------------------ 
    int64_t n, i, j, p, n_col, target = -1, sigma; 
    int sgn; 
    SPEX_options *option = NULL; 
    SPEX_Chol_analysis *analysis = NULL, *best_analysis = NULL; 
    SPEX_matrix *Prob_A = NULL, *Prob_b = NULL, *Prob_c = NULL; 
    SPEX_matrix *A1 = NULL, *L1 = NULL, *rhos1 = NULL, 
                *A2 = NULL, *L2 = NULL, *rhos2 = NULL, 
                *A3 = NULL, *L3 = NULL, *rhos3 = NULL, 
                *Atmp = NULL, *Ltmp = NULL, *rhostmp = NULL; 
    SPEX_matrix *rhos_update = NULL; 
    SPEX_matrix *L_update = NULL, *A0 = NULL; 
    SPEX_vector *w = NULL; 
    mpz_t tmpz; 
    int64_t *P_update_inv = NULL, *P_update = NULL; 
    int64_t *basis = NULL, *used= NULL; 
    clock_t start, end, start1, start2, end1, end2; 
    char file_name[1000]; 
    FILE *result_file = NULL; 
    glp_prob *LP; 
    glp_smcp parm; 
    LP = glp_create_prob(); 
    glp_init_smcp(&parm); 
    parm.it_lim = 1; 
    OK(SPEX_create_default_options(&option)); 
    OK(SPEX_mpz_init(tmpz)); 
 
    //-------------------------------------------------------------------------- 
    // open output file 
    //-------------------------------------------------------------------------- 
#ifdef PRINT_TO_FILE 
    sprintf(file_name, "Results/LPnetlib_CholUpdate/%s_chol.txt",prob_name); 
    result_file = fopen(file_name, "w"); 
    if (result_file == NULL) 
    { 
        printf("Could not open file: %s!\n", file_name); 
        FREE_WORKSPACE; 
        return 0; 
    } 
#endif 
 
    //-------------------------------------------------------------------------- 
    // read in matrix 
    //-------------------------------------------------------------------------- 
    sprintf(file_name, "TestMats/LPnetlib/%s/%s",prob_name,prob_name); 
    // read matrix A, arrays b, lb, ub and c to LP 
    double z0 = 0; 
    SPEX_construct_LP(LP, &Prob_A, &Prob_b, &Prob_c, &z0, file_name, option); 
    n = Prob_A->m; 
    n_col = Prob_A->n; // number of columns in matrix
     
    //-------------------------------------------------------------------------- 
    // generate initial basis matrix 
    //-------------------------------------------------------------------------- 
    basis = (int64_t*) SPEX_malloc(n_col * sizeof(int64_t)); 
    used  = (int64_t*) SPEX_calloc(n_col,  sizeof(int64_t)); 
    if (!basis || !used) 
    { 
        FREE_WORKSPACE; 
        return 0; 
    } 
 
    printf("getting initial basic variables....\n"); 
    glp_adv_basis(LP, 0); 
    glp_factorize(LP); // in order to use glp_get_bhead
    parm.msg_lev = GLP_MSG_OFF;
    while (glp_get_status(LP) != GLP_FEAS)
    {
        glp_simplex(LP, &parm);
    }
    if (glp_get_status(LP) == GLP_OPT)
    {
        printf("LP is optimized!\n");
        FREE_WORKSPACE;
        return 0;
    }
    for (i = 0; i < n; i++)
    {
        basis[i] = glp_get_bhead(LP, i+1)-1;
#if 0
        if (basis[i] < n && glp_get_mat_row(LP, basis[i]+1, NULL, NULL) != 0)
        {
            //printf("basis[%ld]=%ld\n",i,basis[i]);
            GOTCHA;
            return 0;;
        }
#endif
    }

    OK(SPEX_mpq_get_den(tmpz, Prob_A->scale));
    OK(SPEX_mpz_cmp_ui(&sgn, tmpz, 1));
    if (sgn != 0)
    {
        OK(SPEX_gmp_printf("scale is %Qd, whose den is not 1\n",
            Prob_A->scale));
        FREE_WORKSPACE;
        return 0;
    }

    // allocate A0 with n sparse vectors with initially n nnz
    OK(SPEX_matrix_allocate(&A0, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, n, 0, false,
        true, option));
    for (i = 0; i < n; i++)
    {
        OK(SPEX_vector_realloc(A0->v[i], n, option));
    }


    //--------------------------------------------------------------------------
    // compute A1 as B*B^T 
    //--------------------------------------------------------------------------
    printf("set of basic variables found, now computing A=B*B^T....\n");
    start = clock();
    for (i = 0; i < n; i++)
    {
        if (basis[i] < n)
        {
            j = basis[i];
            target = -1;
            for (p = 0; p < A0->v[j]->nz; p++)
            {
                if (A0->v[j]->i[p] == j) {target = p; break;}
            }
            if (target == -1)
            {
                target = A0->v[j]->nz;
                A0->v[j]->i[target] = j;
                A0->v[j]->nz++;
            }
            mpq_get_num(tmpz, Prob_A->scale);
            OK(SPEX_mpz_addmul(A0->v[j]->x[target], tmpz, tmpz));
        }
        else
        {
            j = basis[i]-n;
            used[j] = 1;

            OK(SPEX_A_plus_vvT(A0, Prob_A, j));
        }
    }
    OK(SPEX_mpq_mul(A0->scale, Prob_A->scale, Prob_A->scale));
    OK(SPEX_matrix_copy(&A1, SPEX_CSC, SPEX_MPZ, A0, option));

    //--------------------------------------------------------------------------
    // find the most sparse and dense vectors from remaining cols
    //--------------------------------------------------------------------------
    int64_t sparsest = -1, sparse_nz = n, densest = -1, dense_nz = 0;
    for (i = 0; i < n_col; i++)
    {
        if (used[i] != 1)
        {
            int64_t Ai_nz = Prob_A->p[i+1] - Prob_A->p[i];
            if (Ai_nz < sparse_nz)
            {
                sparsest = i;
                sparse_nz = Ai_nz;
            }
            if (Ai_nz > dense_nz)
            {
                densest = i;
                dense_nz = Ai_nz;
            }
        }

    }
 
    if (sparsest == -1 || densest == -1) 
    { 
        printf("cannot find any unused column\n"); 
        return 0; 
    } 
 
    //-------------------------------------------------------------------------- 
    // compute A2 as A1+v_sparsest*v_sparsest^T 
    //-------------------------------------------------------------------------- 
    OK(SPEX_A_plus_vvT(A0, Prob_A, sparsest)); 
    OK(SPEX_matrix_copy(&A2, SPEX_CSC, SPEX_MPZ, A0, option)); 
 
    //-------------------------------------------------------------------------- 
    // compute A3 as A2+v_densest*v_densest^T 
    //--------------------------------------------------------------------------
    OK(SPEX_A_plus_vvT(A0, Prob_A, densest));
    OK(SPEX_matrix_copy(&A3, SPEX_CSC, SPEX_MPZ, A0, option));
    end = clock();
    printf("time to compute BB^T: %lf\n", (double) (end - start)/CLOCKS_PER_SEC);

    //--------------------------------------------------------------------------
    // perform Cholesky factorization for A1, A2, A3 with same pre-ordering
    //--------------------------------------------------------------------------
    printf("compute initial Cholesky factorization for A....\n");
    bool left_looking = true;// True = left, false = up
    OK(SPEX_Chol_preorder(&analysis, A1, option));
    OK(SPEX_Chol_permute_A(&Atmp, A1, analysis));
    OK(SPEX_Chol_Factor(&L1, &rhos1, analysis, Atmp, left_looking, option));

    SPEX_matrix_free(&Atmp, NULL);
    OK(SPEX_Chol_permute_A(&Atmp, A2, analysis));
    OK(SPEX_Chol_Factor(&L2, &rhos2, analysis, Atmp, left_looking, option));

    SPEX_matrix_free(&Atmp, NULL);
    OK(SPEX_Chol_permute_A(&Atmp, A3, analysis));
    OK(SPEX_Chol_Factor(&L3, &rhos3, analysis, Atmp, left_looking, option));

    //--------------------------------------------------------------------------
    // generate initial inputs for LU update
    //--------------------------------------------------------------------------
    printf("generating inputs for Cholesky rank-1 update....\n");
    // generate permutation vectors P, Q, P_inv, Q_inv and vectors sd
    P_update     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P_update_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    if (!P_update || !P_update_inv)
    {
        FREE_WORKSPACE;
        return 0;
    }
    OK(SPEX_matrix_allocate(&rhos_update, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));
    for (i = 0; i < n; i++)
    {
        P_update[i] = analysis->q[i];
        P_update_inv[P_update[i]] = i;
        OK(SPEX_mpz_set(SPEX_1D(rhos_update, i, mpz), SPEX_1D(rhos1, i, mpz)));
    }

    // convert the factorization to SPEX_mat to be used in the update process
    OK(SPEX_matrix_copy(&L_update, SPEX_DYNAMIC_CSC, SPEX_MPZ, L1, option));
    OK(SPEX_Update_permute_row(L_update, P_update, option));
    OK(SPEX_Update_matrix_canonicalize(L_update, P_update, option));

    // allocate space for scattered vector w
    OK(SPEX_vector_allocate(&w, 0, option));
    OK(SPEX_vector_realloc(w, n, option));

    for (int iter = 0; iter < 4; iter++)
    {
        if (iter < 2) sigma = 1; // the first two iterations are updates
        else sigma = -1;// while the last two iterations are downdates

        //----------------------------------------------------------------------
        // get the vector w
        //----------------------------------------------------------------------
        int64_t new_col;
        if (iter == 1 || iter == 2) new_col = densest;
        else new_col = sparsest;

        j = 0;
        for (p = Prob_A->p[new_col]; p < Prob_A->p[new_col+1]; p++)
        {
            OK(SPEX_mpz_set(w->x[j], SPEX_1D(Prob_A, p, mpz)));
            w->i[j] = Prob_A->i[p];
            j++;
        }
        w->nz = j;
        // used[new_col] = 1;

        //----------------------------------------------------------------------
        // perform update
        //----------------------------------------------------------------------
        printf("computing Cholesky rank-1 update if col %ld is %s B...\n",
            new_col, iter < 2?"added to":"removed from");
        start1 = clock();
        OK(SPEX_Update_Chol_Rank1(L_update, rhos_update, P_update, P_update_inv,
            w, sigma, option));
        end1 = clock();

        //----------------------------------------------------------------------
        // compute updated Cholesky using direct factorization
        //----------------------------------------------------------------------
        printf("computing updated Cholesky using direct factorization...\n");
        SPEX_Chol_analysis_free(&best_analysis);
        SPEX_matrix_free(&Atmp, NULL);
        SPEX_matrix_free(&rhostmp, NULL);
        SPEX_matrix_free(&Ltmp, NULL);

        // perform Cholesky factorization
        start2 = clock();
        if (iter == 0)
        {
            OK(SPEX_Chol_preorder(&best_analysis, A2, option));
            OK(SPEX_Chol_permute_A(&Atmp, A2, best_analysis));
        }
        else if (iter == 1)
        {
            OK(SPEX_Chol_preorder(&best_analysis, A3, option));
            OK(SPEX_Chol_permute_A(&Atmp, A3, best_analysis));
        }
        else if(iter == 3)
        {
            OK(SPEX_Chol_preorder(&best_analysis, A1, option));
            OK(SPEX_Chol_permute_A(&Atmp, A1, best_analysis));
        }
        if (iter != 2)// this will be the same as when iter==0
        {
            OK(SPEX_Chol_Factor(&Ltmp, &rhostmp, best_analysis, Atmp,
                left_looking, option));
        }
        end2 = clock();

        bool Isequal;
        size_t s, L_sum_size = 0, L_update_sum_size = 0;
        int64_t L_nnz, L_update_nnz = 0;
        if (iter == 0 || iter == 2)
        {
            SPEX_matrix_equal(&Isequal, L2, L_update, P_update);
            for (p = 0; p < L2->p[n]; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, L2->x.mpz[p], 2));
                L_sum_size += s;
            }
            L_nnz = L2->p[n];
        }
        else if (iter == 1)
        {
            SPEX_matrix_equal(&Isequal, L3, L_update, P_update);
            for (p = 0; p < L3->p[n]; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, L3->x.mpz[p], 2));
                L_sum_size += s;
            }
            L_nnz = L3->p[n];
        }
        else // if(iter == 3)
        {
            SPEX_matrix_equal(&Isequal, L1, L_update, P_update);
            for (p = 0; p < L1->p[n]; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, L1->x.mpz[p], 2));
                L_sum_size += s;
            }
            L_nnz = L1->p[n];
        }
        if (!Isequal) return 0;
        for (i = 0; i < n; i++)
        {
            L_update_nnz += L_update->v[i]->nz;
            for (p = 0; p < L_update->v[i]->nz; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, L_update->v[i]->x[p], 2));
                L_update_sum_size += s;
            }
        }
        //----------------------------------------------------------------------
        // print results
        //----------------------------------------------------------------------
        // Timing stats
        double t1 = (double) (end1 - start1)/CLOCKS_PER_SEC;
        double t2 = (double) (end2 - start2) / CLOCKS_PER_SEC;

        printf("test succeeded!!");
        printf("\n\t\t time \t\tL_nnz \tsum(Lb)");
        printf("\nFactorization: \t%lf \t%ld \t%ld", t2,L_nnz,L_sum_size);
        printf("\nrank-1 Update: \t%lf \t%ld \t%ld\n\n", t1,L_update_nnz,
            L_update_sum_size);
#ifdef PRINT_TO_FILE
        if (iter == 2)
        {
            fprintf(result_file,"\t%lf \t%ld \t%ld",t1,L_update_nnz,
                L_update_sum_size);
        }
        else
        {
            fprintf(result_file,"\t%lf \t%ld \t%ld",t2,L_nnz,L_sum_size);
            fprintf(result_file,"\t%lf \t%ld \t%ld",t1,L_update_nnz,
                L_update_sum_size);
        }
#endif
    }

    FREE_WORKSPACE;
    return 0;
}

