/*------------------------------------------------------------------------*/
/*  File: kQKP_inst.c from L. Letocart                                    */
/*------------------------------------------------------------------------*/


/***************************************************************************/
/* exact solution method for solving the k-item quadratic knapsack problem */
/***************************************************************************/

// note that we consider only integer weights & costs!

#include <string.h>
#include <stdlib.h>
#include <iostream>  
#include <stdio.h>


using namespace std;

// #include "generate_instance_kQKP.h"

#define srand(x)     srand48(x)
#define randm(x)    (lrand48() % (long) (x))


// to be used by qsort
typedef struct {
    int valeur;
    int indice;
} tabint_t;

static int compare_tabint_t(const void *a, const void *b)
{
    const tabint_t *ra = (const tabint_t *)a;
    const tabint_t *rb = (const tabint_t *)b;
    if (ra->valeur < rb->valeur)
        return(-1);
    else if (ra->valeur == rb->valeur)
        return(0);
    else
        return(1);
}


/* ======================================================================
			      generate_instance
   ====================================================================== */

int generate_instance(int n, int r, int pct, int v, int *p,int *w,int *c,int *nbobj)
{
  /* generate instance */
  int i, j;
  int wsum;
  FILE *inst;
  tabint_t *tabInt = (tabint_t *) malloc( n * sizeof(tabint_t) ); //used by qsort

  srand(v+n+r+pct);
  char buffer[100];
  sprintf(buffer, "./instances/QKP_%d_%d_%d", n, r, pct);
  inst = fopen(buffer, "w");
  fprintf(inst,"n=%d\n", n);
  fprintf(inst,"density=%d\n", pct);

  /* generate profits and weights */
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      p[i*n+j] = p[j*n+i] = (randm(100) >= pct ? 0 : randm(r)+1);
      // fprintf(inst,"%d\t", p[i*n+j]);
    }
    // fprintf(inst,"\n");

    w[i] = randm(r/2)+1;
  }


  fprintf(inst,"Q1=[");
  for (i = 0; i < n-1; i++) {
    for (j = 0; j < n; j++) {
      (j<=i) ? fprintf(inst,"%d ", p[i*n+j]) : fprintf(inst,"%d ", 0);
    }
    fprintf(inst,";\n");
  }
  i=n-1;
  for (j = 0; j < n; j++) {
      fprintf(inst,"%d ", p[i*n+j]);
  }
  fprintf(inst,"]\n");



  fprintf(inst,"w=[");
  for (i = 0; i < n; i++) {
    tabInt[i].valeur=w[i];
    tabInt[i].indice=i;
  }
  for (i = 0; i < n-1; i++) fprintf(inst,"%d, ", w[i]);
  fprintf(inst,"%d]\n", w[n-1]);

  wsum = 0;
  for (i = 0; i < n; i++) 
    wsum += w[i];
  if (wsum - 50 <= 0) {
    printf("too small weight sum\n");
    printf("PROGRAM TERMINATED!\n\n");
    exit(-1);
  }
  *nbobj = randm(n/4) + 1; /* generate number of items k */
  if (*nbobj==1) *nbobj=2;// if k==1 then k=2
  /* qsort on weights w */
  qsort(tabInt, n, sizeof(tabint_t), compare_tabint_t);
//  printf("After sort\n");
  wsum = 0;
  for (i = 0; i < *nbobj; i++) wsum += tabInt[i].valeur;  
//   c    = randm(wsum-50) + 50;
  *c    = randm(((*nbobj)*30)-50) + 50; /* generate capacity */
//decrease k until wsum>c
  while (wsum>*c) {
    wsum-=tabInt[*nbobj-1].valeur;
    *nbobj=*nbobj-1;
  }
  fprintf(inst,"W=%d \n", *c);
  fprintf(inst,"k=%d \n", *nbobj);

  fclose(inst);
  return 0;
}


/* main */

int main(int argc, char* argv[])
{
  int m_n,m_r,m_pct;
  int m_v;               // instance number
  int *m_C;             // matrix P
  int *m_a;              // weight vector a
  int m_b;               // capacity, i.e. rhs of knapsack constraint
  int m_k;               // number of objects to be selected


  // Check input
  if (argc == 5) {
    m_n = atoi(argv[1]); // nb objets
    m_r = atoi(argv[2]); // intervalle coeff ex 100
    m_pct = atoi(argv[3]); // density matrix cij, ex 25, 50, 100
    m_v = atoi(argv[4]);  // alea
  } else {
    cerr<< "Usage: " << argv[0] << " n r d v" <<endl;
    return 1;
  }

  // generate instance, i.e., cost matrix C, weights a, capacity b, number of items k

  m_C = new int[m_n*m_n];
  m_a = new int[m_n];
  for(int i = 0; i < m_n*m_n; i++) m_C[i] = 0;
  for(int i = 0; i < m_n; i++) m_a[i] = 0;

  generate_instance(m_n,m_r,m_pct,m_v,m_C,m_a,&m_b,&m_k);


  delete[] m_C;
  delete[] m_a;
  return 0;
}
