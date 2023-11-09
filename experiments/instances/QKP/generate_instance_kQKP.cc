
/*------------------------------------------------------------------------*/
/*  File: generate_instance_kQKP.cc from L. Letocart                      */
/*------------------------------------------------------------------------*/


/* This code generates instances for the exact k-item quadratic knapsack problem */

extern "C" {

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

 
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
  inst = fopen("Instance.txt", "w");
  fprintf(inst,"#E-kQKP Instance:%d\n\n%d #n: number of items\n", v, n);

  /* generate profits and weights */
  fprintf(inst,"\n#Profits:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      p[i*n+j] = p[j*n+i] = (randm(100) >= pct ? 0 : randm(r)+1);
      fprintf(inst,"%d\t", p[i*n+j]);
    }
    fprintf(inst,"\n");
//     printf("i=%d pii=%d  ", i, p[i][i]);
    w[i] = randm(r/2)+1;
  }
  fprintf(inst,"\n#Weights:\n");
//   printf("Before sort\n");
  for (i = 0; i < n; i++) {
    tabInt[i].valeur=w[i];
    tabInt[i].indice=i;
  }
  fprintf(inst,"\n");
  for (i = 0; i < n; i++) fprintf(inst,"%d\t", w[i]);
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
  fprintf(inst,"\n%d #capacity\n", *c);
  fprintf(inst,"\n%d #k: cardinality\n\n", *nbobj);

  fclose(inst);
  return 0;
}


}
