#include <math.h>
// main.h has both type definitions, including simple.h from typedef.h, and network.h which includes MPI communication
#include "main.h"
// gmx wrappers for omp stuff
#include "../utility/gmxomp.h"

#define UMB_MAX_OMP 256

typedef struct Umbrella_Data {
  real k_Q;
  real Q_0;
  real Q_init;
  real Q_steps;
  gmx_int64_t step;
  t_commrec *cr;
  int n_omp;
  int freq;
  FILE *fp;
} Umbrella_Data;

static Umbrella_Data udata;

void Init_Umbrella_Communicate(t_commrec *cr)
{
  int i;
  // int i_mpi;
  char fnm[256];
  FILE *fpin;

  fpin=fopen("umbrella_params","r");

  if (fpin!=NULL) {
    fprintf(stderr,"Found umbrella parameter file\n");
    // MPI_Comm_rank(MPI_COMM_WORLD,&i_mpi);
    i=fscanf(fpin,"freq_out %d",&udata.freq);
    // if (i==0 || i_mpi!=0) {
    if (i==0 || gmx_node_rank()!=0) {
      udata.freq=-1;
      udata.fp=NULL;
    } else {
      i=0;
      udata.fp=NULL;
      do {
        if (udata.fp) {
          fclose(udata.fp);
        }
        i+=1;
        sprintf(fnm,"Qumbrella.part%04d.dat",i);
      } while (udata.fp=fopen(fnm,"r"));
      fprintf(stderr,"Writing umbrella coordinates to %s\n",fnm);
      udata.fp=fopen(fnm,"w");
      if (udata.fp==NULL) {
        fprintf(stderr,"Warning, outfile %s open failed\n",fnm);
      }
    }

    udata.k_Q=0;
    udata.Q_0=0;
    udata.Q_init=0;
    udata.Q_steps=1;
    udata.cr=cr;
    i=fscanf(fpin,"%f %f %f %f",&udata.k_Q,&udata.Q_0,&udata.Q_init,&udata.Q_steps);
    if (i<4) {
      fprintf(stderr,"Warning, umbrella parameters read incorrectly\n");
    }

    fclose(fpin);
  } else {
    udata.k_Q=0;
    udata.Q_0=0;
    udata.Q_init=0;
    udata.Q_steps=1;
    udata.cr=cr;
    udata.freq=-1;
    udata.fp=NULL;
  }
}

void Umbrella_Set_Step(int n_omp,gmx_int64_t step)
{
  udata.n_omp=n_omp;
  udata.step=step;
}

real Umbrella_Communicate(real Q_local,real *k_Q,real *Q_0)
{
  int i;
  int i_omp=gmx_omp_get_thread_num();
  // int n_omp=gmx_omp_get_num_procs(); // Wrong number
  int n_omp=udata.n_omp;
  // int i_mpi;
  static real Q_semilocal[UMB_MAX_OMP];
  real Q_global;

  if (n_omp>UMB_MAX_OMP) {
    fprintf(stderr,"Seg fault is probably about to happen because Q_semilocal is not big enough to accommodate %d omp threads. See %d in %s.\n",n_omp,__LINE__,__FILE__);
  }

  Q_semilocal[i_omp]=Q_local;
  #pragma omp barrier
  // #pragma omp master
  // {
  if (i_omp==0) {
    Q_local=0;
    // fprintf(stderr,"%f %f %f %f %f\n",Q_semilocal[0],Q_semilocal[1],Q_semilocal[2],Q_semilocal[3],Q_semilocal[4]);
    for (i=0; i<n_omp; i++) {
      Q_local+=Q_semilocal[i];
    }
    // fprintf(stderr,"%f %f %f %f %f\n",Q_semilocal[0],Q_semilocal[1],Q_semilocal[2],Q_semilocal[3],Q_semilocal[4]);
//    #ifdef GMX_DOUBLE
//    MPI_Allreduce(&Q_local,&Q_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    #else
//    MPI_Allreduce(&Q_local,&Q_global,1,MPI_FLOAT, MPI_SUM,MPI_COMM_WORLD);
//    #endif
//    for (i=0; i<n_omp; i++) {
//      Q_semilocal[i]=Q_global;
//    }
    // gmx_sum declared in src/gromacs/legacyheaders/network.h (in main.h)
    #ifdef GMX_MPI
    gmx_sum(1,&Q_local,udata.cr);
    #endif
    for (i=0; i<n_omp; i++) {
      Q_semilocal[i]=Q_local;
    }
    // fprintf(stderr,"%f %f %f %f %f\n",Q_semilocal[0],Q_semilocal[1],Q_semilocal[2],Q_semilocal[3],Q_semilocal[4]);
  }
  #pragma omp barrier
  Q_global=Q_semilocal[i_omp];
  *k_Q=udata.k_Q;
  *Q_0=udata.Q_0+(udata.Q_init-udata.Q_0)*exp(-udata.step/udata.Q_steps);

  // #pragma omp master
  // {
  if (i_omp==0) {
    // MPI_Comm_rank(MPI_COMM_WORLD,&i_mpi);
    // if (i_mpi==0) { // udata.fp=NULL on other mpi processes.
      if ((udata.step % udata.freq)==0 && udata.fp!=NULL) {
        fprintf(udata.fp,"%d %g %g\n",udata.step,Q_global,0.5*(*k_Q)*(Q_global-(*Q_0))*(Q_global-(*Q_0)));
      }
    // }
  }

  return Q_global;
}

void Free_Umbrella_Communicate(void)
{
  // int i_mpi;

  // MPI_Comm_rank(MPI_COMM_WORLD,&i_mpi);
  // if (i_mpi==0) {
    if (udata.fp!=NULL) {
      fclose(udata.fp);
    }
  // }
}
