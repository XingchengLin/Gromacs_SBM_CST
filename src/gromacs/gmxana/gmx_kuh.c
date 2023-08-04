/*
 * $Id: gmx_bond.c,v 1.3.2.2 2007/09/19 10:54:49 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 *
 * ifdef HAVE_CONFIG_H
 *  include <config.h>
 *  endif
 *  include <math.h>
 *  include <string.h>
 *  include "sysstuff.h"
 *  include "typedefs.h"
 *  include "smalloc.h"
 *  include "macros.h"
 *  include "vec.h"
 *  include "pbc.h"
 *  include "xvgr.h"
 *  include "copyrite.h"
 *  include "gmx_fatal.h"
 *  include "futil.h"
 *  include "statutil.h"
 *  include "index.h"
 *  include "gstat.h"
 *  include "confio.h"
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <ctype.h>

#include "sysstuff.h"
#include <string.h>
#include "gromacs/utility/cstringutil.h"
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/commandline/pargs.h"
#include "index.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "physics.h"
#include "gmx_ana.h"
#include "macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/3dview.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/fileio/futil.h"


#include "gromacs/legacyheaders/gmx_fatal.h"
enum eqifmt { eqiNONE, eqiINDEX, eqiPAIR, eqiLIST } ;

typedef struct {
   int nmax1, nmax2 ;
   int size1, size2 ;
   const atom_id *members1, *members2 ; // [size] numbers 0..N-1
   int *ids1, *ids2 ; // [N] numbers 1 .. size, NOTE 0 is not present
   const char *name1, *name2 ;
   gmx_bool onegrp;
   real *dist ;
} GrpDists ;

static real group_dist( const GrpDists *dists, int i, int j )
{
   return dists->dist[i*dists->size2+j] ;
}

static void set_group_dist( GrpDists *dists, int i, int j, real d )
{
   dists->dist[i*dists->size2+j] = d ;
   if (dists->onegrp)
      dists->dist[j*dists->size1+i] = d ;
}

static void init_group_dists( GrpDists *dists, int nmax1, int nmax2,
      int size1, const atom_id idx1[], const char *name1,
      int size2, const atom_id idx2[], const char *name2 )
{ // just copying redundant data if onegrp == TRUE
   int i ;
   dists->nmax1 = nmax1 ;
   dists->nmax2 = nmax2 ;
   dists->size1  = size1 ;
   dists->size2  = size2 ;
   dists->members1 = idx1 ;
   dists->members2 = idx2 ;
   snew(dists->ids1,nmax1);
   snew(dists->ids2,nmax2);
   dists->name1 = name1 ;
   dists->name2 = name2 ;
   for (i=0;i<size1;++i)
      dists->ids1[idx1[i]]=i+1;
   for (i=0;i<size2;++i)
      dists->ids2[idx2[i]]=i+1;
   snew(dists->dist,size1*size2) ;
   dists->onegrp = (strcmp(name1,name2)==0 && nmax1 == nmax2) ;
}

static void clone_group_dists( const GrpDists *tpl, GrpDists *dists )
{
   int i ;
   dists->nmax1 = tpl->nmax1 ;
   dists->nmax2 = tpl->nmax2 ;
   dists->size1  = tpl->size1 ;
   dists->size2  = tpl->size2 ;
   dists->members1 = tpl->members1 ;
   dists->members2 = tpl->members2 ;
   snew(dists->ids1,dists->nmax1);
   snew(dists->ids2,dists->nmax2);
   dists->name1 = tpl->name1 ;
   dists->name2 = tpl->name2 ;
   for (i=0;i<dists->nmax1;++i)
      dists->ids1[i]=tpl->ids1[i];
   for (i=0;i<dists->nmax2;++i)
      dists->ids2[i]=tpl->ids2[i];
   snew(dists->dist,dists->size1*dists->size2) ;
   dists->onegrp = tpl->onegrp ;
}

//static void free_group_dists( GrpDists *dists )
//  {
//  // FIXME
//  }

static void calc_all_dists( GrpDists *dists, rvec x[], const t_pbc *pbc )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const atom_id *idx1 = dists->members1 ;
   const atom_id *idx2 = dists->members2 ;
   real *mat = dists->dist ;

   int i,j ;
   rvec    dx;
   real d ;

   for (i=0;i<size1;++i)
      for (j=0;j<size2;++j) {
	 pbc_dx(pbc,x[i],x[j],dx);
	 d = norm(dx) ;
	 mat[i*size2+j] = d ;
      }
}

static void calc_all_one_dists( GrpDists *dists, rvec x[], const t_pbc *pbc )
{
   int size = dists->nmax1 ;
   real *mat = dists->dist ;
   int i,j ;
   rvec    dx;
   real d ;

   for (i=0;i<size;++i) {
      mat[i*size+i] = 0.0 ;
      for (j=i+1;j<size;++j) {
	 pbc_dx(pbc,x[i],x[j],dx);
	 d = norm(dx) ;
	 mat[i*size+j] = d ;
	 mat[j*size+i] = d ;
      }
   }
}

static void calc_one_group_dists( GrpDists *dists, rvec x[], const t_pbc *pbc )
{
   int size = dists->size1 ;
   const atom_id *idx = dists->members1 ;
   real *mat = dists->dist ;
   int i,j,ai,aj ;
   rvec    dx;
   real d ;  

   for (i=0;i<size;++i) {
      mat[i*size+i] = 0.0 ;
      ai = idx[i] ;
      for (j=i+1;j<size;++j) {
	 aj = idx[j] ;
	 pbc_dx(pbc,x[ai],x[aj],dx);
	 d = norm(dx) ;
	 mat[i*size+j] = d ;
	 mat[j*size+i] = d ;
      }
   }
}

static void calc_two_group_dists( GrpDists *dists, rvec x[], const t_pbc *pbc )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const atom_id *idx1 = dists->members1 ;
   const atom_id *idx2 = dists->members2 ;
   real *mat = dists->dist ;

   int i,j,ai,aj ;
   rvec    dx;
   real d ;

   for (i=0;i<size1;++i) {
      ai = idx1[i] ;
      for (j=0;j<size2;++j) {
	 aj = idx2[j] ;
	 pbc_dx(pbc,x[ai],x[aj],dx);
	 d = norm(dx) ;
	 mat[i*size2+j] = d ;
      }
   }
}

static void calc_group_dists( GrpDists *dists, rvec x[], int ePBC, matrix box )
{
   t_pbc   pbc;

   set_pbc(&pbc,ePBC,box);
   if (dists->onegrp)
   {
      if (dists->size1 == dists->nmax1 )
	 calc_all_one_dists( dists, x, &pbc ) ;
      else
	 calc_one_group_dists( dists, x, &pbc ) ;
   }
   else
   { // not optimized for only one group with size == max
      if (dists->size1 == dists->nmax1 && dists->size2 == dists->nmax2)
	 calc_all_dists( dists, x, &pbc ) ;
      else
	 calc_two_group_dists( dists, x, &pbc ) ;
   } // ! onegrp
}

static gmx_bool read_cont_dists( const char* fname, int *gnx,
      atom_id **index, real **dists )
{
   int ncont ;
   int idxlen ;
   atom_id *idx ;
   real *d ;
   int ai, aj ;
   float rij ;
   int i, ci ;

   FILE *fp = gmx_ffopen(fname,"r") ;
   if (!fp) {
      fprintf(stderr,"FAILED to open file") ;
      exit(1) ;
   }

   fscanf(fp,"%d",&ncont);
   idxlen = 2*ncont ;
   snew(idx,idxlen) ;
   snew(d,ncont) ;

   i=ci=0 ;
   while ( fscanf(fp,"%d %d %f",&ai,&aj,&rij) == 3 && ci < ncont )
   {
      idx[i++] = ai-1 ; // index uses 0..N-1
      idx[i++] = aj-1 ;
      d[ci++] = rij ;
   }

   if (ci != ncont)
   {
      sfree(idx) ;
      sfree(d) ;
      return FALSE ;
   }
   else
   {
      *gnx = idxlen ;
      *index = idx ;
      *dists = d ;
      return TRUE ;
   }
}


// from gmx_rmsdist.c
static void calc_cont_dists(int ncont, const atom_id index[], rvec x[],
      int ePBC,matrix box, real *d)
{
   int     i,ci;
   rvec    dx;
   t_pbc   pbc;

   set_pbc(&pbc,ePBC,box);
   i=0;
   for(ci=0; (ci<ncont); ++ci ) {
      pbc_dx(&pbc,x[index[i]],x[index[i+1]],dx);
      d[ci]=norm(dx);
      //    fprintf(stderr,"%6d%8d%8d%14.5f\n",ci+1,index[i]+1,index[i+1]+1,d[ci]);
      i+=2;
   }
}

static void get_contact_resnums( const t_atoms *atoms, const atom_id index[],
      int ncont, int *resnums )
{
   int i, res, ni=2*ncont ;
   for ( i=0; i<ni; ++i )
   {
      res = atoms->atom[index[i]].resind ;
      resnums[i] = res ;
   }
   // atoms.resinfo[resind].nr would be the residue number,
   // the index itself also serves
}

static void get_all_resnums( const t_atoms *atoms, int natoms,
      int *resnums, int* resatoms )
{
   int i, res ;
   for (i=0;i<natoms;++i) {
      res = atoms->atom[i].resind ;
      resnums[i] = res ;
      ++resatoms[res] ;
   }
}

static void get_group_resnums( const t_atoms *atoms, const atom_id index[],
      int size, int *resnums, int* resatoms )
{                           // resnums against group index
   int i,res;
   for (i=0;i<size;++i) {
      res = atoms->atom[index[i]].resind ;
      resnums[i] = res ;
      ++resatoms[res] ;  
   }
}

static void make_residue_index( const int resatoms[], int nres,
      int *resindex[], int *size )
{
   int i, n=0;
   for (i=0;i<nres;++i)
      if (resatoms[i]>0)
	 ++n ;
   snew(*resindex,n) ;
   n=0 ;
   for (i=0;i<nres;++i)
      if (resatoms[i]>0)
	 (*resindex)[n++]=i ;
   *size = n ;
}


static void avg_res_res_dists( const GrpDists *dists,
      const int *resnum1, const int *resnum2,
      GrpDists *resdists )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;

   int rsize1 = resdists->size1 ;
   int rsize2 = resdists->size2 ;
   const int *rids1 = resdists->ids1 ;
   const int *rids2 = resdists->ids2 ;

   real *rdist = resdists->dist ;

   int i,j,ridi,ridj ;
   real d ;
   int *count ;

   memset(rdist, 0, rsize1*rsize2*sizeof(real));
   snew(count,rsize1*rsize2);
   for (i=0;i<size1;++i) {
      ridi = rids1[resnum1[i]]-1 ;
      for (j=0;j<size2;++j) {
	 ridj = rids2[resnum2[j]]-1 ;
	 d = dist[i*size2+j] ;
	 rdist[ridi*rsize2+ridj] += d ;
	 ++count[ridi*rsize2+ridj] ;
      }
   }

   for (i=0;i<rsize1;++i)
      for (j=0;j<rsize2;++j)
	 if (count[i*rsize2+j]) {
	    rdist[i*rsize2+j] /= count[i*rsize2+j] ;
	 } 

   sfree(count) ;
}

static void avg_atom_res_dists( const GrpDists *dists, const int *resnum2,
      GrpDists *resdists )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;

   int rsize2 = resdists->size2 ;
   const int *rids2 = resdists->ids2 ;

   real *rdist = resdists->dist ;
   int i,j,ridj ;
   real d ;
   int *count ;

   memset(rdist, 0, size1*rsize2*sizeof(real));
   snew(count,size1*rsize2);
   for (i=0;i<size1;++i) {
      for (j=0;j<size2;++j) {
	 ridj = rids2[resnum2[j]]-1 ;
	 d = dist[i*size2+j] ;
	 rdist[i*rsize2+ridj] += d ;
	 ++count[i*rsize2+ridj] ;
      }
   }

   for (i=0;i<size1;++i)
      for (j=0;j<rsize2;++j)
	 if (count[i*rsize2+j]) {
	    rdist[i*rsize2+j] /= count[i*rsize2+j] ;
	 }

   sfree(count) ;
}

static void avg_res_atom_dists( const GrpDists *dists, const int *resnum1,
      GrpDists *resdists )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;

   int rsize1 = resdists->size1 ;
   const int *rids1 = resdists->ids1 ;

   real *rdist = resdists->dist ;
   int i,j,ridi ;
   real d ;
   int *count ;

   memset(rdist, 0, rsize1*size2*sizeof(real));
   snew(count,rsize1*size2);
   for (i=0;i<size1;++i) {
      ridi = rids1[resnum1[i]]-1 ;
      for (j=0;j<size2;++j) {
	 d = dist[i*size2+j] ;
	 rdist[ridi*size2+j] += d ;
	 ++count[ridi*size2+j] ;
      }
   }

   for (i=0;i<rsize1;++i)
      for (j=0;j<size2;++j)
	 if (count[i*size2+j]) {
	    rdist[i*size2+j] /= count[i*size2+j] ;
	 }

   sfree(count) ;
}


static void min_res_res_dists( const GrpDists *dists,  
      const int *resnum1, const int* resnum2,
      GrpDists *resdists )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;

   int rsize1 = resdists->size1 ;
   int rsize2 = resdists->size2 ;
   const int *rids1 = resdists->ids1 ;
   const int *rids2 = resdists->ids2 ;

   real *mindist = resdists->dist ;

   int i,j,ridi,ridj ;
   real d ;

   j=rsize1*rsize2 ;
   for (i=0;i<j;++i)
      mindist[i] = -1 ;

   for (i=0;i<size1;++i) {
      ridi = rids1[resnum1[i]]-1 ;
      for (j=0;j<size2;++j) {
	 ridj = rids2[resnum2[j]]-1 ;
	 d = dist[i*size2+j] ;    
	 if (mindist[ridi*rsize2+ridj]<0 || d<mindist[ridi*rsize2+ridj] ) {
	    mindist[ridi*rsize2+ridj] = d ;
	 }
      }
   }
}

static void min_atom_res_dists( const GrpDists *dists, const int* resnum2,
      GrpDists *resdists )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;  

   int rsize2 = resdists->size2 ;
   const int *rids2 = resdists->ids2 ;

   real *mindist = resdists->dist ;

   int i,j,ridj ;
   real d ;

   j=size1*rsize2 ;
   for (i=0;i<j;++i)
      mindist[i] = -1 ;

   for (i=0;i<size1;++i) {
      for (j=0;j<size2;++j) {
	 ridj = rids2[resnum2[j]]-1 ;
	 d = dist[i*size2+j] ;
	 if (mindist[i*rsize2+ridj]<0 || d<mindist[i*rsize2+ridj] ) {
	    mindist[i*rsize2+ridj] = d ;
	 }
      }
   }
}

static void min_res_atom_dists( const GrpDists *dists, const int* resnum1,
      GrpDists *resdists )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;

   int rsize1 = resdists->size1 ;
   const int *rids1 = resdists->ids1 ;

   real *mindist = resdists->dist ;

   int i,j,ridi ;
   real d ;

   j=rsize1*size2 ;
   for (i=0;i<j;++i)
      mindist[i] = -1 ;

   for (i=0;i<size1;++i) {
      ridi = rids1[resnum1[i]]-1 ;
      for (j=0;j<size2;++j) {
	 d = dist[i*size2+j] ;
	 if (mindist[ridi*size2+j]<0 || d<mindist[ridi*size2+j] ) {
	    mindist[ridi*size2+j] = d ;
	 }
      }
   }
}


typedef enum eavgsel { eAvgAtmAtm, eAvgAtmRes, eAvgResAtm, eAvgResRes } eAvgSel;
typedef enum eresval { eResMin, eResAvg } eResVal;
typedef enum edfunc  { eStep, eGauss } eDensFunc ;
typedef enum eweight { eWNone, eWAbs } eWFunc ;

// if eAvgAtmAtm make resdists = dists before or after this
static void aggregate( const GrpDists *dists, eAvgSel avsel, eResVal rval,
      const int resnums1[], const int resnums2[],
      GrpDists *resdists )
{
   switch (avsel) {
      case eAvgResRes :
	 if (eResMin == rval)
	    min_res_res_dists( dists, resnums1, resnums2, resdists ) ;
	 else // eResAvg
	    avg_res_res_dists( dists, resnums1, resnums2, resdists ) ;
	 break ;
      case eAvgResAtm :
	 if (eResMin == rval)
	    min_res_atom_dists( dists, resnums1, resdists ) ;
	 else // eResAvg
	    avg_res_atom_dists( dists, resnums1, resdists ) ;  
	 break ;
      case eAvgAtmRes :
	 if (eResMin == rval)
	    min_atom_res_dists( dists, resnums2, resdists ) ;
	 else // eResAvg
	    avg_atom_res_dists( dists, resnums2, resdists ) ;
	 break ;
      case eAvgAtmAtm : // nothing to be done
	 break ;
   }
}

static void calc_dens_cut ( const GrpDists *dists, const GrpDists *natdists,
      real rmax, real rnatmin, real *density )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;
   const real *natdist = natdists ? natdists->dist : NULL ;

   int i,j ;
   real v = 1.0 ;

   memset( density, 0, size1*sizeof(real) ) ;

   if (natdist)
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	    if ( dist[i*size2+j] < rmax && natdist[i*size2+j] > rnatmin )
	       density[i] += v ;
   }
   else // ! natdist
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	    if ( dist[i*size2+j] < rmax )
	       density[i] += v ;
   }

}

static void calc_weighted_dens_cut ( const GrpDists *dists,
      const GrpDists *natdists,
      const real weights[],
      real rmax, real rnatmin,
      real *density )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;
   const real *natdist = natdists ? natdists->dist : NULL ;

   int i,j ;
   real v ;

   memset( density, 0, size1*sizeof(real) ) ;

   if (natdist)
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	    if ( dist[i*size2+j] < rmax && natdist[i*size2+j] > rnatmin )
	    {
	       density[i] += weights[j] ;
	    }
   }
   else // ! natdist
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	    if ( dist[i*size2+j] < rmax )
	       density[i] += weights[j] ;
   }

}



static void calc_dens_gauss( const GrpDists *dists, const GrpDists *natdists,
      real sig, real rnatmin, real *density )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;
   const real *natdist = natdists ? natdists->dist : NULL ;

   int i,j ;
   real by2s2 = 0.5/(sig*sig) ;
   real siglim = 6.0 * sig ;
   real r ;
   real v ;

   memset( density, 0, size1*sizeof(real) ) ;

   if (natdist)
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	    if ( natdist[i*size2+j] > rnatmin )
	    {
	       r = dist[i*size2+j] ;
	       if (r<siglim) {
		  v = r*r*by2s2 ;
		  v = exp(-v) ;
		  density[i] += v ;
	       }
	    }
   }
   else // ! natdist
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	 { 
	    r  = dist[i*size2+j] ;
	    if (r<siglim) {
	       v = r*r*by2s2 ;
	       v = exp(-v) ;
	       density[i] += v ;
	    }
	 }
   }

}

static void calc_weighted_dens_gauss( const GrpDists *dists,
      const GrpDists *natdists,
      const real weights[],
      real sig, real rnatmin,
      real *density )
{
   int size1 = dists->size1 ;
   int size2 = dists->size2 ;
   const real *dist = dists->dist ;
   const real *natdist = natdists ? natdists->dist : NULL ;

   int i,j ;
   real by2s2 = 0.5/(sig*sig) ;
   real siglim = 6.0 * sig ;
   real r ;
   real v ;

   memset( density, 0, size1*sizeof(real) ) ;

   if (natdist)
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	    if ( natdist[i*size2+j] > rnatmin )
	    {
	       r = dist[i*size2+j] ;
	       if (r<siglim) {
		  v = r*r*by2s2 ;
		  v = exp(-v) ;
		  v *= weights[j] ;
		  density[i] += v ;
	       }
	    }
   }
   else // ! natdist
   {
      for (i=0;i<size1;++i)
	 for (j=0;j<size2;++j)
	 {
	    r  = dist[i*size2+j] ;
	    if (r<siglim) {
	       v = r*r*by2s2 ;
	       v = exp(-v) ;
	       v *= weights[j] ;
	       density[i] += v ;
	    }
	 }
   }

}

static void sum_to_residue( const real density[], const GrpDists *dists,
      const int resnums1[], real resdensity[] )
{
   int i, s1, ri ;

   s1 = dists->size1 ;
   for (i=0;i<s1;++i) {
      (resdensity)[resnums1[i]] += density[i] ;
   }
}


typedef struct {

   const t_atoms *atoms ;
   int nres ;

   eAvgSel avgsel ;
   eResVal resval ;
   eDensFunc dfunc ;
   eWFunc wfunc ;
   real drad,
	nrad ;
   gmx_bool bExNative ;
   gmx_bool bSumRes ;

   GrpDists *gdists,
	    *resdists,
	    *gdistsN,
	    *resdistsN;

   int *resnums1,
       *resatoms1,
       *resnums2,
       *resatoms2;

   real *rawdensity ;

   real *density ;
   const int *index ;
   int size ;

} DensData ;

static void map_to_opt_residues( const GrpDists *dists, const t_atoms *atoms,
      int nres, eAvgSel avgsel,
      int *resnums1[], int *resatoms1[],
      int *resnums2[], int *resatoms2[],
      GrpDists *resdists )
{
   int *resindex1 = NULL, *resindex2 = NULL ;
   int rsize1 = 0, rsize2 = 0 ;

   if (eAvgResAtm == avgsel || eAvgResRes == avgsel )
   {
      snew(*resnums1,dists->size1);
      snew(*resatoms1,nres);

      get_group_resnums( atoms, dists->members1, dists->size1,
	    *resnums1, *resatoms1 ) ;

      make_residue_index( *resatoms1, nres, &resindex1, &rsize1 ) ;
   }
   else
      *resnums1 = *resatoms1 = NULL ;

   if (eAvgAtmRes == avgsel || eAvgResRes == avgsel)
   {
      snew(*resnums2,dists->size2);
      snew(*resatoms2,nres);

      get_group_resnums( atoms, dists->members2, dists->size2,
	    *resnums2, *resatoms2 ) ;

      make_residue_index( *resatoms2, nres, &resindex2, &rsize2 ) ;
   }
   else
      *resnums2 = *resatoms2 = NULL ;

   switch (avgsel) {
      case eAvgResRes :
	 init_group_dists( resdists, nres, nres,
	       rsize1, resindex1, dists->name1,
	       rsize2, resindex2, dists->name2 ) ;
	 break ;
      case eAvgAtmRes :
	 init_group_dists( resdists, dists->nmax1, nres,
	       dists->size1, dists->members1, dists->name1,
	       rsize2, resindex2, dists->name2) ;  
	 break ;
      case eAvgResAtm :
	 init_group_dists( resdists, nres, dists->nmax2,
	       rsize1, resindex1, dists->name1,
	       dists->size2, dists->members2, dists->name2) ;
	 break ;
      case eAvgAtmAtm :
	 break ; // nothing
   }
}

static void init_density_data( DensData *data, const t_atoms *atoms,
      int nres, int natoms,
      int size1, const int idx1[], const char *name1,
      int size2, const int idx2[], const char *name2,
      eAvgSel avgsel, eResVal resval,
      eDensFunc dfunc, eWFunc wfunc,
      real drad, real nrad, gmx_bool bExNative,
      gmx_bool bSumRes )
{
   int *tmpindex ;
   data->atoms = atoms ;
   data->nres = nres ;
   data->avgsel = avgsel ;
   data->resval = resval ;
   data->dfunc = dfunc ;
   data->wfunc = wfunc ;
   data->drad = drad ;
   data->nrad = nrad ;
   data->bExNative = bExNative ;
   data->bSumRes = bSumRes ;

   snew(data->gdists,1);
   init_group_dists( data->gdists, natoms, natoms,
	 size1, idx1, name1,
	 size2, idx2, name2 ) ;
   if (eAvgAtmAtm==data->avgsel)
      data->resdists = data->gdists ;
   else
   {
      snew(data->resdists,1) ;
      map_to_opt_residues( data->gdists, data->atoms, data->nres, data->avgsel,
	    &data->resnums1, &data->resatoms1,
	    &data->resnums2, &data->resatoms2,
	    data->resdists ) ;
   }

   if (data->bExNative)
   {
      snew(data->gdistsN,1) ;
      clone_group_dists( data->gdists, data->gdistsN ) ;
      if (eAvgAtmAtm==data->avgsel)
	 data->resdistsN = data->gdistsN ;
      else
      {
	 snew(data->resdistsN,1);
	 clone_group_dists( data->resdists, data->resdistsN ) ;
      }
   }
   else // ! bExNative
   {
      data->gdistsN = NULL ;
      data->resdistsN = NULL ;
   }

   snew(data->rawdensity,data->resdists->size1) ;

   if (data->bSumRes)
   { // if bSumRes, resdists->size1 should still be atom-based
      snew(data->resnums1,data->resdists->size1) ;
      snew(data->resatoms1,data->nres) ;

      get_group_resnums( data->atoms,
	    data->resdists->members1, data->resdists->size1,
	    data->resnums1, data->resatoms1 ) ; 

      make_residue_index( data->resatoms1, data->nres,
	    &tmpindex, &data->size ) ;
      data->index = tmpindex ; // needed to keep data->index as const int*

      snew(data->density, data->size) ;
   }
   else
   {
      data->density = data->rawdensity ;
      data->index = data->resdists->members1 ;
      data->size = data->resdists->size1 ;
   }
} 

static void do_contacts(const char *fn,
      const char *fqval,const char *fqi, const char *fqnum,
      const char *fqaa, const char *fqiaa, const char *fqnaa,
      const char *fdens,
      enum eqifmt qifmt, int gnx, atom_id index[],
      int ePBC, int ncont, real *natdist,
      real cutval, gmx_bool bAbs, gmx_bool bSym,
      gmx_bool bTimes, real kappa,
      const DensData* ddata ,  
      int *resnums, int nres,
      output_env_t *oenv )
{
   FILE *qout, *iout = NULL, *nout = NULL,
	*rout = NULL, *riout = NULL, *rnout = NULL,
	*dout = NULL ;
   real t ;
   rvec *x ;
   matrix box ;
   real *aktdist ;
   int natoms, nframes,
       ci, cj, nformed ;
   real arg, tanhsum ;
   gmx_bool formed = FALSE;

   real *resrescont = NULL ;
   real *rescont = NULL ;
   real *atmcont = NULL ;
   real *denscont = NULL ;
   int *resindex = NULL ;
   int nrescont ;
   gmx_bool bResResCont = FALSE ;
   gmx_bool bResCont    = FALSE ;
   gmx_bool bAtmCont    = FALSE ;

   int ri, rj, sumi, sum ;
   real qloc ;
   t_trxstatus  *status ;

   qout = gmx_ffopen(fqval,"w") ;
   if ( qifmt != eqiNONE )
      iout = gmx_ffopen(fqi,"w") ;
   if ( fqnum )
      nout = gmx_ffopen(fqnum,"w" ) ;

   if (fqaa) {
      rout = gmx_ffopen(fqaa,"w") ;
      if ( qifmt != eqiNONE )
	 riout = gmx_ffopen(fqiaa,"w") ;
   }
   if ( fqnaa ) {
      rnout = gmx_ffopen(fqnaa,"w") ;
   }

   if ( ddata ) {
      dout = gmx_ffopen(fdens,"w") ;
   }

   natoms=read_first_x(*oenv, &status,fn,&t,&x,box);
   if (natoms == 0)
      gmx_fatal(FARGS,"No atoms in trajectory!");

   if (nout)
      bAtmCont = TRUE ;
   if (rout||rnout)
      bResCont = TRUE ;
   if (riout)
      bResResCont = TRUE ;

   if ( ddata && ddata->wfunc != eWNone )
   {
      snew(denscont,ddata->resdists->size2);
      if ( eAvgAtmRes == ddata->avgsel || eAvgResRes == ddata->avgsel )
	 bResCont = TRUE ;
      else // eAvgAtmAtm || eAvgResAtm 
	 bAtmCont = TRUE ;
   } 

   if (bResCont)
      bResResCont = TRUE ;

   snew( aktdist, ncont ) ;
   if (bResResCont)
      snew( resrescont, nres*nres ) ;
   if (bResCont)
      snew( rescont, nres ) ;
   if (bAtmCont)
      snew( atmcont, natoms ) ;

   if ( riout && ( eqiINDEX == qifmt || eqiLIST == qifmt ) )
   {
      memset(resrescont,0,nres*nres*sizeof(real)) ;
      for (ci=0;ci<ncont;++ci)
      {
	 ri = resnums[ci+ci] ;
	 rj = resnums[ci+ci+1] ;
	 resrescont[ri*nres+rj] += 1 ;
	 resrescont[rj*nres+ri] += 1 ;
      }
      nrescont=0;
      for (ri=0;ri<nres;++ri)
	 for (rj=ri+1;rj<nres;++rj)
	    if(resrescont[ri*nres+rj]>0)
	       ++nrescont ;
      snew(resindex,2*nrescont) ;
      ci=0;
      for (ri=0;ri<nres;++ri)
	 for (rj=ri+1;rj<nres;++rj)
	    if(resrescont[ri*nres+rj]>0) {
	       resindex[ci++] = ri ;
	       resindex[ci++] = rj ;
	    }
   }

   nframes = 0 ;
   do{
      ++nframes ;
      calc_cont_dists( ncont, index, x, ePBC, box, aktdist ) ;

      tanhsum = 0.0 ; 
      nformed = 0 ;
      if ( resrescont )
	 memset( resrescont, 0, nres*nres*sizeof(real) ) ;
      if ( rescont )
	 memset( rescont, 0, nres*sizeof(real) ) ;
      if ( atmcont )
	 memset( atmcont,0, natoms*sizeof(real) ) ;

      for ( ci=0 ; ci<ncont ; ++ci )
      {
	 if ( kappa > 0 )
	 {
	    arg = kappa * ( aktdist[ci] - 1.2 * natdist[ci] ) ;
	    tanhsum += 0.5 * ( 1.0 - tanh( arg ) ) ;
	 }
	 else { // kappa == 0, normal behavior
	    if ( bAbs )
	    {
	       formed = aktdist[ci] - natdist[ci] <= cutval ;
	       if ( bSym )
	       {
		  if ( formed && aktdist[ci] - natdist[ci] < -cutval )
		     formed = FALSE ;
	       }
	    }
	    else // relative
	    {
	       formed = aktdist[ci] <= (1.0+cutval)*natdist[ci] ;
	       if ( bSym )
	       {
		  if ( formed && aktdist[ci] < (1.0-cutval)*natdist[ci] )
		     formed = FALSE ;
	       }
	    } // ! if bAbs
	    if ( formed )
	    {
	       ++nformed ;

	       if ( resrescont )
	       {
		  ri = resnums[ci+ci] ;
		  rj = resnums[ci+ci+1] ;
		  ++resrescont[ri*nres+rj] ; // using resind, starts at 0
		  ++resrescont[rj*nres+ri] ;
	       }
	       if ( atmcont )
	       {
		  ++atmcont[index[ci+ci]] ;
		  ++atmcont[index[ci+ci+1]] ;
	       }

	       if ( qifmt == eqiPAIR )
	       {  // I mean 2*ci and 2*ci+1 of course
		  fprintf( iout, "%d %d\n", index[ci+ci]+1, index[ci+ci+1]+1 ) ;
	       }
	       else if ( qifmt == eqiINDEX )
		  fprintf( iout, "%d\n",ci+1 ) ;
	    }
	    if ( qifmt == eqiLIST )
	       fprintf( iout, " %d", formed ) ;
	 } // kappa == 0
      } // for ci

      if ( kappa > 0 ) {
	 if ( bTimes )
	    fprintf( qout, "%f %.3f\n",t, tanhsum ) ;
	 else
	    fprintf( qout, "%.3f\n", tanhsum ) ;
      }
      else { // kappa == 0
	 if ( bTimes )
	    fprintf( qout, "%f %d\n", t, nformed ) ;
	 else    
	    fprintf( qout, "%d\n", nformed ) ;
	 if ( qifmt == eqiPAIR )
	    fprintf( iout, "0 0\n" );
	 else if ( qifmt == eqiLIST )
	    fprintf( iout, "\n" );
	 if ( nout ) {
	    for (ci=0;ci<natoms;++ci)
	       fprintf( nout, " %.0f", atmcont[ci] ) ;
	    fprintf( nout, "\n" ) ;
	 }

	 if(rout||rnout) {
	    sum=0 ;
	    for (ri=0;ri<nres;++ri) {
	       sumi=0;
	       for (rj=0;rj<nres;++rj) { // FIXME check BUGFIX 030813, was ri+1
		  if (resrescont[ri*nres+rj]>0)
		     ++sumi ;
	       }
	       if (rescont)
		  if (sumi>0) ++rescont[ri] ; // FIXME check! addition 030713
	       if (rnout)
		  fprintf(rnout," %d",sumi) ;
	       sum += sumi ;
	    }
	    if (rnout)
	       fprintf(rnout,"\n") ;
	    if (rout)
	       fprintf(rout," %d\n",sum/2) ; // FIXME BUGFIX 030813, added /2
	 }
	 if (riout)
	 {
	    if ( eqiINDEX == qifmt ) {
	       for (ci=0;ci<nrescont;++ci) {
		  ri = resindex[ci+ci] ;
		  rj = resindex[ci+ci+1] ;
		  if ( resrescont[ri*nres+rj]>0 )
		     fprintf(riout,"%d\n",ci+1) ;
	       }
	    }
	    else if ( eqiPAIR == qifmt ) {
	       for (ri=0;ri<nres;++ri)
		  for (rj=ri+1;rj<nres;++rj)
		     if (resrescont[ri*nres+rj]>0)
			fprintf(riout,"%d %d\n",ri+1,rj+1) ;
	       fprintf(riout,"0 0\n") ;
	    }
	    else if ( eqiLIST == qifmt ) {
	       for (ci=0;ci<nrescont;++ci) {
		  formed = 0 ;
		  ri = resindex[ci+ci] ;
		  rj = resindex[ci+ci+1] ;
		  if ( resrescont[ri*nres+rj]>0 )
		     formed = 1 ;
		  fprintf(riout," %d",formed) ;
	       }
	       fprintf(riout,"\n") ;
	    }
	 }//if riout

      } // kappa == 0

      if (ddata)
      {
	 calc_group_dists( ddata->gdists, x, ePBC, box ) ;

	 if (eAvgAtmAtm != ddata->avgsel)
	    aggregate( ddata->gdists, ddata->avgsel, ddata->resval,
		  ddata->resnums1, ddata->resnums2,
		  ddata->resdists ) ;

	 if (denscont) {
	    if (eAvgAtmRes == ddata->avgsel || eAvgResRes == ddata->avgsel)
	    {
	       for (ri=0;ri<ddata->resdists->size2;++ri)
		  denscont[ri] = rescont[ddata->resdists->members2[ri]] ;
	    }
	    else //atmcont
	    {
	       for (ri=0;ri<ddata->resdists->size2;++ri)
		  denscont[ri] = atmcont[ddata->resdists->members2[ri]] ;
	    }
	 }

	 switch (ddata->dfunc)
	 {
	    case eStep:
	       if (denscont)
		  calc_weighted_dens_cut( ddata->resdists, ddata->resdistsN,
			denscont, ddata->drad, ddata->nrad,
			ddata->rawdensity ) ;
	       else
		  calc_dens_cut( ddata->resdists, ddata->resdistsN,
			ddata->drad, ddata->nrad, ddata->rawdensity ) ;
	       break ;
	    case eGauss:
	       if (denscont)
		  calc_weighted_dens_gauss( ddata->resdists, ddata->resdistsN,
			denscont, ddata->drad, ddata->nrad,
			ddata->rawdensity ) ;
	       else
		  calc_dens_gauss( ddata->resdists, ddata->resdistsN,
			ddata->drad, ddata->nrad, ddata->rawdensity ) ;
	       break ;
	 }
	 if (ddata->bSumRes)
	    sum_to_residue( ddata->rawdensity,
		  ddata->resdists,
		  ddata->resnums1,
		  ddata->density ) ;

	 if (bTimes)
	    fprintf(dout,"%f",t) ;
	 for (ci=0;ci<ddata->size;++ci)
	    fprintf(dout," %.3f",ddata->density[ci]) ;
	 fprintf(dout,"\n") ;

      } // ddata

   }
   while( read_next_x(*oenv,status,&t,x,box) ) ;
   close_trj(status);

   gmx_ffclose(qout);
   if (iout)
      gmx_ffclose(iout);
   if (nout)
      gmx_ffclose(nout);
   if (rout)
      gmx_ffclose(rout);
   if (riout)
      gmx_ffclose(riout);
   if (rnout)
      gmx_ffclose(rnout);
   if(dout)
      gmx_ffclose(dout);

   if (resindex) sfree(resindex) ;
   if (denscont) sfree(denscont) ;
   if (resrescont) sfree(resrescont) ;
   if (rescont) sfree(rescont) ;
   if (atmcont) sfree(atmcont) ;
}

void filter_cmap( atom_id filtgrp1[], int n1, atom_id filtgrp2[], int n2,
      int natoms, atom_id **pairs, int *size, real **dists )
{
   atom_id *newpairs = NULL ;
   real *newdists = NULL ;
   int newsize = 0 ;
   int *inone, *intwo, *use ;
   snew(inone,natoms);
   snew(intwo,natoms);
   snew(use,*size/2);
   int i,j,ci,cj ;

   for (i=0;i<n1;++i)
      inone[filtgrp1[i]] = 1 ;

   for (i=0;i<n2;++i)
      intwo[filtgrp2[i]] = 1 ;

   i=0 ;
   for (ci=0;ci<*size/2;++ci)
   {
      if ((inone[(*pairs)[i]]   && intwo[(*pairs)[i+1]] )
	    ||(inone[(*pairs)[i+1]] && intwo[(*pairs)[i]]   ))
      {
	 use[ci]=1 ;
	 newsize+=2 ;
      }
      i+=2;
   }

   i=0 ; j=0 ; cj = 0 ;
   snew(newpairs,newsize);
   if (dists)
      snew(newdists,newsize/2) ;
   for (ci=0;ci<*size/2;++ci)
   {
      if (use[ci])
      {
	 newpairs[j] = (*pairs)[i] ;
	 newpairs[j+1] = (*pairs)[i+1] ;
	 if (dists)
	    newdists[cj] = (*dists)[ci] ;
	 j+=2 ;
	 ++cj ;
      }
      i+=2 ;
   }

   sfree(*pairs) ;
   *size = newsize ;
   *pairs = newpairs ;
   if (dists)
      *dists = newdists ;

   sfree(use) ;
   sfree(inone);
   sfree(intwo);
}


int gmx_kuh(int argc,char *argv[])
{
   static const char *desc[] = {
      "g_kuh is very neat."
   };
   static const char *bugs[] = {
      "r u kidding?"
   };

   static const char *qifmtsel[] = { 
      NULL, "none", "index", "pair", "list", NULL
   };

   static const char *avgsels[] = {
      NULL, "atmatm", "atmres", "resatm", "resres", NULL
   };

   static const char *resvals[] = {
      NULL, "min", "avg", NULL
   }; 

   static const char *densfunc[] = {
      NULL, "step", "gauss", NULL
   };

   static const char *weightfunc[] = {
      NULL, "none", "abs", NULL
   } ;

   static gmx_bool bAbs = TRUE, bSym = TRUE, bTimes = FALSE ;
   static real cutval ;
   static real kappa_tanh = 0.0 ;
   static gmx_bool bFilt = FALSE ;
   static gmx_bool bPBC = FALSE ;

   static gmx_bool bGroups = FALSE, bG1std = TRUE, bG2std = TRUE ;
   static gmx_bool bSumRes = FALSE ;
   static gmx_bool bExNative = FALSE ;
   static real drad = 0.6, nrad = 0.0 ;

   //  static gmx_bool bNormRes = FALSE ;

   t_pargs pa[] = {
      { "-cut", FALSE, etREAL, {&cutval},
	 "Contact cutoff" },
      { "-abscut", FALSE, etBOOL, {&bAbs},
	 "use absolute cutoff (instead of relative)" },
      { "-shortcut", FALSE, etBOOL, {&bSym},
	 "use cutoff also at short distances" },
      { "-qiformat", FALSE, etENUM, {qifmtsel},
	 "for individual contacts" },
      { "-times", FALSE, etBOOL, {&bTimes},
	 "print times into output file"   },
      { "-pbc", FALSE, etBOOL, {&bPBC}, "force PBC (XYZ), or NONE" },
      { "-kappa", FALSE, etREAL, {&kappa_tanh}, "kappa for continuous Q" },
      { "-groups", FALSE, etBOOL, {&bGroups}, "define groups" },
      { "-g1std", FALSE, etBOOL, {&bG1std}, "use built-in groups for group 1" },
      { "-g2std", FALSE, etBOOL, {&bG2std}, "use built-in groups for group 2" },
      { "-filtmap", FALSE, etBOOL, {&bFilt}, "filter cmap by groups" },
      { "-avgdens", FALSE, etENUM, {avgsels}, "aggregating choice for density" },
      { "-resval", FALSE, etENUM, {resvals}, "aggregation method" },
      { "-weight", FALSE, etENUM, {weightfunc}, "weight with contacts?" },
      { "-sum2res", FALSE, etBOOL, {&bSumRes}, "sum up residue densities" },
      { "-dfunc", FALSE, etENUM, {densfunc}, "density function" },
      { "-drad", FALSE, etREAL, {&drad}, "radius for density" },
      { "-nrad", FALSE, etREAL, {&nrad}, "native radius for density" },
      { "-exnative", FALSE, etBOOL, {&bExNative}, "exclude native neighbors" }
      //    ,{ "-normres", FALSE, etBOOL, {&bNormRes}, "normalize Qres" }
   };
   FILE      *fp;
   char      *grpname ;
   int       gnx;
   atom_id   *index;
   matrix    box;

   char      title[STRLEN]; 
   int       natoms, ePBC ;
   t_atoms   atoms ;
   rvec      *xn = NULL ;
   int       ncont ;
   real      *natdist ;
   gmx_bool  bCalcDist = TRUE ;

   enum      eqifmt qifmt = eqiNONE ;

   output_env_t oenv;

   int       isize_ctr, isize_ext ;
   atom_id   *idx_ctr, *idx_ext ;
   char      *grpname_ctr, *grpname_ext ;
   int i ;
   eAvgSel   avgsel = eAvgAtmAtm ;
   eResVal   resval = eResMin ;
   eDensFunc dfunc = eStep ;
   eWFunc    wfunc = eWNone ;
   DensData  ddata ; 
   gmx_bool  bDensity = FALSE;

   const char *qaaname, *resname ;
   gmx_bool  bRes ;
   int       *resnums = NULL ;

   t_filenm fnm[] = {
      { efSTX, "-s", "native", ffREAD  },
      { efTRX, "-f", NULL, ffREAD  },
      { efNDX, "-n", NULL, ffOPTRD},
      { efDAT, "-nc", "contacts", ffOPTRD},
      { efNDX, "-g",  "groups",   ffOPTRD },
      { efOUT, "-o",  "qvals",    ffWRITE },
      { efOUT, "-i",  "qimap",    ffOPTWR },
      { efOUT, "-on", "qatom",    ffOPTWR },
      { efOUT, "-r",  "qaa",      ffOPTWR },
      { efOUT, "-ri", "qiaa",     ffOPTWR },
      { efOUT, "-rn", "qres",     ffOPTWR },
      { efOUT, "-d",  "density",  ffOPTWR }
   };
#define NFILE asize(fnm)

   //  CopyRight(stderr,argv[0]);
   parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_BE_NICE , // PCA_CAN_VIEW
	 NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs,
	 &oenv);

   switch (qifmtsel[0][0]) {
      case 'n':
	 qifmt = eqiNONE ;
	 break;
      case 'i':
	 qifmt = eqiINDEX ;
	 break;
      case 'p':
	 qifmt = eqiPAIR ;
	 break;
      case 'l':
	 qifmt = eqiLIST ;
	 break;
   }

   get_stx_coordnum(ftp2fn(efSTX,NFILE,fnm),&natoms);
   init_t_atoms(&atoms,natoms,TRUE);
   snew(xn,natoms);
   read_stx_conf( ftp2fn(efSTX,NFILE,fnm), title, &atoms, xn, NULL, &ePBC, box);
   if (opt2parg_bSet("-pbc",asize(pa),pa)) {
      if (bPBC) {
	 fprintf(stderr,"manually setting PBC to XYZ\n");
	 ePBC = epbcXYZ ;
      }
      else {
	 fprintf(stderr,"manually setting PBC to NONE\n");
	 ePBC = epbcNONE ;
      }
   }

   bDensity = opt2bSet("-d",NFILE,fnm) ;
   if (bGroups && !(bDensity || bFilt)) {
      fprintf(stderr,"no use for -groups without -d (density) or -filtmap, deactivating\n" );
      bGroups = FALSE ;
   }
   if (bFilt && ! bGroups) {
      fprintf(stderr,"-filtmap selected without -groups, activating\n" );
      bGroups = TRUE ;
   }

   if ( bGroups )
   {
      get_index( &atoms, bG1std ? NULL : opt2fn("-g",NFILE,fnm), 1,
	    &isize_ctr, &idx_ctr, &grpname_ctr );

      get_index( &atoms, bG2std ? NULL : opt2fn("-g",NFILE,fnm), 1,
	    &isize_ext, &idx_ext, &grpname_ext );
   }
   else if ( bDensity )// only necessary for density, not for filtering cmap
   {
      snew(idx_ctr,natoms) ;
      for (i=0;i<natoms;++i)
	 idx_ctr[i]=i;
      idx_ext = idx_ctr ;
      isize_ctr = isize_ext = natoms ;
      grpname_ctr = grpname_ext = strdup("all"); // same name ...
   }

   if ( bDensity )
   {

      if ('a' == avgsels[0][0])
	 if ('a' == avgsels[0][3])
	    avgsel = eAvgAtmAtm ;
	 else
	    avgsel = eAvgAtmRes ;
      else // 'r'
	 if ('a' == avgsels[0][3])
	    avgsel = eAvgResAtm ;
	 else
	    avgsel = eAvgResRes ;

      switch (resvals[0][0]) {
	 case 'm':
	    resval = eResMin ;
	    break ;
	 case 'a':
	    resval = eResAvg ;
	    break ;
      }

      switch (densfunc[0][0]) {
	 case 's':
	    dfunc = eStep ;
	    break ;
	 case 'g':
	    dfunc = eGauss ;
	    break ;
      }

      switch (weightfunc[0][0]) {
	 case 'n':
	    wfunc = eWNone ;
	    break ;
	 case 'a':
	    wfunc = eWAbs ;
	    break ;
      }

      if (bSumRes)
	 if (eAvgResAtm == avgsel || eAvgResRes == avgsel) {
	    fprintf( stderr,
		  "already aggregating to residue, deactivating -sum2res\n");
	    bSumRes = FALSE ;
	 }  

      init_density_data( &ddata, &atoms, atoms.nres, natoms,
	    isize_ctr, idx_ctr, grpname_ctr,
	    isize_ext, idx_ext, grpname_ext,
	    avgsel, resval, dfunc, wfunc,
	    drad, nrad, bExNative, bSumRes ) ;

      if (bExNative)
      {
	 calc_group_dists( ddata.gdistsN, xn, ePBC, box ) ;

	 if (eAvgAtmAtm != avgsel)
	    aggregate( ddata.gdistsN, avgsel, resval,
		  ddata.resnums1, ddata.resnums2,
		  ddata.resdistsN ) ;
      } // bExNative

      fprintf(stdout,"# DENSITY for indices [ %d ]:\n",ddata.size);
      for (i=0;i<ddata.size;++i)
	 fprintf(stdout," %5d", ddata.index[i]+1);
      fprintf(stdout,"\n");
   } // if bDensity

   qaaname = NULL ;
   resname = NULL ;
   bRes = opt2bSet( "-r", NFILE, fnm )
      || opt2bSet( "-ri", NFILE, fnm )
      || opt2bSet( "-rn", NFILE, fnm ) ;
   bRes = bRes || (eWNone != wfunc) ;

   bCalcDist = TRUE ;
   if ( opt2bSet( "-nc", NFILE, fnm ) )
      if ( opt2bSet( "-n", NFILE, fnm ) ) {
	 fprintf(stderr, "both -n and -nc given, using -n with calc. dists.\n");
      }
      else {
	 fprintf(stderr, "-nc given, reading contact distances\n") ;
	 bCalcDist = FALSE ;
      }
   else
      if ( ! opt2bSet( "-n", NFILE, fnm ) )
      {
	 fprintf(stderr, "no contacts defined, quitting\n") ;
	 exit(1) ; // could actually be tolerated for -d only
      }

   if (bCalcDist)
   {
      rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
      if ( !even(gnx) )
	 fprintf(stderr,"WARNING: odd number of atoms (%d) in group!\n",gnx);
   }
   else
   {
      if (! read_cont_dists( opt2fn("-nc",NFILE,fnm), &gnx, &index, &natdist ) )
      {
	 fprintf(stderr,"ERROR reading contacs, quitting\n") ;
	 exit(1);
      }
   }
   if (bFilt)
   {
      fprintf(stderr,"filtering contact map with groups\n");
      filter_cmap( idx_ctr, isize_ctr, idx_ext, isize_ext,
	    natoms, &index, &gnx, (bCalcDist ? NULL : &natdist) ) ;
   }
   ncont = gnx/2 ;
   fprintf(stderr,"Will gather information on %d contacts\n",ncont);
   if (bCalcDist) {
      snew( natdist, ncont ) ;
      fprintf(stderr,"calculating contact distances from structure\n") ;
      calc_cont_dists( ncont, index, xn, ePBC, box, natdist ) ;
   }


   if ( kappa_tanh > 0 ) {
      fprintf(stderr, "kappa = %f given for continuous Q\nOther options reset to defaults:\n -shortcut=FALSE,\n -abscut=FALSE,\n -cut=1.2,\n -qimap=FALSE\n", kappa_tanh ) ;
      bSym=FALSE ;
      bAbs=FALSE ;
      cutval=1.2 ;
      qifmt = eqiNONE ;
   }

   if (bRes) {
      snew( resnums, 2*ncont ) ;
      get_contact_resnums( &atoms, index, ncont, resnums ) ;
   }

   do_contacts(ftp2fn(efTRX,NFILE,fnm),
	 opt2fn("-o",NFILE,fnm), opt2fn("-i",NFILE,fnm),
	 opt2fn_null("-on",NFILE,fnm),
	 opt2fn_null("-r",NFILE,fnm), opt2fn("-ri",NFILE,fnm),
	 opt2fn_null("-rn",NFILE,fnm),
	 opt2fn("-d",NFILE,fnm),
	 qifmt, gnx, index, ePBC, ncont, natdist, cutval,
	 bAbs, bSym, bTimes, kappa_tanh,
	 (bDensity ? &ddata : NULL),
	 resnums, atoms.nres, &oenv );

   gmx_thanx(stderr);

   return 0;
}

