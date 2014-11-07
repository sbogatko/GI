
#include	<math.h>
#include	<stdio.h>
#include	"three_center.h"

static double pi=3.1415926535897931159979635;
static double xiA,xiB,xiC,xi,zi,xi_abc,zi_abc;
static double oneoverxi,oneoverxi_abc;
static double A[3],B[3],C[3],P[3],G[3];

/* three center integral (s_A|s_C|s_C) = (0_A|0_B|0_C) */
static double sa_sc_sb;


extern double three_center_sub();

double three_center(int n1[], double xi1, double r1[], 
                    int n2[], double xi2, double r2[], 
                    int n3[], double xi3, double r3[])
{
    int i;
    double sumAB,sumPC,c1;

    double value = 0.0;


    /* set up basic parameters */
    xiA = xi1; xiC = xi2; xiB = xi3;
    for (i=0;i<3;++i) A[i] = r1[i];
    for (i=0;i<3;++i) C[i] = r2[i];
    for (i=0;i<3;++i) B[i] = r3[i];

    xi            = xiA + xiB;
    xi_abc        = xi  + xiC;
    oneoverxi     = 1.0/xi;
    oneoverxi_abc = 1.0/xi_abc;

    zi     = xiA*xiB*oneoverxi;
    zi_abc = xi *xiC*oneoverxi_abc;

    for (i=0; i<3; ++i) P[i] = (xiA*A[i]+xiB*B[i])*oneoverxi;
    for (i=0; i<3; ++i) G[i] = (xi *P[i]+xiC*C[i])*oneoverxi_abc;


    /* compute basis three center integral (s_A|s_C|s_B) */
    sumAB = 0.0; for (i=0; i<3; ++i) sumAB += (A[i]-B[i])*(A[i]-B[i]); sumAB *= zi;
    sumPC = 0.0; for (i=0; i<3; ++i) sumPC += (P[i]-C[i])*(P[i]-C[i]); sumPC *= zi_abc;
    c1 = (pi*oneoverxi_abc); c1 = sqrt(c1*c1*c1);
    sa_sc_sb = c1*exp(-sumAB)*exp(-sumPC);


    /* call recursive three center integral formula */
    value = three_center_sub(n1[0],n1[1],n1[2],
                             n2[0],n2[1],n2[2],
                             n3[0],n3[1],n3[2]);

    return value;
}


double three_center_sub(int nax, int nay, int naz,
                        int ncx, int ncy, int ncz,
                        int nbx, int nby, int nbz)
{
   int    nsumx,nsumy,nsumz,nsumall;
   double value;
   
   nsumx = nax + nbx + ncx;
   nsumy = nay + nby + ncy;
   nsumz = naz + nbz + ncz;
   nsumall = nsumx + nsumy + nsumz;
   printf("hera n1=(%d,%d,%d) n2=(%d,%d,%d) n3=(%d,%d,%d)\n",nax,nay,naz,ncx,ncy,ncz,nbx,nby,nbz);

   if ((nax < 0) || (nay < 0) || (naz < 0) ||
       (nbx < 0) || (nby < 0) || (nbz < 0) ||
       (ncx < 0) || (ncy < 0) || (ncz < 0))
   {
      value = 0.0;
   }
   else if (nsumall == 0)
   {
      value = sa_sc_sb;
   }

   /* decrement by (1,0,0) (0,0,0) (0,0,0) */
   else if (nax > 0)
   {
      value  = (G[0] - A[0])              *three_center_sub(nax-1,nay,naz, ncx,ncy,ncz,   nbx,nby,nbz);
      value += 0.5*oneoverxi_abc*( (nax-1)*three_center_sub(nax-2,nay,naz, ncx,ncy,ncz,   nbx,nby,nbz)
                                 + (nbx)  *three_center_sub(nax-1,nay,naz, ncx,ncy,ncz,   nbx-1,nby,nbz)
                                 + (ncx)  *three_center_sub(nax-1,nay,naz, ncx-1,ncy,ncz, nbx,nby,nbz));
   }

   /* decrement by (0,1,0) (0,0,0) (0,0,0) */
   else if (nay > 0)
   {
      value  = (G[1] - A[1])              *three_center_sub(nax,nay-1,naz, ncx,ncy,ncz,   nbx,nby,nbz);
      value += 0.5*oneoverxi_abc*( (nay-1)*three_center_sub(nax,nay-2,naz, ncx,ncy,ncz,   nbx,nby,nbz)
                                 + (nby)  *three_center_sub(nax,nay-1,naz, ncx,ncy,ncz,   nbx,nby-1,nbz)
                                 + (ncy)  *three_center_sub(nax,nay-1,naz, ncx,ncy-1,ncz, nbx,nby,nbz));
   }



   /* decrement by (0,0,1) (0,0,0) (0,0,0) */
   else if (naz > 0)
   {
      value  = (G[2] - A[2])              *three_center_sub(nax,nay,naz-1, ncx,ncy,ncz,   nbx,nby,nbz);
      value += 0.5*oneoverxi_abc*( (naz-1)*three_center_sub(nax,nay,naz-2, ncx,ncy,ncz,   nbx,nby,nbz)
                                 + (nbz)  *three_center_sub(nax,nay,naz-1, ncx,ncy,ncz,   nbx,nby,nbz-1)
                                 + (ncz)  *three_center_sub(nax,nay,naz-1, ncx,ncy,ncz-1, nbx,nby,nbz));
   }


   /* decrement by (0,0,0) (1,0,0) (0,0,0) */
   else if (ncx > 0)
   {
      value  = (G[0] - C[0])              *three_center_sub(nax,nay,naz,   ncx-1,ncy,ncz,   nbx,nby,nbz);
      value += 0.5*oneoverxi_abc*( (nax)  *three_center_sub(nax-1,nay,naz, ncx-1,ncy,ncz,   nbx,nby,nbz)
                                 + (nbx-1)*three_center_sub(nax,nay,naz,   ncx-1,ncy,ncz,   nbx-1,nby,nbz)
                                 + (ncx-1)*three_center_sub(nax,nay,naz,   ncx-2,ncy,ncz,   nbx,nby,nbz));
   }

   /* decrement by (0,0,0) (0,1,0) (0,0,0) */
   else if (ncy > 0)
   {
      value  = (G[1] - C[1])              *three_center_sub(nax,nay,naz,   ncx,ncy-1,ncz,   nbx,nby,nbz);
      value += 0.5*oneoverxi_abc*( (nay)  *three_center_sub(nax,nay-1,naz, ncx,ncy-1,ncz,   nbx,nby,nbz)
                                 + (nby)  *three_center_sub(nax,nay,naz,   ncx,ncy-1,ncz,   nbx,nby-1,nbz)
                                 + (ncy-1)*three_center_sub(nax,nay,naz,   ncx,ncy-2,ncz,   nbx,nby,nbz));
   }

   /* decrement by (0,0,0) (0,0,1) (0,0,0) */
   else if (ncz > 0)
   {
      value  = (G[2] - C[2])              *three_center_sub(nax,nay,naz,   ncx,ncy,ncz-1, nbx,nby,nbz);
      value += 0.5*oneoverxi_abc*( (naz)  *three_center_sub(nax,nay,naz-1, ncx,ncy,ncz-1, nbx,nby,nbz)
                                 + (nbz)  *three_center_sub(nax,nay,naz,   ncx,ncy,ncz-1, nbx,nby,nbz-1)
                                 + (ncz-1)*three_center_sub(nax,nay,naz,   ncx,ncy,ncz-2, nbx,nby,nbz));
   }

   /* decrement by (0,0,0) (0,0,0) (1,0,0) */
   else if (nbx > 0)
   {
      value  = (G[0] - B[0])              *three_center_sub(nax,nay,naz,   ncx,ncy,ncz,   nbx-1,nby,nbz);
      value += 0.5*oneoverxi_abc*( (nax)  *three_center_sub(nax-1,nay,naz, ncx,ncy,ncz,   nbx-1,nby,nbz)
                                 + (nbx-1)*three_center_sub(nax,nay,naz,   ncx,ncy,ncz,   nbx-2,nby,nbz)
                                 + (ncx)  *three_center_sub(nax,nay,naz,   ncx-1,ncy,ncz, nbx-1,nby,nbz));
   }

   /* decrement by (0,0,0) (0,0,0) (0,1,0) */
   else if (nby > 0)
   {
      value  = (G[1] - B[1])              *three_center_sub(nax,nay,naz,   ncx,ncy,ncz,   nbx,nby-1,nbz);
      value += 0.5*oneoverxi_abc*( (nay)  *three_center_sub(nax,nay-1,naz, ncx,ncy,ncz,   nbx,nby-1,nbz)
                                 + (nby-1)*three_center_sub(nax,nay,naz,   ncx,ncy,ncz,   nbx,nby-2,nbz)
                                 + (ncy)  *three_center_sub(nax,nay,naz,   ncx,ncy-1,ncz, nbx,nby-1,nbz));
   }

   /* decrement by (0,0,0) (0,0,0) (0,0,1) */
   else if (nbz > 0)
   {
      value  = (G[2] - B[2])              *three_center_sub(nax,nay,naz,   ncx,ncy,ncz,   nbx,nby,nbz-1);
      value += 0.5*oneoverxi_abc*( (naz)  *three_center_sub(nax,nay,naz-1, ncx,ncy,ncz,   nbx,nby,nbz-1) 
                                 + (nbz-1)*three_center_sub(nax,nay,naz,   ncx,ncy,ncz,   nbx,nby,nbz-2) 
                                 + (ncz)  *three_center_sub(nax,nay,naz,   ncx,ncy,ncz-1, nbx,nby,nbz-1));
   }


   else
   {
      printf("should not be here!!!\n");
   }

   return value;
}

