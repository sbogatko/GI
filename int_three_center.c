
#include	<math.h>
#include	<stdio.h>
#include	"int_three_center.h"

static double pi=3.1415926535897931159979635;
static double xi,zi,xi_abc,zi_abc;
static double oneoverxi,oneoverxi_abc;
static double R[3][3],G[3];

/* three center integral (s_A|s_C|s_C) = (0_A|0_B|0_C) */
static double sa_sc_sb;


extern double three_center_sub();

double int_three_center(int nA[], double xiA, double A[], 
                        int nC[], double xiC, double C[], 
                        int nB[], double xiB, double B[])
{
    int i;
    int n[3][3];
    double sumAB,sumPC,c1,P[3];

    double value = 0.0;


    /* set up basic parameters */
    for (i=0;i<3;++i) R[i][0] = A[i];
    for (i=0;i<3;++i) R[i][1] = C[i];
    for (i=0;i<3;++i) R[i][2] = B[i];

    for (i=0;i<3;++i) n[i][0] = nA[i];
    for (i=0;i<3;++i) n[i][1] = nC[i];
    for (i=0;i<3;++i) n[i][2] = nB[i];


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
    value = three_center_sub(n);

    return value;
}


double three_center_sub(int n[3][3])
{
   int    xyz,acb,i,nsumall,done;
   double value,value2;
   
   nsumall = 0;
   for (acb=0; acb<3; ++acb)
   for (xyz=0; xyz<3; ++xyz)
      nsumall += n[xyz][acb];
/*
   printf("hera nA=(%d,%d,%d) nC=(%d,%d,%d) nB=(%d,%d,%d)\n",n[0][0],n[1][0],n[2][0],
                                                             n[0][1],n[1][1],n[2][1],
                                                             n[0][2],n[1][2],n[2][2]);
*/

   if (nsumall == 0)
   {
      value = sa_sc_sb;
   }
   else
   {
      done = 0;
      acb  = 0;
      xyz  = 0;
      while (!done)
      {
         if (acb>2)
         {
             printf("Should not happen\n\n");
             return 0.0;
         }

         if (n[xyz][acb] > 0)
         {
            n[xyz][acb] -= 1;
            value  = (G[xyz] - R[xyz][acb])*three_center_sub(n);

            value2 = 0.0;
            for (i=0; i<3; ++i)
            {
               n[xyz][i] -= 1;
               if (n[xyz][i] >=0) 
                  value2 += (n[xyz][i]+1)*three_center_sub(n);
               n[xyz][i] += 1;
            }
            value += 0.5*oneoverxi_abc*value2;

            n[xyz][acb] += 1;
            done = 1;
         }

         xyz = (xyz+1)%3;
         if (xyz==0) ++acb;
      }
   }

   return value;
}

