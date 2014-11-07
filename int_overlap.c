
#include	<math.h>
#include	<stdio.h>
#include	"int_overlap.h"

static double pi=3.1415926535897931159979635;
static double xi,zi;
static double oneoverxi;
static double R[3][2],P[3];

/* two center integral (s_A|s_B) = (0_A|0_B) */
static double sa_sb;


extern double overlap_sub();

double int_overlap(int nA[], double xiA, double A[], 
                   int nB[], double xiB, double B[]) 
{
    int i,n[3][2];
    double sumAB,sumPC,c1;

    double value = 0.0;


    /* set up basic parameters */
    for (i=0;i<3;++i) R[i][0] = A[i];
    for (i=0;i<3;++i) R[i][1] = B[i];

    for (i=0;i<3;++i) n[i][0] = nA[i];
    for (i=0;i<3;++i) n[i][1] = nB[i];

    xi            = xiA + xiB;
    oneoverxi     = 1.0/xi;

    zi     = xiA*xiB*oneoverxi;

   
    for (i=0;i<3;++i) P[i] = (xiA*A[i] + xiB*B[i])*oneoverxi;


    /* compute basis two center integral (s_A|s_B) */
    sumAB = 0.0; for (i=0; i<3; ++i) sumAB += (A[i]-B[i])*(A[i]-B[i]); sumAB *= zi;
    c1 = (pi*oneoverxi); c1 = sqrt(c1*c1*c1);
    sa_sb = c1*exp(-sumAB);


    /* call recursive two center integral formula */
    value = overlap_sub(n);
                   
    return value;
}


double overlap_sub(int n[3][2])
{
   int    i,xyz,ab,nsumall,done;
   double value,value2;
   
   nsumall = 0;
   for (ab=0;  ab<2;  ++ab)
   for (xyz=0; xyz<3; ++xyz)
      nsumall += n[xyz][ab];

   if (nsumall == 0)
   {
      value = sa_sb;
   }
   else 
   {
      done = 0;
      xyz  = 0;
      ab   = 0;
      while (!done)
      {
         if (ab>1)
         {
            printf("Should not happen\n\n");
            return 0.0;
         }

         if (n[xyz][ab] > 0)
         {
            n[xyz][ab] -= 1;
            value = (P[xyz] - R[xyz][ab])*overlap_sub(n);

            value2 = 0.0;
            for (i=0; i<2; ++i)
            {
               n[xyz][i] -= 1;
               if (n[xyz][i] >= 0)
                  value2 += (n[xyz][i]+1)*overlap_sub(n);
               n[xyz][i] += 1;
            }
            value += 0.5*oneoverxi*value2;

            n[xyz][ab] += 1;
            done = 1;
         }
         
         xyz = (xyz+1)%3;
         if (xyz==0) ++ab;
      }
   }

   return value;
}

