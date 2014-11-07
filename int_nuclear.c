
#include	<math.h>
#include	<stdio.h>
#include	"igamma.h"
#include	"int_nuclear.h"

#define	M_MAX	200

static double pi=3.1415926535897931159979635;
static double xi,zi,xiarray[2];
static double oneoverxi,U,Uconst;
static double R[3][3],P[3];

/* two center integral (s_A|s_B) = (0_A|0_B) */
static double sa_sb;

/* two center integral (s_A|1/|r-C||s_B)^(m) */
static double sa_V_sb[M_MAX];
static int    mset;


extern double nuclear_sub();

double int_nuclear(int nA[], double xiA, double A[], 
                                         double C[],
                   int nB[], double xiB, double B[]) 
{
    int i,n[3][2];
    double sumAB,sumPC,c1;

    double value = 0.0;


    /* set up basic parameters */
    for (i=0;i<3;++i) R[i][0] = A[i];
    for (i=0;i<3;++i) R[i][1] = B[i];
    for (i=0;i<3;++i) R[i][2] = C[i];

    for (i=0;i<3;++i) n[i][0] = nA[i];
    for (i=0;i<3;++i) n[i][1] = nB[i];

    xiarray[0] = xiA;
    xiarray[1] = xiB;
    xi            = xiA + xiB;
    oneoverxi     = 1.0/xi;

    zi     = xiA*xiB*oneoverxi;

   
    for (i=0;i<3;++i) P[i] = (xiA*A[i] + xiB*B[i])*oneoverxi;

    U = 0.0; for (i=0; i<3; ++i) U += (P[i]-C[i])*(P[i]-C[i]); U *= xi;
    Uconst = 2.0*sqrt(xi/pi);

    /* compute basis two center integral (s_A|s_B) */
    sumAB = 0.0; for (i=0; i<3; ++i) sumAB += (A[i]-B[i])*(A[i]-B[i]); sumAB *= zi;
    c1 = (pi*oneoverxi); c1 = sqrt(c1*c1*c1);
    sa_sb = c1*exp(-sumAB);

    mset = 0; sa_V_sb[mset] = Uconst*sa_sb*igamma(mset,U);
//    ++mset;   sa_V_sb[mset] = Uconst*sa_sb*igamma(mset,U);


    /* call recursive two center integral formula */
    value = nuclear_sub(0,n);
                   
    return value;
}


double nuclear_sub(int m, int n[3][2])
{
   int    i,xyz,ab,nsumall,done;
   double value,value2;

   nsumall = 0;
   for (ab=0;  ab<2;  ++ab)
   for (xyz=0; xyz<3; ++xyz)
      nsumall += n[xyz][ab];

   if (nsumall == 0)
   {
      while (mset<m)
      {
        ++mset;
        sa_V_sb[mset] = Uconst*sa_sb*igamma(mset,U);
      }
      value = sa_V_sb[m];
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
            printf("nuclear_sub: Should not happen\n\n");
            return 0.0;
         }
         if (n[xyz][ab] > 0)
         {

            n[xyz][ab] -= 1; 
            value  = (P[xyz] - R[xyz][ab])*nuclear_sub(m,  n);
            value -= (P[xyz] - R[xyz][2] )*nuclear_sub(m+1,n);

            value2 = 0.0;
            for (i=0; i<2; ++i)
            {
               n[xyz][i] -= 1;
               if (n[xyz][i] >= 0)
                  value2 += (n[xyz][i]+1)*(nuclear_sub(m,n)-nuclear_sub(m+1,n));
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


