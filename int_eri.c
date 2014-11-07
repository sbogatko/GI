
#include	<math.h>
#include	<stdio.h>
#include	"igamma.h"
#include	"int_nuclear.h"

#define	M_MAX	200

static double pi=3.1415926535897931159979635;
static double xi[2],zi[2];
static double rho;
static double oneoverxi[2],oneoverxiall,T,Tconst;
static double R[3][2][2],P[3][2],W[3];

/* two center integrals (s_A|s_B) = (0_A|0_B), (s_C|s_D) = (0_C|0_D) */
static double sa_sb[2];

/* two center integral (s_A|1/|r-C||s_B)^(m) */
static double sasb_V_scsd[M_MAX];
static int    mset;


extern double eri_sub();

double int_eri(int nA[], double xiA, double A[], 
               int nB[], double xiB, double B[], 
               int nC[], double xiC, double C[], 
               int nD[], double xiD, double D[])
{
    int i,n[3][2][2];
    double sumAB,sumCD,c1,tmp;

    double value = 0.0;


    /* set up basic parameters */
    for (i=0;i<3;++i) R[i][0][0] = A[i];
    for (i=0;i<3;++i) R[i][1][0] = B[i];
    for (i=0;i<3;++i) R[i][0][1] = C[i];
    for (i=0;i<3;++i) R[i][1][1] = D[i];


    for (i=0;i<3;++i) n[i][0][0] = nA[i];
    for (i=0;i<3;++i) n[i][1][0] = nB[i];
    for (i=0;i<3;++i) n[i][0][1] = nC[i];
    for (i=0;i<3;++i) n[i][1][1] = nD[i];

    xi[0]         = xiA + xiB;
    xi[1]         = xiC + xiD;
    oneoverxi[0]  = 1.0/xi[0];
    oneoverxi[1]  = 1.0/xi[1];
    oneoverxiall  = 1.0/(xi[0]+xi[1]);

    zi[0]    = xiA*xiB*oneoverxi[0];
    zi[1]    = xiC*xiD*oneoverxi[1];

    rho = (xi[0]*xi[1])/(xi[0]+xi[1]);
   
    for (i=0;i<3;++i) P[i][0] = (xiA*A[i] + xiB*B[i])*oneoverxi[0];
    for (i=0;i<3;++i) P[i][1] = (xiC*C[i] + xiD*D[i])*oneoverxi[1];
    for (i=0;i<3;++i) W[i]    = (xi[0]*P[i][0] + xi[1]*P[i][1])*oneoverxiall;

    T = 0.0; for (i=0; i<3; ++i) T += (P[i][0]-P[i][1])*(P[i][0]-P[i][1]); T *= rho;

    /* compute basis two center integrals (s_A|s_B) and (s_C|s_D)*/
    sumAB = 0.0; for (i=0; i<3; ++i) sumAB += (A[i]-B[i])*(A[i]-B[i]); sumAB *= zi[0];
    tmp = exp(-sumAB);
    c1 = (pi*oneoverxi[0]); c1 = sqrt(c1*c1*c1);
    sa_sb[0] = c1*tmp;

    sumCD = 0.0; for (i=0; i<3; ++i) sumCD += (C[i]-D[i])*(C[i]-D[i]); sumCD *= zi[1];
    tmp = exp(-sumCD);
    c1 = (pi*oneoverxi[1]); c1 = sqrt(c1*c1*c1);
    sa_sb[1] = c1*tmp;

    Tconst = 2.0*sqrt(rho/pi);
    Tconst *= sa_sb[0]*sa_sb[1];
    mset = 0; sasb_V_scsd[mset] = Tconst*igamma(mset,T);
    //printf("Tconst=%le  T=%le,  sasb_V_scsd=%le\n",Tconst,T,sasb_V_scsd[0]);
    //printf("igamma(%d, %lf)=%le\n",mset,T,igamma(mset,T));
    //printf("rho=%lf\n",rho);


    /* call recursive two center integral formula */
    value = eri_sub(0,n);
                   
    return value;
}


double eri_sub(int m, int n[3][2][2])
{
   int    i,xyz,ab,lr,rl,nsumall,done;
   double value,value2;

   nsumall = 0;
   for (lr=0;  lr<2;  ++lr)
   for (ab=0;  ab<2;  ++ab)
   for (xyz=0; xyz<3; ++xyz)
      nsumall += n[xyz][ab][lr];

   if (nsumall == 0)
   {
      while (mset<m)
      {
        ++mset;
        sasb_V_scsd[mset] = Tconst*igamma(mset,T);
      }
      value = sasb_V_scsd[m];
   }
   else
   {
      done = 0;
      xyz  = 0;
      ab   = 0;
      lr   = 0;
      rl   = 1;
      while (!done)
      {
         if (lr>1)
         {
            printf("eri_sub: Should not happen\n\n");
            return 0.0;
         }
         if (n[xyz][ab][lr] > 0)
         {

            n[xyz][ab][lr] -= 1; 
            value  = (P[xyz][lr] - R[xyz][ab][lr])*eri_sub(m,  n);
            value -= (P[xyz][lr] - W[xyz]        )*eri_sub(m+1,n);

            value2 = 0.0;
            for (i=0; i<2; ++i)
            {
               n[xyz][i][lr] -= 1;
               if (n[xyz][i][lr] >= 0)
                  value2 += (n[xyz][i][lr]+1)*(eri_sub(m,n)-rho*oneoverxi[lr]*eri_sub(m+1,n));
               n[xyz][i][lr] += 1;
            }
            value += 0.5*oneoverxi[lr]*value2;

            value2 = 0.0;
            for (i=0; i<2; ++i)
            {
               n[xyz][i][rl] -= 1;
               if (n[xyz][i][rl] >= 0)
                  value2 += (n[xyz][i][rl]+1)*(eri_sub(m+1,n));
               n[xyz][i][rl] += 1;
            }
            value += 0.5*oneoverxiall*value2;


            n[xyz][ab][lr] += 1;
            done = 1;
         }

         xyz = (xyz+1)%3;
         if (xyz==0) 
         { 
            ab = (ab+1)%2;
            if (ab==0) 
            {
               ++lr; 
               --rl;
            } 
         }

      }
   }

   return value;
}


