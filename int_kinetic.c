
#include	<math.h>
#include	<stdio.h>
#include	"int_kinetic.h"

static double pi=3.1415926535897931159979635;
static double xi,zi,xiarray[2];
static double oneoverxi;
static double R[3][2],P[3];

/* two center integral (s_A|s_B) = (0_A|0_B) */
static double sa_sb;

/* two center integral (s_A|ke|s_B) = (0_A|ke|0_B) */
static double sa_ke_sb;


extern double ke_overlap_sub();
extern double ke_sub();

double int_kinetic(int nA[], double xiA, double A[], 
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

    xiarray[0] = xiA;
    xiarray[1] = xiB;
    xi            = xiA + xiB;
    oneoverxi     = 1.0/xi;

    zi     = xiA*xiB*oneoverxi;

   
    for (i=0;i<3;++i) P[i] = (xiA*A[i] + xiB*B[i])*oneoverxi;


    /* compute basis two center integral (s_A|s_B) */
    sumAB = 0.0; for (i=0; i<3; ++i) sumAB += (A[i]-B[i])*(A[i]-B[i]); sumAB *= zi;
    c1 = (pi*oneoverxi); c1 = sqrt(c1*c1*c1);
    sa_sb = c1*exp(-sumAB);

    /* compute basis two center integral (s_A|ke|s_B) */
    sa_ke_sb = zi*(3.0-2.0*sumAB)*sa_sb;



    /* call recursive two center integral formula */
    value = ke_sub(n);
                   
    return value;
}




double ke_sub(int n[3][2])
{
   int    i,xyz,ab,nsumall,done;
   double value,value2;

   nsumall = 0;
   for (ab=0;  ab<2;  ++ab)
   for (xyz=0; xyz<3; ++xyz)
      nsumall += n[xyz][ab];

   if (nsumall == 0)
   {
      value = sa_ke_sb;
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
            printf("ke_sub: Should not happen\n\n");
            return 0.0;
         }

         if (n[xyz][ab] > 0)
         {
                                                    // current state (a)(b)
            /* 1st overlap term 2*zi*(a|b) */
            value = 2*zi*ke_overlap_sub(n);

            /* 1st ke term  (Pi-Ai)*(a-1|ke|b) */
            n[xyz][ab] -= 1; 
                                                    // current state (a-1)(b)
            value += (P[xyz] - R[xyz][ab])*ke_sub(n);

            /* 2nd and third ke terms  1/(2*xi)*Nxyz(a-1)*(a-2|ke|b) + 1/(2*xi)*Nxyz(b)*(a-1|ke|b-1) */
            value2 = 0.0;
            for (i=0; i<2; ++i)
            {
               n[xyz][i] -= 1;
                                                    // current state (a-2)(b) or (a-1)(b-1)
               if (n[xyz][i] >= 0)
                  value2 += (n[xyz][i]+1)*ke_sub(n);
               n[xyz][i] += 1;
                                                    // current state (a-1)(b)
            }
            value += 0.5*oneoverxi*value2;

            /* 2nd overlap term  1/(2*xi_xyz)*Nxyz(a-1)*(a-2|b) */
            n[xyz][ab] -= 1;
                                                    // current state (a-2)(b)
            if (n[xyz][ab] >= 0)
               value -= (zi/xiarray[ab])*(n[xyz][ab]+1)*ke_overlap_sub(n);
            n[xyz][ab] += 1;
                                                    // current state (a-1)(b)


            n[xyz][ab] += 1;
                                                    // current state (a)(b)

            done = 1;
         }

         xyz = (xyz+1)%3;
         if (xyz==0) ++ab;
      }
   }

   return value;
}




double ke_overlap_sub(int n[3][2])
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
            printf("ke_overlap_sub:Should not happen\n\n");
            return 0.0;
         }

         if (n[xyz][ab] > 0)
         {
            n[xyz][ab] -= 1;
            value = (P[xyz] - R[xyz][ab])*ke_overlap_sub(n);

            value2 = 0.0;
            for (i=0; i<2; ++i)
            {
               n[xyz][i] -= 1;
               if (n[xyz][i] >= 0)
                  value2 += (n[xyz][i]+1)*ke_overlap_sub(n);
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

