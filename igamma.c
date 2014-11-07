#include	<math.h>

#define	M_MAX	200
#define T_MAX	30.0
#define	EXPLIM	100.0
#define	ACCY	1.0e-14
#define	PI	3.1415926535897931159979635

#define	min(a,b)	((a<b) ? a : b)

static double Fj[M_MAX];


double igamma(int m, double T)
{
   int    j;
   double t2,term,sum,t2inv,expt;

   if (T>T_MAX)
   {
      t2inv = 0.5/T;
      expt  = exp(-min(EXPLIM,T));
      Fj[0] = sqrt(0.25*PI/T);
      for (j=1; j<=m; ++j) Fj[j] = ((2*j-1)*Fj[j-1] - expt)*t2inv;
   }
   else
   {
      t2 = 2.0*T;
      expt  = exp(-min(EXPLIM,T));
      term = 1.0/((double) (2*m+1));

      sum = term;
      j   = 1;
      while ((term >= ACCY)&& (j<1000))
      {
         term = term *t2/((double) (2*(m+j)+1));
         sum += term;
         ++j;
      }
      Fj[m] = expt*sum;

   }
   return Fj[m];
}
