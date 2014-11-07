
#include 	<stdio.h>
#include	<math.h>

#include	"int_three_center.h"
#include	"int_overlap.h"
#include	"int_kinetic.h"
#include	"int_nuclear.h"
#include	"int_eri.h"
#include	"igamma.h"

#define LMAX 1

main()
{

   int i1,j1,k1;
   int i2,j2,k2;
   int i3,j3,k3;
   int i4,j4,k4;
   int n1[3],n2[3],n3[3],n4[3];
   double xi1,xi2,xi3,xi4;
   double r1[3],r2[3],r3[3],r4[3];
   double value,value2,value3,value4;
   double       kalue2,kalue3,kalue4;
   double v1,eri;

   xi1 = 0.62391373;
   xi2 = 0.62391373;
   xi3 = 0.1688554;
   xi4 = 0.1688554;
   r1[0] =  0.40; r1[1] =  0.10;   r1[2] =  0.20;
   r2[0] =  0.20; r2[1] =  0.30;   r2[2] =  0.40;
   r3[0] =  0.30; r3[1] =  0.40;   r3[2] =  0.10;
   r4[0] =  0.10; r4[1] =  0.20;   r4[2] =  0.30;

/*
   for (k1=0; k1<=LMAX; ++k1) for (j1=0; j1<=LMAX; ++j1) for (i1=0; i1<=LMAX; ++i1) 
   for (k2=0; k2<=LMAX; ++k2) for (j2=0; j2<=LMAX; ++j2) for (i2=0; i2<=LMAX; ++i2) 
   for (k3=0; k3<=LMAX; ++k3) for (j3=0; j3<=LMAX; ++j3) for (i3=0; i3<=LMAX; ++i3) 
*/
   i1 = 0; j1 = 1; k1 = 0;
   i2 = 0; j2 = 0; k2 = 1;
   i3 = 1; j3 = 0; k3 = 0;
   i4 = 0; j4 = 1; k4 = 0;
   {
      n1[0] = i1; n1[1] = j1; n1[2] = k1;
      n2[0] = i2; n2[1] = j2; n2[2] = k2;
      n3[0] = i3; n3[1] = j3; n3[2] = k3;
      n4[0] = i4; n4[1] = j4; n4[2] = k4;

      value = int_three_center(n1,xi1,r1, n2,xi2,r2, n3,xi3,r3);

      value2 = int_overlap(n1,xi1,r1, n2,xi2,r2);
      value3 = int_overlap(n1,xi1,r1, n3,xi3,r3);
      value4 = int_overlap(n2,xi2,r2, n3,xi3,r3);

      kalue2 = int_kinetic(n1,xi1,r1, n2,xi2,r2);
      kalue3 = int_kinetic(n1,xi1,r1, n3,xi3,r3);
      kalue4 = int_kinetic(n2,xi2,r2, n3,xi3,r3);

      v1 = int_nuclear(n1,xi1,r1, r2, n3,xi3,r3);

      eri = int_eri(n1,xi1,r1, n2,xi2,r2, n3,xi3,r3, n4,xi4,r4);

      printf("(a:[n1=(%d,%d,%d), xi1= %lf, r1=(%lf %lf %lf)]\n",n1[0],n1[1],n1[2],xi1,r1[0],r1[1],r1[2]);
      printf("(c:[n2=(%d,%d,%d), xi2= %lf, r2=(%lf %lf %lf)]\n",n2[0],n2[1],n2[2],xi2,r2[0],r2[1],r2[2]);
      printf("(b:[n3=(%d,%d,%d), xi3= %lf, r3=(%lf %lf %lf)]\n",n3[0],n3[1],n3[2],xi3,r3[0],r3[1],r3[2]);
      printf("(d:[n4=(%d,%d,%d), xi4= %lf, r4=(%lf %lf %lf)]\n",n4[0],n4[1],n4[2],xi4,r4[0],r4[1],r4[2]);
      printf("(a|c|b)   =%lf\n\n",value);
      printf("(a|c)     =%lf\n",  value2);
      printf("(a|b)     =%lf\n",  value3);
      printf("(c|b)     =%lf\n\n",value4);
      printf("(a|ke|c)  =%lf\n",  kalue2);
      printf("(a|ke|b)  =%lf\n",  kalue3);
      printf("(c|ke|b)  =%lf\n\n",kalue4);

      printf("(a|V(C)|b)=%lf\n\n",v1);
      
      printf("(ac||bd)   =%le\n\n",eri);
   }

}
