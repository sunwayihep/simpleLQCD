/************************************************************
 *
 * Luscher's random number generator
 *
 ***********************************************************/
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include "random.h"


#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static double one_bit;
static vec_t carry;

static union
{
   dble_vec_t vec[12];
   int num[96];
} x;

#define STEP(pi,pj) \
      d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
      (*pi).c2.c1+=(d<0); \
      d+=BASE; \
      (*pi).c1.c1=d&MASK; \
      d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
      (*pi).c2.c2+=(d<0); \
      d+=BASE; \
      (*pi).c1.c2=d&MASK; \
      d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
      (*pi).c2.c3+=(d<0); \
      d+=BASE; \
      (*pi).c1.c3=d&MASK; \
      d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
      (*pi).c2.c4+=(d<0); \
      d+=BASE; \
      (*pi).c1.c4=d&MASK; \
      d=(*pj).c2.c1-(*pi).c2.c1; \
      carry.c1=(d<0); \
      d+=BASE; \
      (*pi).c2.c1=d&MASK; \
      d=(*pj).c2.c2-(*pi).c2.c2; \
      carry.c2=(d<0); \
      d+=BASE; \
      (*pi).c2.c2=d&MASK; \
      d=(*pj).c2.c3-(*pi).c2.c3; \
      carry.c3=(d<0); \
      d+=BASE; \
      (*pi).c2.c3=d&MASK; \
      d=(*pj).c2.c4-(*pi).c2.c4; \
      carry.c4=(d<0); \
      d+=BASE; \
      (*pi).c2.c4=d&MASK

static int err_no, err_flg=0;
static char prog_name[128], err_msg[512];

#define MAX_TAG 32767
#define MAX_PERMANENT_TAG MAX_TAG/2

static int error_loc(int test,int no,char *name,char *message)
{
   if ((test!=0)&&(err_flg==0))
   {
      err_no=no;
      strncpy(prog_name,name,128);
      strncpy(err_msg,message,512);
      prog_name[127]='\0';
      err_msg[511]='\0';
      err_flg=1;
   }

   return test;
}

static void update(void)
{
   int k,kmax,d;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];

   for (k=0;k<kmax;k++)
   {
      STEP(pi,pj);
      pi+=1;
      pj+=1;
      if (pi==pmax)
         pi=pmin;
      if (pj==pmax)
         pj=pmin;
   }

   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


static void define_constants(void)
{
   int k;

   one_bit=ldexp(1.0,-24);

   for (k=0;k<96;k++)
   {
      next[k]=(k+1)%96;
      if ((k%4)==3)
         next[k]=(k+5)%96;
   }
}


void rlx_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   error_loc((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||
             (DBL_MANT_DIG<48),1,"rlxd_init [ranlxd.c]",
             "Arithmetic on this machine is not suitable for ranlxd");

   define_constants();

   error_loc((level<1)||(level>2),1,"rlxd_init [ranlxd.c]",
             "Bad choice of luxury level (should be 1 or 2)");

   if (level==1)
      pr=202;
   else if (level==2)
      pr=397;

   i=seed;

   for (k=0;k<31;k++)
   {
      xbit[k]=i%2;
      i/=2;
   }

   error_loc((seed<=0)||(i!=0),1,"rlxd_init [ranlxd.c]",
             "Bad choice of seed (should be between 1 and 2^31-1)");

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++)
         {
            iy=xbit[ibit];
            ix=2*ix+iy;

            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)!=i)
            ix=16777215-ix;

         x.num[4*k+i]=ix;
      }
   }

   carry.c1=0;
   carry.c2=0;
   carry.c3=0;
   carry.c4=0;

   ir=0;
   jr=7;
   is=91;
   is_old=0;
   prm=pr%12;
   init=1;
}


void ranlux(Float r[],int n)
{
   int k;

   if (init==0)
      rlx_init(1,1);

   for (k=0;k<n;k++)
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=one_bit*((Float)(x.num[is+4])+one_bit*(Float)(x.num[is]));
   }
}


static double twopi;

static int init2 =0;

void gauss(Float rd[],int n)
{
   int k;
   Float ud[2];
   Float x1,x2,rho,y1,y2;

   if (init2==0)
   {
      twopi=8.0*atan(1.0);
      init=1;
   }

   for (k=0;k<n;)
   {
      ranlux(ud,2);
      x1=ud[0];
      x2=ud[1];

      rho=-log(1.0-x1);
      rho=sqrt(rho);
      x2*=twopi;
      y1=rho*sin(x2);
      y2=rho*cos(x2);

      rd[k++]=y1;
      if (k<n)
         rd[k++]=y2;
   }
}

static void normalize(su3_vector *v)
{
   int i;
   Float *r,fact;

   r=(Float*)(v);
   fact=0.0f;

   for (i=0;i<6;i++)
      fact+=r[i]*r[i];

   fact=1.0f/(Float)sqrt((double)(fact));

   for (i=0;i<6;i++)
      r[i]*=fact;
}

/*
* v.c1=(w.c2*z.c3-w.c3*z.c2)^*
* v.c2=(w.c3*z.c1-w.c1*z.c3)^*
* v.c3=(w.c1*z.c2-w.c2*z.c1)^*
*/

#define _vector_cross_prod(v,w,z) \
   (v).c1.re= (w).c2.re*(z).c3.re-(w).c2.im*(z).c3.im  \
             -(w).c3.re*(z).c2.re+(w).c3.im*(z).c2.im; \
   (v).c1.im= (w).c3.re*(z).c2.im+(w).c3.im*(z).c2.re  \
             -(w).c2.re*(z).c3.im-(w).c2.im*(z).c3.re; \
   (v).c2.re= (w).c3.re*(z).c1.re-(w).c3.im*(z).c1.im  \
             -(w).c1.re*(z).c3.re+(w).c1.im*(z).c3.im; \
   (v).c2.im= (w).c1.re*(z).c3.im+(w).c1.im*(z).c3.re  \
             -(w).c3.re*(z).c1.im-(w).c3.im*(z).c1.re; \
   (v).c3.re= (w).c1.re*(z).c2.re-(w).c1.im*(z).c2.im  \
             -(w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \
   (v).c3.im= (w).c2.re*(z).c1.im+(w).c2.im*(z).c1.re  \
             -(w).c1.re*(z).c2.im-(w).c1.im*(z).c2.re

static void cross_prod(su3_vector *v1,su3_vector *v2,su3_vector *v3)
{
   _vector_cross_prod(*v3,*v1,*v2);
}

 void random_su3_vector(su3_vector *v)
{
   int i;
   Float *r,norm,fact;

   r=(Float*)(v);
   norm=0.0f;

   while (norm<=FLT_EPSILON)
   {
      gauss(r,6);
      norm=0.0f;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=(Float)sqrt((double)norm);
   }

   fact=1.0f/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;
}

void random_su3(su3 *u)
{
   int i;
   Float *r,norm,fact;
   su3_vector *v1,*v2,*v3;

   v1=(su3_vector*)(u);
   v2=v1+1;
   v3=v1+2;
   r=(Float*)(v3);

   random_su3_vector(v1);
   norm=0.0f;

   while (norm<=FLT_EPSILON)
   {
      random_su3_vector(v2);
      cross_prod(v1,v2,v3);
      norm=0.0f;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=(Float)sqrt((double)norm);
   }        

   fact=1.0f/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;   

   cross_prod(v3,v1,v2);
}

void random_spinor(spinor* ss, Float sigma)
{
  Float r[24], *s;

  gauss(r, 24);
  s=(Float*)(ss);

  for(int i=0; i<24; i++)
    s[i] = sigma*r[i];
}




