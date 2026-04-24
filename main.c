#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"fftw3.h"

double
*array1(int nx)
{
  double *m;
    m=(double *) fftw_malloc( (size_t) nx*sizeof(double ) );  
    return m;
}

double
**array2(int nx, int ny)
{
    int i;
    double **m;

    m=(double **) fftw_malloc( (size_t) nx*sizeof(double *) );
    m[0]=(double *) fftw_malloc( (size_t) nx*ny*sizeof(double) );

    for (i=1; i<nx; i++)
    {
      m[i]=m[i-1]+ny;
    }
    
    return m;
}

int main(int argc, char *argv[])
{
    int m=37,n=1000,i,j,tem[8];
    double *lambda, *segmaexp, *p, *pnew, *deltap, *b, *miu, **theta,segma0[8];
    double **ki;
    long int itt=100000000,it;
    double alpha,erro,z,sum,pmax1,miumax1,n1,pmax2,miumax2,n2;
    FILE *fp1, *fp2,*fp3,*fp4;
    int a=atof(argv[1]);

    char name1[500];
    char name2[500];
    
    alpha=0.2;

    b=array1(m);
    miu=array1(n);

    lambda=array1(2*m);
    deltap=array1(n);
    segmaexp=array1(2*m);
    p=array1(n);
    pnew=array1(n);
    ki=array2(n,2*m);
    theta=array2(n,m);
    
    fp1=fopen("segma0.txt","r");
    for(j=0; j<8; j++)
    {
      fscanf(fp1,"%d  %lf\n",&tem[j], &segma0[j]);
      //printf("%d  %lf\n", tem[j], segma0[j]);
    }
    fclose(fp1);

    sprintf(name1,"%dK.txt", a);
    fp2=fopen(name1,"r");
    for(j=0; j<m; j++)
    {
      fscanf(fp2,"%lf  %lf  %lf\n",&b[j], &segmaexp[m+j], &segmaexp[j]);
      //printf("%lf  %lf  %lf\n", b[j], segmaexp[m+j], segmaexp[j]);
    }
    fclose(fp2);


    for(i=0; i<n; i++)
    {
      miu[i]=3*i/n-1.5;   //miu:-10~10
      //printf("%f\n", miu[i]);
    }

    for(j=0; j<2*m; j++)
    {
      lambda[j]=0.0;
      //printf("%f ", lambda[j]);
    }
    //printf("\n");
    //printf("****************\n");
    
    for(i=0; i<n; i++)
    {
      for(j=0; j<m; j++)
      {
        theta[i][j]=atan(miu[i]*b[j]);
        //printf("%f ", theta[i][j]);
      }
      //printf("\n");
    }

    for(i=0; i<n; i++)
    {
      for(j=0; j<2*m; j++)
      {
        if(j<m)
        {
          ki[i][j]=cos(theta[i][j])*cos(theta[i][j]);
        }
        else
        {
          ki[i][j]=sin(2.0*theta[i][j-m])/2.0;
        }
        //printf("%f ", ki[i][j]);
      }
      //printf("\n");
    }

    for (it=0; it<itt; it++)
    {
      //printf("****************************\n");
      z=0.0;
      for (i=0; i<n; i++)
      {
        sum=0.0;
        for(j=0; j<m; j++)
        {
          sum += lambda[j]*ki[i][j]+lambda[m+j]*ki[i][m+j];
        }
        p[i]=exp(-sum);
        //printf("%f\n", p[i]);
        z += exp(-sum);
      }
      //printf("***************************\n");
      sum=0.0;
      for(i=0; i<n; i++)
      {
        p[i] /= z;
        sum += p[i];
      }
      //printf("sum_pi_%f\n", sum);

      for(j=0;j<2*m;j++)
      {
        sum=0.0;
        for(i=0;i<n;i++)
        {
          sum += ki[i][j]*p[i];
        }
        lambda[j] += -alpha*(segmaexp[j]-sum);
      }

      sum=0.0;
      for(i=0; i<n; i++)
      {
        pnew[i]=0;
        for(j=0; j<m; j++)
        {
          pnew[i] += ki[i][j]*lambda[j]+ki[i][j+m]*lambda[j+m];
        }
        pnew[i] = exp(-pnew[i]);
        sum += pnew[i];
      }
      for(i=0; i<n; i++)
      {
        pnew[i] /= sum;
      }

      erro=0.0;
      for(i=0; i<n; i++)
      {
        deltap[i]=fabs(p[i]-pnew[i]);
        if(erro<=deltap[i]) erro=deltap[i];
      }
      

      if(erro<0.00000001) break;
    }

    pmax1=0.0;
    miumax1=0.0;
    for(i=0;i<n/2-2;i++)
    {
      if((p[i+1]-p[i])/(miu[i+1]-miu[i])>=0 && (p[i+2]-p[i+1])/(miu[i+2]-miu[i+1])<=0)
      {
        if(pmax1<=p[i+1])
        {
          pmax1=p[i+1];
          miumax1=miu[i+1];
          printf("%f  %f\n", pmax1, miumax1);
        }
      }
    }

    pmax2=0.0;
    miumax2=0.0;
    for(i=n/2-2;i<n-2;i++)
    {
      if((p[i+1]-p[i])/(miu[i+1]-miu[i])>=0 && (p[i+2]-p[i+1])/(miu[i+2]-miu[i+1])<=0)
      {
        if(pmax2<=p[i+1])
        {
          pmax2=p[i+1];
          miumax2=miu[i+1];
          printf("%f  %f\n", pmax2, miumax2);
        }
      }
    }

    for(j=0;j<8;j++)
    {
      if(tem[j]==a)
      {
        n1=0.0;
        for(i=0;i<n/2-9;i++)
        {
          n1+=fabs((miu[i+1]-miu[i]))*p[i+1]/fabs(miu[i+1]);
        }
        n1*=segma0[j];

        n2=0.0;
        for(i=n/2+10;i<n-1;i++)
        {
          n2+=fabs((miu[i+1]-miu[i]))*p[i+1]/fabs(miu[i+1]);
        }
        n2*=segma0[j];

      }
    }

    fp3=fopen("n_miumax_pmax.txt", "a");
    fprintf(fp3, "%d  %f  %f  %f  %f  %f  %f\n", a, pmax1, miumax1, n1, pmax2, miumax2, n2);
    fclose(fp3);    

    sprintf(name2,"miu_p_%dK.txt", a);
    fp4=fopen(name2, "w");
    for(i=0; i<n; i++)
    {
      for(j=0;j<8;j++)
      {
        if(tem[j]==a)
        {
          fprintf(fp4, "%f  %f\n", miu[i], p[i]*segma0[j]/fabs(miu[i]));
        }
      }
    }
    fclose(fp4);


    fftw_free(lambda);
    fftw_free(segmaexp);
    fftw_free(p);
    fftw_free(pnew);
    fftw_free(deltap);
    fftw_free(b);
    fftw_free(miu);
    fftw_free(ki[0]);
    fftw_free(ki);
    fftw_free(theta[0]);
    fftw_free(theta);

    return 1;
}
