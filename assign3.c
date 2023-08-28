#include<stdio.h>
#include<math.h>
#include<conio.h>

int main()
{
    FILE *fp;
    FILE *P;
    FILE *Q;
    fp=fopen("result.dat","w");
    P=fopen("central_u.dat","w");
     Q=fopen("central_v.dat","w");
    int i,j,m,n,p,re,iteration;
    float dx,dy,x,y,b;
    printf("enter re= ");
    scanf("%d",&re);
     printf("\nenter m= ");
    scanf("%d",&m);
     printf("\nenter n= ");
    scanf("%d",&n);

    p=(m-2)*(n-2);            //total no. of interior points
    dx=1.0/(m-1);
    dy=1.0/(n-1);

    b=(dx/dy);               //beta

    float s[m][n],s_pr[m][n],w[m][n],w_pr[m][n],u[m][n],v[m][n];

    float s_error=0.0,w_error=0.0;

    iteration=0;

    //boundary condition initialization
    for(j=0;j<n;j++)
     {
         for(i=0;i<m;i++)
         {

             v[i][j]=0.0;
             s[i][j]=0.0;
             if(j==(n-1)){
              u[i][j]=1.0;
             }else
               u[i][j]=0.0;
         }
     }

     for(j=0;j<n;j++)
     {
          for(i=0;i<m;i++)
          {
              if(j==0)
              w[i][j]=(2.0/dy*dy)*(s[i][j]-s[i][j+1]);
              else if(i==0)
                w[i][j]=(2.0/dx*dx)*(s[i][j]-s[i+1][j]);
                else if(i==(m-1))
                w[i][j]=(2.0/dx*dx)*(s[i][j]-s[i-1][j]);
                 else if(j==(n-1))
                w[i][j]=(2.0/dy*dy)*(s[i][j]-s[i][j-1])-((2.0/dy)*u[i][j]);
                else
                    w[i][j]=0.0;

          }
     }

     //GAUSS-SIEDEL METHOD

     do
     {
         for(j=0;j<n;j++)
         {
             for(i=0;i<m;i++)
             {
                 s_pr[i][j]=s[i][j];
                 w_pr[i][j]=w[i][j];
             }
         }
         //SOLVING STREAM FUNCTION
          for(j=1;j<(n-1);j++)
         {
             for(i=1;i<(m-1);i++)
             {
                 s[i][j]=(1.0/(2*(1.0+b*b)))*(s[i+1][j]+s[i-1][j]+b*b*(s[i][j+1]+s[i][j-1])+dx*dx*w[i][j]);
             }
         }

         //SOLVING FOR VORTICITY
         for(j=1;j<(n-1);j++)
         {
             for(i=1;i<(m-1);i++)
             {
                 w[i][j] = (0.5/(1.0+pow(b,2)))*((1.0-(s[i][j+1]-s[i][j-1])*((b*re)/4.0))*w[i+1][j]+ (1.0+(s[i][j+1]-s[i][j-1])*((b*re)/4.0))*w[i-1][j]
                     + (1.0+(s[i+1][j]-s[i-1][j])*(re/(4.0*b)))*(pow(b,2)*w[i][j+1])+ (1.0-(s[i+1][j]-s[i-1][j])*(re/(4.0*b)))*(pow(b,2)*w[i][j-1]));;
             }
         }
         // UPDATE VORTICITY AT BOUNDARIES

          for(j=0; j<n; j++)
       {
           w[0][j]=(-2*(s[1][j]-s[0][j]))/(dx*dx);
           w[m-1][j]=(-2*(s[m-2][j]-s[m-1][j]))/(dx*dx);
       }
    for(i=0; i<m; i++)
       {
            w[i][0]=(-2*(s[i][1]-s[i][0]))/(dy*dy);
            w[i][n-1]=(-2*(s[i][n-2]-s[i][n-1]+dy))/(dy*dy);
       }
            //ERROR CALCULATION
            s_error=0.0;
            w_error=0.0;
            for(j=1;j<(n-1);j++)
         {
             for(i=1;i<(m-1);i++)
        {

              s_error=s_error+pow((s[i][j]-s_pr[i][j]),2.0);
              w_error=w_error+pow((w[i][j]-w_pr[i][j]),2.0);
          }
      }
      s_error=sqrt(s_error/p);
      w_error=sqrt(w_error/p);



      printf("iteration=%d\t",iteration);
      printf("error_psi=%.9lf\terror_omega=%.9lf\n",s_error,w_error);
      iteration++;
    }while(s_error>pow(10,-6) || w_error>pow(10,-6));

     //UPDATING VELOCITIES

         for (j=1;j<(n-1);j++)
         {
             for(i=1;i<(m-1);i++)
     {
             u[i][j]=(s[i][j+1]-s[i][j-1])/(2*dy);
             v[i][j]=(s[i+1][j]-s[i-1][j])/(-2.0*dx);
         }
     }
    fprintf(fp,"ZONE I=%d, J=%d\n",m,n);
for(i = 0; i < m; i++)
{
   x=i*dx;
for(j = 0; j < n; j++)
{
   y=j*dy;
fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\n",x,y,u[i][j],v[i][j],s[i][j],w[i][j]);
}
}
for(i=0;i<m;i++)
fprintf(P,"%f \t %f \n",u[n/2][i],i*dy);
for(j=0;j<n;j++)
fprintf(Q,"%lf \t %lf \n",j*dx,v[j][m/2]);
fclose(fp);
fclose(P);
fclose(Q);
return 0;
}
