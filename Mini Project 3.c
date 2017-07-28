#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define G 6.67*pow(10,-11)
#define M1 0.08*1.989*pow(10,30)
#define M2 0.85*5.9722*pow(10,24)

double dist(double r1x,double r2x,double r1y,double r2y)
{
    return(pow(pow((r1x-r2x),2)+pow((r1y-r2y),2),1.5));
}
double f1(double u1,double u2)
{
        return(-1*G*M2*(u1-u2));
}
double f2(double u1,double u2)
{
        return(-1*G*M1*(u2-u1));
}
void rk()
{
    double t,dt=60.0,tMAX=1.509*24*60*60*2,del,desdel,phi;
    double x1=0, x2= 1.9889*pow(10,8), y1=0, y2=-1.651*pow(10,9);
    double v1=0, v2=79348.8, u1=0, u2=9078.256;
    double i, u, epsilon=0.0000000001, r, r1, r2, r3;
    double x11,x12,x13,x21,x22,x23,y11,y12,y13,y21,y22,y23;
    double sx1,sx2,sy1,sy2,sv1,sv2,su1,su2,v11,v12,v13,v21,v22,v23,u11,u12,u13,u21,u22,u23;
    FILE *fptr;
    fptr=fopen("orbitdata.csv","w");
    fprintf (fptr,"t,x1,y1,x2,y2\n");

for(t=0;t<tMAX;t+=dt)
{
    for(i=0;i<2;i++)
    {
        r=dist(x1,x2,y1,y2);
        x11=x1+.5*dt*v1;
        v11=v1+.5*dt*f1(x1,x2)/r;
        y11=y1+.5*dt*u1;
        u11=u1+.5*dt*f1(y1,y2)/r;
        x21=x2+.5*dt*v2;
        v21=v2+.5*dt*f2(x1,x2)/r;
        y21=y2+.5*dt*u2;
        u21=u2+.5*dt*f2(y1,y2)/r;

        r1=dist(x11,x21,y11,y21);
        x12=x1+.5*dt*v11;
        v12=v1+.5*dt*f1(x11,x21)/r1;
        y12=y1+.5*dt*u11;
        u12=u1+.5*dt*f1(y11,y21)/r1;
        x22=x2+.5*dt*v21;
        v22=v2+.5*dt*f2(x11,x21)/r1;
        y22=y2+.5*dt*u21;
        u22=u2+.5*dt*f2(y11,y21)/r1;

        r2=dist(x12,x22,y12,y22);
        x13=x1+dt*v12;
        v13=v1+dt*f1(x12,x22)/r2;
        y13=y1+dt*u12;
        u13=u1+dt*f1(y12,y22)/r2;
        x23=x2+dt*v22;
        v23=v2+dt*f2(x12,x22)/r2;
        y23=y2+dt*u22;
        u23=u2+dt*f2(y12,x22)/r2;

        r3=dist(x13,x23,y13,y23);
        sx1=x1+dt/6*(v1+2*v11+2*v12+v13);
        sv1=v1+dt/6*(f1(x1,x2)/r+2*f1(x11,x21)/r1+2*f1(x12,x22)/r2+f1(x13,x23)/r3);
        sy1=y1+dt/6*(u1+2*u11+2*u12+u13);
        su1=u1+dt/6*(f1(y1,y2)/r+2*f1(y11,y21)/r1+2*f1(y12,y22)/r2+f1(y13,y23)/r3);
        sx2=x2+dt/6*(v2+2*v21+2*v22+v23);
        sv2=v2+dt/6*(f2(x1,x2)/r+2*f2(x11,x21)/r1+2*f2(x12,x22)/r2+f2(x13,x23)/r3);
        sy2=y2+dt/6*(u2+2*u21+2*u22+u23);
        su2=u2+dt/6*(f2(y1,y2)/r+2*f2(y11,y21)/r1+2*f2(y12,y22)/r2+f2(y13,y23)/r3);
    }
        printf("Time: %.lf\tdt: %.lf\tx1: %e\ty1: %e\tx2: %e\ty2: %e\n",t,dt,x1,y1,x2,y2);
        fprintf(fptr,"%lf, %lf, %lf, %lf, %lf\n",t,x1,y1,x2,y2);

        r=dist(x1,x2,y1,y2);
        x11=x1+.5*dt*v1;
        v11=v1+.5*dt*f1(x1,x2)/r;
        y11=y1+.5*dt*u1;
        u11=u1+.5*dt*f1(y1,y2)/r;
        x21=x2+.5*dt*v2;
        v21=v2+.5*dt*f2(x1,x2)/r;
        y21=y2+.5*dt*u2;
        u21=u2+.5*dt*f2(y1,y2)/r;

        r1=dist(x11,x21,y11,y21);
        x12=x1+.5*dt*v11;
        v12=v1+.5*dt*f1(x11,x21)/r1;
        y12=y1+.5*dt*u11;
        u12=u1+.5*dt*f1(y11,y21)/r1;
        x22=x2+.5*dt*v21;
        v22=v2+.5*dt*f2(x11,x21)/r1;
        y22=y2+.5*dt*u21;
        u22=u2+.5*dt*f2(y11,y21)/r1;

        r2=dist(x12,x22,y12,y22);
        x13=x1+dt*v12;
        v13=v1+dt*f1(x12,x22)/r2;
        y13=y1+dt*u12;
        u13=u1+dt*f1(y12,y22)/r2;
        x23=x2+dt*v22;
        v23=v2+dt*f2(x12,x22)/r2;
        y23=y2+dt*u22;
        u23=u2+dt*f2(y12,x22)/r2;

        r3=dist(x13,x23,y13,y23);
        x1=x1+dt/6*(v1+2*v11+2*v12+v13);
        v1=v1+dt/6*(f1(x1,x2)/r+2*f1(x11,x21)/r1+2*f1(x12,x22)/r2+f1(x13,x23)/r3);
        y1=y1+dt/6*(u1+2*u11+2*u12+u13);
        u1=u1+dt/6*(f1(y1,y2)/r+2*f1(y11,y21)/r1+2*f1(y12,y22)/r2+f1(y13,y23)/r3);
        x2=x2+dt/6*(v2+2*v21+2*v22+v23);
        v2=v2+dt/6*(f2(x1,x2)/r+2*f2(x11,x21)/r1+2*f2(x12,x22)/r2+f2(x13,x23)/r3);
        y2=y2+dt/6*(u2+2*u21+2*u22+u23);
        u2=u2+dt/6*(f2(y1,y2)/r+2*f2(y11,y21)/r1+2*f2(y12,y22)/r2+f2(y13,y23)/r3);

        if(v1>v2&&v1>u1&&v1>u2)
        {
            u=v1;
            desdel=epsilon*(fabs(u)+dt*fabs(f1(x1,x2)/r));
            del=fabs(sv1-u);
            phi=desdel/del;
            if (phi>1)
            dt*=pow(phi,0.2);
            if (phi<1)
            dt*=pow(phi,0.25);
        }
        else if (v2>v1&&v2>u1&&v2>u2)
        {
            u=v2;
            desdel=epsilon*(fabs(u)+dt*fabs(f2(x1,x2)/r));
            del=fabs(sv2-u);
            phi=desdel/del;
            if (phi>1)
            dt*=pow(phi,0.2);
            if (phi<1)
            dt*=pow(phi,0.25);
        }
        else if (u1>v2&&u1>v1&&u1>u2)
        {
            u=u1;
            desdel=epsilon*(fabs(u)+dt*fabs(f1(x1,x2)/r));
            del=fabs(su1-u);
            phi=desdel/del;
            if (phi>1)
            dt*=pow(phi,0.2);
            if (phi<1)
            dt*=pow(phi,0.25);
        }
        else if (u2>v2&&u2>u1&&u2>v1)
        {
            u=u2;
            desdel=epsilon*(fabs(u)+dt*fabs(f2(x1,x2)/r));
            del=fabs(su2-u);
            phi=desdel/del;
            if (phi>1)
            dt*=pow(phi,0.2);
            if (phi<1)
            dt*=pow(phi,0.25);
        }
}
fclose(fptr);
}

void main()
{
rk();
}
