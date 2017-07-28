#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.67*pow(10, -11)
#define M1 .08*1.989*pow(10, 30)
#define M2 .85*5.9722*pow(10, 24)

double f1(double r1, double r2);
double f2(double r1, double r2);
double r12(double x1, double x2, double y1, double y2);
void rungekutta();

int main()
{
    rungekutta();

    return 0;
}

void rungekutta()
{
    double t, dt=60.0, tMax=1.509*24*60*60;
    double x1=0, x2= 1.9889*pow(10,8), y1=0, y2=-1.651*pow(10,9);
    double v1=0, v2=79348.8, u1=0, u2=9078.256;
    double r, r1, r2, r3;
    double x11,x12,x13,x21,x22,x23,y11,y12,y13,y21,y22,y23;
    double v11,v12,v13,v21,v22,v23,u11,u12,u13,u21,u22,u23;
    FILE *fptr = fopen("Orbit.csv", "w");
    fprintf(fptr,"Time, X1, Y1, X2, Y2\n");
    for(t=0; t<tMax; t+=dt)
    {
        printf("Time: %.lf\tx1: %e\ty1: %e\tx2: %e\ty2: %e\n",t,x1,y1,x2,y2);
        fprintf(fptr,"%lf, %lf, %lf, %lf, %lf\n",t,x1,y1,x2,y2);
        r=r12(x1,x2,y1,y2);
        x11=x1+.5*dt*v1;
        v11=v1+.5*dt*f1(x1,x2)/r;
        y11=y1+.5*dt*u1;
        u11=u1+.5*dt*f1(y1,y2)/r;
        x21=x2+.5*dt*v2;
        v21=v2+.5*dt*f2(x1,x2)/r;
        y21=y2+.5*dt*u2;
        u21=u2+.5*dt*f2(y1,y2)/r;

        r1=r12(x11,x21,y11,y21);
        x12=x1+.5*dt*v11;
        v12=v1+.5*dt*f1(x11,x21)/r1;
        y12=y1+.5*dt*u11;
        u12=u1+.5*dt*f1(y11,y21)/r1;
        x22=x2+.5*dt*v21;
        v22=v2+.5*dt*f2(x11,x21)/r1;
        y22=y2+.5*dt*u21;
        u22=u2+.5*dt*f2(y11,y21)/r1;

        r2=r12(x12,x22,y12,y22);
        x13=x1+dt*v12;
        v13=v1+dt*f1(x12,x22)/r2;
        y13=y1+dt*u12;
        u13=u1+dt*f1(y12,y22)/r2;
        x23=x2+dt*v22;
        v23=v2+dt*f2(x12,x22)/r2;
        y23=y2+dt*u22;
        u23=u2+dt*f2(y12,x22)/r2;

        r3=r12(x13,x23,y13,y23);
        x1=x1+dt/6*(v1+2*v11+2*v12+v13);
        v1=v1+dt/6*(f1(x1,x2)/r+2*f1(x11,x21)/r1+2*f1(x12,x22)/r2+f1(x13,x23)/r3);
        y1=y1+dt/6*(u1+2*u11+2*u12+u13);
        u1=u1+dt/6*(f1(y1,y2)/r+2*f1(y11,y21)/r1+2*f1(y12,y22)/r2+f1(y13,y23)/r3);
        x2=x2+dt/6*(v2+2*v21+2*v22+v23);
        v2=v2+dt/6*(f2(x1,x2)/r+2*f2(x11,x21)/r1+2*f2(x12,x22)/r2+f2(x13,x23)/r3);
        y2=y2+dt/6*(u2+2*u21+2*u22+u23);
        u2=u2+dt/6*(f2(y1,y2)/r+2*f2(y11,y21)/r1+2*f2(y12,y22)/r2+f2(y13,y23)/r3);
    }
    fclose(fptr);
}

double f1(double r1, double r2)
{
    double a = -1*G*M2*(r1-r2);
    return a;
}

double f2(double r1, double r2)
{
    double a = -1*G*M1*(r2-r1);
    return a;
}

double r12(double x1, double x2, double y1, double y2)
{
    double r = pow((pow((x1-x2),2) + pow((y1-y2),2)),1.5);
    return r;
}
