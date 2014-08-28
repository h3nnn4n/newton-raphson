#include <stdio.h>
#include <stdlib.h>

#include "complex.h"

#define ee 1e-60
#define MAX_ROOTS 50

typedef struct{
    _complex root[MAX_ROOTS];
    int nor;
    int grad;
} _roots;

typedef struct{
	double r,g,b;
}_color;

int process_root(_complex,_roots*,double);
void root_init(_roots *);
_complex f(_complex, int, double []);
_complex df(_complex, int, double[]);

int main(){
    int     ix,iy,radius,i,tries,
            grad,nit;
    double	zx,zy,zxn,zyn,cx,cy,theta,
            x,y,eps,
            poly[MAX_ROOTS+1];
    _roots roots;
    _complex z,w,fz,dfz;

    root_init(&roots);

    nit=100;

    radius=6;

    eps=1E-10;

    scanf("%d",&grad);

    roots.grad=grad;

    for(i=grad;i>=0;i--){
        scanf("%lf",&poly[i]);
    }

    printf("Coefficients:\n");
    for(i=grad;i>=0;i--){
        printf(" a%d=%+.2f\n",i,poly[i]);
    }

    printf("\n f(0+0i)=");
    complex_print(f(complex_init(0,0),grad, poly));
    printf("df(0+0i)=");
    complex_print(df(complex_init(0,0),grad, poly));

    tries=0;

    do{
        tries++;
        theta=drand48()*2*M_PI;

        x=radius*cos(theta);
        y=radius*sin(theta);

        z=complex_init(x,y);
        for(i=0;i<=nit;i++){
            fz = f(z,roots.grad,poly);
            dfz=df(z,roots.grad,poly);
            if(complex_abs(dfz)<ee){
                break;
            }
            w=z;
            z=complex_sub(z,complex_div(fz,dfz));
            if(complex_abs(complex_sub(z,w))<=eps){
                process_root(z,&roots,eps);
                break;
            }
        }
    }while(roots.nor<grad);

    printf("\nTook %d tries to get all %d roots\n",tries,grad);

    printf("\nZeroes and their images:\n\n");
    for(i=0;i<grad;i++){
        printf("Root Z%d=%+lf %+lfi \tf(z%d)=",i+1,roots.root[i].x,roots.root[i].y,i+1);
        complex_print(f(roots.root[i],grad, poly));
    }

    return EXIT_SUCCESS;
}

int process_root(_complex z, _roots *p, double eps){
    int i;
    _color c;
    for(i=0;i<p->nor;i++){
        if(complex_abs(complex_sub(z,p->root[i]))<2*eps){
            return i+1;
        }
    }

    p->root[p->nor]=z;
    p->nor+=1;

    return p->nor;
}

void root_init(_roots *roots){
    int i;
    for(i=0;i<MAX_ROOTS;i++){
        roots->root[i]=complex_init(0,0);
    }
    roots->nor=0;
    roots->grad=0;
}

_complex f(_complex z, int grad, double poly[]){
    int i;
    _complex f;
    f=complex_init(poly[grad],0);
    for(i=grad-1;i>=0;i--){
        f=complex_sum(complex_mult(f,z),complex_init(poly[i],0));
    }
    return f;
}

_complex df(_complex z, int grad, double poly[]){
    int i;
    _complex df;
    df=complex_init(grad*poly[grad],0);
    for(i=grad-1;i>0;i--){
        df=complex_sum(complex_mult(df,z),complex_init(i*poly[i],0));
    }
    return df;
}
