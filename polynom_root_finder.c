#include <stdio.h>
#include <stdlib.h>

#include "complex.h"

#define ee 1e-60
#define MAX_ROOTS 50

typedef struct{                 // This struct will contain at the end of the program execution
    _complex root[MAX_ROOTS];   // all roots.
    int nor;    // Number of Roots found;
    int grad;   // Degree
} _roots;

// Prototype of the functions
int process_root(_complex,_roots*,double);
void root_init(_roots *);
_complex f(_complex, int, double []);
_complex df(_complex, int, double[]);

int main(){
    int     ix,iy,radius,i,tries,nit;
    double	zx,zy,zxn,zyn,cx,cy,theta,
            x,y,eps,
            poly[MAX_ROOTS+1];
    _roots roots;
    _complex z,w,fz,dfz;

    root_init(&roots);

    // Number of iteration for each pointer.
    // The bigger, the less prone to error the program will be.
    nit=100;

    // The radius where the points will be looked for.
    radius=6;

    // Precision for the roots
    eps=1E-10;

    scanf("%d",&roots.grad);

    for(i=roots.grad;i>=0;i--){
        scanf("%lf",&poly[i]);
    }

    printf("Coefficients:\n");
    for(i=roots.grad;i>=0;i--){
        printf(" a%d=%+.2f\n",i,poly[i]);
    }

    printf("\n f(0+0i)=");
    complex_print(f(complex_init(0,0),roots.grad, poly));
    printf("df(0+0i)=");
    complex_print(df(complex_init(0,0),roots.grad, poly));

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
    }while(roots.nor<roots.grad);

    printf("\nTook %d tries to get all %d roots\n",tries,roots.grad);

    printf("\nZeroes and their images:\n\n");
    for(i=0;i<roots.grad;i++){
        printf("Root Z%d=%+lf %+lfi \tf(z%d)=",i+1,roots.root[i].x,roots.root[i].y,i+1);
            complex_print(f(roots.root[i],roots.grad, poly));
    }

    return EXIT_SUCCESS;
}

// When a root is found this function is called.
// It looks in a vector to find if the given root
// was already found, if yes it returns the position
// in the vector of the root, if not it simple adds
// the root in the first avaliable position.
int process_root(_complex z, _roots *p, double eps){
    int i;
    for(i=0;i<p->nor;i++){
        if(complex_abs(complex_sub(z,p->root[i]))<2*eps){
            return i+1;
        }
    }

    p->root[p->nor]=z;
    p->nor+=1;

    return p->nor;
}

// This function initialize the data structure for the roots
void root_init(_roots *roots){
    int i;
    for(i=0;i<MAX_ROOTS;i++){
        roots->root[i]=complex_init(0,0);
    }
    roots->nor=0;
    roots->grad=0;
}

// Here the polynomial function is evaluated at a point Z
// in the complex plane.
_complex f(_complex z, int grad, double poly[]){
    int i;
    _complex f;
    f=complex_init(poly[grad],0);
    for(i=grad-1;i>=0;i--){
        f=complex_sum(complex_mult(f,z),complex_init(poly[i],0));
    }
    return f;
}

// The evaluation of the derivative happens here in
// the same way as the funtion bove.
_complex df(_complex z, int grad, double poly[]){
    int i;
    _complex df;
    df=complex_init(grad*poly[grad],0);
    for(i=grad-1;i>0;i--){
        df=complex_sum(complex_mult(df,z),complex_init(i*poly[i],0));
    }
    return df;
}
