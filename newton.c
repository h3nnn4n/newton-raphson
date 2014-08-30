#include <stdio.h>
#include <stdlib.h>

#include "complex.h"

#define ee 1e-10
#define MAX_ROOTS 10

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
void newthon_raphson(_complex,_roots*,int,int*,double*,int,int,double,int,int);
_complex f(_complex, int, double []);
_complex df(_complex, int, double[]);


int main(int argc,char *argv[]){
    int 	screenx,screeny,
            ix,iy,
            i,j,k,
            kit,itest,grad,nit;
    double	zx,zy,zxn,zyn,cx,cy,
            sens,
            minx,maxx,miny,maxy,
            xcenter,ycenter,size,
            dy,xin,yin,x,y,
            eps,x0,y0,dx,
            poly[MAX_ROOTS+1];
    char	answ;
    _roots roots;
    _complex z,w,fz,dfz,*f_space,*df_space,*z_space;

    root_init(&roots);

    FILE *img=fopen("mandel.ppm","wt");

#ifdef DEBB

    puts("Using debug files...");

    FILE *z_space_real=fopen("z_space_real.dat","wt");
    FILE *z_space_imag=fopen("z_space_imag.dat","wt");
    FILE *f_space_real=fopen("f_space_real.dat","wt");
    FILE *f_space_imag=fopen("f_space_imag.dat","wt");
    FILE *df_space_real=fopen("df_space_real.dat","wt");
    FILE *df_space_imag=fopen("df_space_imag.dat","wt");
#endif
    _color *bitmap,col;
    int *escapetime;

    screenx=1600;
    screeny=1200;

    grad=9;
    nit=100;

    x0=0;
    y0=0;
    dx=2;

    eps=1E-10;

    escapetime=(int*)malloc(sizeof(int)*screenx*screeny);
    z_space=(_complex*)malloc(sizeof(_complex)*screenx*screeny);
    f_space=(_complex*)malloc(sizeof(_complex)*screenx*screeny);
    df_space=(_complex*)malloc(sizeof(_complex)*screenx*screeny);

    fprintf(img,"P3\n%d %d\n255\n",screenx,screeny);

    for(i=0;i<=grad;i++){
        if(i%2==0){
            poly[i]=i+1;
        }else{
            poly[i]=-i+2;
        }
    }

    poly[9]=-4;
    poly[8]=12;
    poly[7]=-4;
    poly[6]=0;
    poly[5]=-12;
    poly[4]=1;
    poly[3]=3;
    poly[2]=-2;
    poly[1]=0;
    poly[0]=-7;

    printf("Coeficients:\n");
    for(i=grad;i>=0;i--){
        printf(" a%d=%+.2f\n",i,poly[i]);
    }


    printf("\n f(0+0i)=");
    complex_print(f(complex_init(0,0),grad, poly));
    printf("df(0+0i)=");
    complex_print(df(complex_init(0,0),grad, poly));

    roots.grad=grad;

    dy=(dx*screeny)/screenx;

    xin=x0-dx;
    yin=y0-dy;

    for(i=0;i<screeny;i++){
        y=yin+2*i*dy/screeny;
        for(j=0;j<screenx;j++){
            x=xin+2*j*dx/screenx;
            z=complex_init(x,y);
            itest=1;

            newthon_raphson(z,&roots,nit,escapetime,poly,screenx,screenx,eps,i,j);

            z_space[i*screenx+j]=z;
            f_space[i*screenx+j]=fz;
            df_space[i*screenx+j]=dfz;
        }
    }

    for(i=0;i<screeny;i++){
        for(j=0;j<screenx;j++){
            col.r=escapetime[i*screenx+j]-1;

            if(col.r==0){
                col.r=255;
                col.g=0;
                col.b=0;
            }else if(col.r==1){
                col.r=0;
                col.g=255;
                col.b=0;
            }else if(col.r==2){
                col.r=0;
                col.g=0;
                col.b=255;
            }else if(col.r==3){
                col.r=255;
                col.g=255;
                col.b=0;
            }else if(col.r==4){
                col.r=255;
                col.g=0;
                col.b=255;
            }else if(col.r==5){
                col.r=0;
                col.g=255;
                col.b=255;
            }else if(col.r==6){
                col.r=0;
                col.g=255;
                col.b=127;
            }else if(col.r==7){
                col.r=255;
                col.g=127;
                col.b=0;
            }else if(col.r==8){
                col.r=255;
                col.g=0;
                col.b=127;
            }else if(col.r==9){
                col.r=127;
                col.g=255;
                col.b=0;
            }else if(col.r==10){
                col.r=127;
                col.g=0;
                col.b=255;
            }else{
                col.r=0;
                col.g=0;
                col.b=0;
            }

#ifdef DEBB
            fprintf(z_space_real,"%lf ",complex_real(z_space[i*screenx+j]));
            fprintf(z_space_imag,"%lf ",complex_imag(z_space[i*screenx+j]));
            fprintf(f_space_real,"%lf ",complex_real(f_space[i*screenx+j]));
            fprintf(f_space_imag,"%lf ",complex_imag(f_space[i*screenx+j]));
            fprintf(df_space_real,"%lf ",complex_real(df_space[i*screenx+j]));
            fprintf(df_space_imag,"%lf ",complex_imag(df_space[i*screenx+j]));
#endif

            fprintf(img,"%d %d %d ",(int)col.r,
                                    (int)col.g,
                                    (int)col.b);
        }
        fputc('\n',img);
#ifdef DEBB
        fputc('\n',z_space_real);
        fputc('\n',z_space_imag);
        fputc('\n',f_space_real);
        fputc('\n',f_space_imag);
        fputc('\n',df_space_real);
        fputc('\n',df_space_imag);
#endif
    }

    printf("\nZeroes and their images:\n\n");
    for(i=0;i<grad;i++){
        printf("Root Z%d=%+lf %+lfi \tf(z%d)=",i,roots.root[i].x,roots.root[i].y,i);
        complex_print(f(roots.root[i],grad, poly));
    }

    fclose(img);
    free(escapetime);
    free(z_space);
    free(f_space);
    free(df_space);

    return EXIT_SUCCESS;
}

void newthon_raphson(_complex z,_roots* roots,int nit,int *escapetime,double poly[],int screenx, int screeny,double eps,int i,int j){
    int k,kit,itest;
    _complex w,fz,dfz;
    itest=1;
    for(k=0;k<=nit;k++){
        kit=k;
        fz = f(z,roots->grad,poly);
        dfz=df(z,roots->grad,poly);
        if(complex_abs(dfz)<ee){
            itest=0;
            break;
        }
        w=z;
        z=complex_sub(z,complex_div(fz,dfz));
        if(complex_abs(complex_sub(z,w))<=eps){
            break;
        }
    }
    if(kit==nit){
        itest=2;
    }
    if(i<screeny && i>=0 && j<screenx && j>=0){
        if(itest!=1){
                escapetime[i*screenx+j]=666;
        }else{
            escapetime[i*screenx+j]=process_root(z,roots,eps);
        }
    }
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
