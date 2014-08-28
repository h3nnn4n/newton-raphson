#include <stdio.h>
#include <math.h>

    typedef struct{
        double x,y;
    } _complex;

    _complex complex_init(double x, double y){
        _complex z;
        z.x=x;
        z.y=y;
        return z;
    }

    double complex_real(_complex z){
        return z.x;
    }

    double complex_imag(_complex z){
        return z.y;
    }

    _complex complex_opposite(_complex z){
        _complex z1;
        z1.x= (-1)*z.x;
        z1.y= (-1)*z.y;
    return z1;
}

_complex complex_conjugate(_complex z){
    _complex z1;
    z1.x=z.x;
    z1.y= (-1)*z.y;
    return z1;
}

double complex_abs(_complex z){
    return sqrt(pow(z.x,2)+pow(z.y,2));
}

_complex complex_sum(_complex z1, _complex z2){
    _complex z;
    z.x= z1.x + z2.x;
    z.y= z1.y + z2.y;
    return z;
}

_complex complex_sub(_complex z1,_complex z2){
    _complex z;
    z.x= z1.x - z2.x;
    z.y= z1.y - z2.y;
    return z;
}

_complex complex_real_mul(double x,_complex z1){
    _complex z;
    z.x= z1.x * x;
    z.y= z1.y * x;
    return z;
}

_complex complex_mult(_complex z1,_complex z2){
    _complex z;
    z.x= z1.x*z2.x - z1.y*z2.y;
    z.y= z1.x*z2.y + z1.y*z2.x;
    return z;
}

_complex complex_div(_complex z1,_complex z2){
    _complex z;
    z=complex_real_mul(pow(complex_abs(z2),-2),complex_mult(z1,complex_conjugate(z2)));
    return z;
}

_complex complex_pow(_complex z, int n){
    int i;
    _complex z1=z;
    for(i=1;i<n;i++){
        z1=complex_mult(z1,z);
    }
    return z1;
}

double complex_arg(_complex z){
    double cos8=(z.x/complex_abs(z));
    double sin8=(z.y/complex_abs(z));
    double arg;

    if(cos8>=0){
        if (sin8>=0) arg=acos(cos8);
    }else{
        arg=asin(sin8);
    }

    if(cos8< 0){
        if(sin8> 0){
            arg=acos(cos8);
        }
    }else{
        arg=2*acos(-1)-acos(cos8);
    }
    return arg;
}

void complex_print(_complex z){
    printf("%+.4f %+.4fi\n", z.x,z.y);
}


void complex_scanf(_complex *z){
    double x,y;
    scanf("%lf,%lf",&x,&y);
    z->x=x;
    z->y=y;
}

