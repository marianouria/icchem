#define VER "jfm2full 1.5, Copyright (c) 2000 - 2005 Alexei Podtelezhnikov\n"

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

#define NOTHING         1.0e+37 /* a huge number */

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef M_SQRT2
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif

/* intentionally dangerous internal macro */
#ifndef WITH_SINCOS
# define sin_cos(si, co, x)     si = sin(x); co = cos(x)
#else
# define sin_cos(si, co, x)     asm ("fsincos" : "=t" (co), "=u" (si) : "0" (x))
#endif

#define within(i, j, k)  ((i <= j && j < k) || (k < i && (i <= j || j < k)))

typedef double vector[3];
typedef vector triplet[3];
typedef double matrix[3][3];

int main (){
	
	return 0;

}
