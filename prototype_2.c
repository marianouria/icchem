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

/* storage locations */
struct Data {
        triplet xi;             /* ground state */
        double kac, kb;         /* half-rigidities */
} *set;
triplet *xx;
vector *rr;
double *erg;
double *diam3;
double *bps;

/* all defaults */
int NS = 0, Lk = 0;
int FREQ = 100000;
//double bps = 3.4829;            /* gives density */
//double diam2 = 0.0;
double rdata[3] = { 0.1, 0.1, 0.3 };
const double relaxer[3] = { 4.0, 4.0, 2.0 };

double bound[3];
double condit[3];
double current[3];

double simulation()
{
        int i, step = 0;
        double v = NOTHING;

       // initiate();

       // return v;
	return 2;
}


void eulerset(triplet x, double alpha, double beta, double gamma)
{
        double sa, ca, sb, cb, sc, cc;

        sin_cos(sa, ca, alpha);
        sin_cos(sb, cb, beta);
        sin_cos(sc, cc, gamma);

        x[0][0] = ca * cc - sa * cb * sc;
        x[0][1] = sa * cc + ca * cb * sc;
        x[0][2] = sb * sc;

        x[1][0] = -ca * sc - sa * cb * cc;
        x[1][1] = -sa * sc + ca * cb * cc;
        x[1][2] = sb * cc;

        x[2][0] = sa * sb;
        x[2][1] = -ca * sb;
        x[2][2] = cb;
}


void allocmemory()
{
        xx = (triplet *) malloc(sizeof(triplet) * (NS + 1));
        rr = (vector *) malloc(sizeof(vector) * (NS + 1));
        erg = (double *) malloc(sizeof(double) * (NS + 1));
        set = (struct Data *) malloc(sizeof(struct Data) * (NS + 1));
	diam3 = (double *) malloc(sizeof(double) * (NS + 1));
	bps = (double *) malloc(sizeof(double) * (NS + 1));

        if (xx == NULL || rr == NULL || erg == NULL || set == NULL || diam3 == NULL || bps == NULL) {
                printf("Not enough memory\n");
                exit(EXIT_FAILURE);
        }
}


void readdata()
{
        int i = 1, numbr;
        char c[83];
        double twist_eq = 0., arc_eq = 0.;
        double alpha, beta, gamma, sac, sb, len, diam;

        while (fgets(c, sizeof(c), stdin) != NULL) { // Verifica que el archivo no este vacio
                if (c[0] == '#' || c[0] == '\n')
                        continue; // Reinicia el bucle

                if (NS == 0 && sscanf(c, "Number of segments %d", &NS) == 1) {
                      	allocmemory();
                        continue;
                }

                if (sscanf(c, "Segments: number alpha beta gamma dac db length diam%*d")
                    == EOF) {
                        i = 1;
                        twist_eq = 0.;
                        arc_eq = 0.;
                        continue;
                }

printf("%s",c);

                if (sscanf(c, "%d %lf %lf %lf %lf %lf %lf %lf",&numbr, &alpha, &beta, &gamma, &sac, &sb, &len, &diam) == 8) {
                        for (; numbr > 0 && i <= NS; numbr--, i++) {
                                eulerset(set[i].xi, alpha, beta, gamma);
                                set[i].kac = 0.5 / (sac * sac);
                                set[i].kb = 0.5 / (sb * sb);
                                twist_eq += alpha + gamma;
                                arc_eq += beta;
				diam3[i]=diam*diam;
				bps[i]=len;
                        }
                        continue;
                }
		
printf("%s",c);

/*		if (sscanf(c, "Diameter-to-length ratio %lf", &diam2) == 1) {
                        diam2 *= diam2;
                        continue;
                }*/

                if (sscanf(c, "Linking number %d", &Lk) == 1)
                        continue;

/*                if (sscanf(c, "Segment length %lf", &bps) == 1)
                        continue;*/

printf("%s",c);
                if (sscanf(c, "Conformation restriction: tw bd r %*d") == EOF) {
                        fscanf(stdin, "%lf %lf %lf",rdata, rdata + 1, rdata + 2);
                        rdata[2] /= bps;        /* to segment length */      // <----------------------------------//??
                        continue;
                }

printf("%s",c);
                if (sscanf(c, "Number of movements %d", &FREQ) == 1)
                        continue;
        }

        NS = i - 1;

        /* report what was read */
        printf("Segments = %d   Arc = %g   Tw = %g   Lk = %d   \n", NS, arc_eq / (2 * M_PI), twist_eq / (2 * M_PI), Lk);
}


int main(int argc, char *argv[])
{
        int i, exps = 9;
        double wolume, v, mv = 0.0, ev = 0.0;

        srand((unsigned) time(NULL));

        if (argc > 1)
                if (freopen(argv[1], "rt", stdin) == NULL) {
                        printf(VER "Can't open %s\nUsage: %s [infile]\n",
                               argv[1], argv[0]);
                        return EXIT_FAILURE;
                }
        readdata();

	wolume = 1.5 / (rdata[0] * (1 - cos(rdata[1])) *rdata[2] * rdata[2] * rdata[2]); //??
        /* to moles per liter */ 
        wolume /= 6.022e23 * 3.4e-9 * 3.4e-9 * 3.4e-9 * bps * bps * bps; //<---------------------------------------//??

	/* reset parameters */
        rdata[2] = rdata[2] * rdata[2]; //<------------------------------------------------------------------------//??
        rdata[1] = 1.0 - cos(rdata[1]); //??
        rdata[0] = 1.0 - cos(rdata[0]); //??

	printf("Simulations:\n");

	for (i = 1; i <= exps; i++) {
		v = simulation();
	}

return EXIT_SUCCESS;
}

