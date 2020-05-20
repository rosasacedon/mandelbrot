/******** ROSA MARÍA SACEDÓN ORTEGA --> GRUPO LABORATORIO C1 *******/

/******* Conjunto Mandelbrot, M.
 Sea c un número complejo, n un número natural y {z(n)} una sucesión de números definida por:
 		z(0)= 0;
 		z(n+1) = z(n)² + c
 c pertenece a M ssi, para cualquier n, |z(n)| <= 2

 Como aproximación, tomaremos los valores de c tales que para cualquier n < MAX_ITER, |z(n)| <= 2

/*******    COMPILACIÓN:
 gcc mandelbrot_paralelo.c -o mandelbrot_paralelo -lm -fopenmp -O1 -Wno-unused-result

/*******    EJECUCIÓN: ./mandelbrot_paralelo <archivo de entrada>
 El archivo de entrada debe ajustarse al siguiente formato:
 	1.ª línea: <número de imágenes que queremos generar>
 	Una línea adicional por cada imagen que queramos generar con uno de los siguientes formatos según
 		la imagen que queramos sea rectangular o cuadrada, respectivamente:
 		1 <menor abscisa> <mayor abscisa> <menor ordenada> <mayor ordenada> <nombre archivo imagen>
 		2 <abscisa centro del cuadrado > <ordenada del centro> <lado del cuadrado> <nombre archivo imagen>
 		El nombre del archivo imagen no debe tener más de 15 caracteres ni punto.

 Cada imagen será de una zona rectangular (o cuadrada) cuya abscisa menor será xmin, su abscisa mayor será xmax, 
su ordenada menor será ymin y 	 su ordenada mayor será ymax.Por tanto, la anchura de la zona es xmax-xmin y la 
altura es ymax-ymin
En el archivo de imagen (que tendrá el formato ppm) el número de puntos en las direcciones horizontal y vertical 
será proporcional a  xmax-xmin e ymax-ymin, respectivamente. */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 5000	// Máximo número de términos calculados en cada sucesión
// Cuanto mayor sea MAX_ITER más calidad de imagen y mayor tiempo de ejecución
#define RK 1000			// RK*RK número aproximado de puntos de la imagen

/* La función mandel_val determina si el número x+i·y
 * pertenece (probablemente) o no (con seguridad) al conjunto de Mandelbrot */
int mandel_val(double x, double y, int max_iter){
	int j = 1;
	double temp_real, cr, p_real = x;
	double ci, p_imag = y;

	cr = x*x; ci = y*y;	// cr+ci es el cuadrado del módulo de x + i·y
	while (((cr+ci)<=4) && (j < max_iter)){
		temp_real = cr - ci + x;
		p_imag = 2 * p_real * p_imag + y;
		p_real = temp_real;
		j++;			// z(j) = p_real + i · p_imag
		cr = p_real * p_real;	ci = p_imag * p_imag;
	}

	if((cr+ci) <= 4)	return -1; else return j;
};

int main(int argc, char *argv[])	{

	int n_imag, imag, t_imag, i, j, value, valuer, valueg, valueb;
	float xmin, xmax, ymin, ymax;
	double x,y,l,t;
	char archivoppm[20];
	FILE* filein;
    int hilos = omp_get_num_procs();

	// A. Abrir archivo de entrada y leer número de imágenes
	if (argc == 2) {
	   filein = fopen(argv[1], "r");
	   if (filein == NULL) {
		  fprintf(stderr, "No puedo abrir %s\n", argv[1]);
		  fprintf(stderr, "uso: %s <archivo de entrada> \n",argv[0]);
		  exit(0);
	   }
	   fscanf(filein,"%d",&n_imag);
	}	else {
			fprintf(stderr, "uso: %s <archivo de entrada> \n",argv[0]);
			exit(0);
		}

	// B. Obtener una a una las imágenes y los tiempos respectivos
	for (imag = 0; imag <n_imag; imag ++)	{

		t = omp_get_wtime();

		// B1. Leer tipo de imagen (1 = rectangular, 2 = cuadrada)
		fscanf(filein,"%d",&t_imag);
		if (t_imag == 1)	{
			fscanf(filein,"%f %f %f %f %s", &xmin, &xmax, &ymin, &ymax, archivoppm);
		}	else if (t_imag == 2)	{
				fscanf(filein,"%lf %lf %lf %s", &x, &y, &l, archivoppm);
				l = l/2; xmin = x-l; xmax = x+l; ymin = y-l; ymax = y+l;
			}	else {fprintf(stderr, "\n Error en formato del archivo de entrada\n");exit(0);}

		// B2. Calcular, en pixels, la anchura, width, y la altura, height, de la imagen
		double k = RK/sqrt((xmax-xmin)*(ymax-ymin));
		int width = k*(xmax-xmin);	// Hay conversión de double a int
		int height = k*(ymax-ymin);	// Íd.
		// width x height = RK² aproximadamente
		// k = width/(xmax-xmin) = height/(ymax-ymin);

		// B3. Se abre y configura el archivo de imagen ppm
		strcat(archivoppm,".ppm"); //archivoppm = archivoppm + ".ppm";
		FILE *f_imag = fopen(archivoppm, "wb");
		fprintf(f_imag, "P3\n%d %d\n255\n", width, height);
		/* Formato P3: cada pixel se representa por sus valores de rojo, verde y azul (entre 0 y 255);
		 * width x height: dimensión imagen en pixels */

		// B4. Se calcula el punto (x,y) y se determina si incluir (x + i·y) en el conjunto de Mandelbrot
		int array [height][width];

        #pragma omp parallel num_threads(hilos) shared(xmin, ymax, k, array)
        #pragma omp for schedule(dynamic,80) private(x, y, j)
		for(i = 0; i < height; i++){
			for(j = 0; j < width; j++) {
				x = xmin + j/k;
				y = ymax - i/k;
                array[i][j] = mandel_val(x, y, MAX_ITER);
			}
		}

		// B5. Se determina el color del pixel del número (x + y·i) y se imprime en el archivo
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                if(array[i][j] == -1)	// se supone que (x+i·y) pertenece a M --> pixel blanco
					fprintf(f_imag, " 255 255 255 ");
				else
				// se sabe que (x+i·y) no pertenece a M --> color del pixel según número de iteraciones
				{	value = ((float)array[i][j]/MAX_ITER) * 16777215;
					// Hacemos que value quede entre 0 y 2²⁴ y extraemos sus componentes rgb
					valueg = value >> 8;	valueb = value % 256;
					valuer = valueg >> 8;   valueg = valueg % 256;
					fprintf(f_imag," %d %d %d ",valuer,valueg,valueb);
				}
            }
        }

		fprintf(f_imag, "\n");	fclose(f_imag);

		// B6. Informamos al usuario
		printf("\n En el archivo %s se ha creado una imagen de la zona explorada", archivoppm);
		printf("\n\n Tiempo empleado %0.2f segundos\n\n",omp_get_wtime()-t);
        
	}	// fin de: for (imag = 0; imag <n; imag ++)
	fclose(filein);
}

/******* PROCESADOR: Intel Core i5-6200U 2.30GHz ********/
/******* NÚCLEOS UTILIZADOS: 4 *******/

/******* RESULTADOS DE COMPILACIÓN *******/

/*********************************************************************************/
/************                                                   ******************/
/************             CUADRADO                              ******************/
/*********************************************************************************/
/************                                                   ******************/
/************   SECUENCIAL  --> 16.40 segundos                  ******************/
/************                                                   ******************/
/************   PARALELO (sin shedule)    --> 5.68 segundos     ******************/
/************                                                   ******************/
/*********************************************************************************/
/************                                                   ******************/
/************            RECTÁNGULO                             ******************/
/*********************************************************************************/
/************                                                   ******************/
/************   SECUENCIAL  -->  8.08 segundos                  ******************/
/************                                                   ******************/
/************   PARALELO (sin shedule)   -->  3.71 segundos     ******************/  
/************                                                   ******************/                                                       
/*********************************************************************************/
/*********************************************************************************/
/************                                                   ******************/
/************   PARALELOS CON SHEDULE                           ******************/
/************                                                   ******************/
/*********************************************************************************/
/************             CUADRADO                              ******************/
/*********************************************************************************/
/************   STATIC,1 -->4.50          DYNAMIC,1 -->4.57     ******************/       
/************                                                   ******************/
/************   STATIC,10 -->4.70         DYNAMIC,10 -->5.23    ******************/
/************                                                   ******************/
/************   STATIC,20 -->4.64         DYNAMIC,20 -->5.70    ******************/
/************                                                   ******************/
/************   STATIC,30 -->5.14         DYNAMIC,30 -->5.33    ******************/
/************                                                   ******************/
/************   STATIC,40 -->4.86         DYNAMIC,40 -->5.23    ******************/
/************                                                   ******************/
/************   STATIC,50 -->4.89         DYNAMIC,50 -->5.47    ******************/
/************                                                   ******************/
/************   STATIC,60 -->5.08         DYNAMIC,60 -->5.07    ******************/
/************                                                   ******************/
/************   STATIC,70 -->4.74         DYNAMIC,70 -->5.47    ******************/
/************                                                   ******************/
/************   STATIC,80 -->5.06         DYNAMIC,80 -->5.11    ******************/
/************                                                   ******************/
/************   STATIC,90 -->5.46         DYNAMIC,90 -->5.11    ******************/
/************                                                   ******************/
/*********************************************************************************/
/************            RECTÁNGULO                             ******************/
/*********************************************************************************/
/************   STATIC,1 -->1.91          DYNAMIC,1 -->1.98     ******************/   
/************                                                   ******************/
/************   STATIC,10 -->2.01         DYNAMIC,10 -->2.16    ******************/
/************                                                   ******************/
/************   STATIC,20 -->2.07         DYNAMIC,20 -->2.03    ******************/
/************                                                   ******************/
/************   STATIC,30 -->2.00         DYNAMIC,30 -->2.16    ******************/
/************                                                   ******************/
/************   STATIC,40 -->2.01         DYNAMIC,40 -->2.51    ******************/
/************                                                   ******************/
/************   STATIC,50 -->2.04         DYNAMIC,50 -->2.21    ******************/
/************                                                   ******************/
/************   STATIC,60 -->1.92         DYNAMIC,60 -->2.24    ******************/
/************                                                   ******************/
/************   STATIC,70 -->2.03         DYNAMIC,70 -->2.07    ******************/
/************                                                   ******************/
/************   STATIC,80 -->2.08         DYNAMIC,80 -->1.90    ******************/
/************                                                   ******************/
/************   STATIC,90 -->2.63         DYNAMIC,90 -->1.93    ******************/
/************                                                   ******************/
/*********************************************************************************/

/************                   CONCLUSIÓN                      ******************/
/*Tras realizar las ejecuciones, observando el codigo y funcionamiento del mismo podemos deducir 
que en las zonas mas proximas a los bordes de las figuras es donde mas tiempo tarda la ejecucion.
Tenemos un cuadrado y un rectangulo, por tanto gracias a la prediccion dinamica obtenemos unos 
resultados con un menor tiempo de ejecucion. Sin embargo, si solo tuviesemos una de las imagenes, 
o bien cuadrado o bien rectangulo, la prediccion estatica daria unos mejores resultados. Estas 
diferencias nacen en el momento en el que se realiza la prediccion, en la estatica se realiza en 
tiempo de compilacion, mientras que en la dinamica se realiza en tiempo de ejecucion. 

Gracias a lo estudiado en el tema, podemos decir que la prediccion dinamica tiene una mayor precision
gracias a un coste mayor del hardware. La prediccion estatica tiene menos precision, un coste menor 
del hardware, pero tiene algunas implementaciones que pueden llegar a incrementar la longitud del ejecutable. 
La efectividad de cada prediccion depende del numero de aciertos, es decir, si ha habido mas fallos puede 
aumentar el tiempo de ejecucion. 
*/

