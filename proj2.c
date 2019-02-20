/*
Boris Burkalo, 2018
Fakulta informačních technologií,
VUT v Brně.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define EPS 0.00000001

double taylor_log(double p, unsigned int n);
double cfrac_log(double p, unsigned int n);
double power(double x, double y, unsigned int n);
double taylor_pow(double p, double y, unsigned int n);
double taylorcf_pow(double p, double y, unsigned int n);
double mylog(double p);
double mypow(double p, double y);
int isItInfOrNan(char *arg, char *arg1);
int isItInfOrNanPow(char *arg, char *arg1, char *arg2);
void printMeLog(double c, unsigned int n);
void printMePow(double c, double d, unsigned int n);


double taylor_log(double p, unsigned int n){
  double e = 1.0;
  double sum = 0.0;
  double x = 1.0;
  double result = 0.0;
  if(p < 1.0 && p > 0){ // pro první případ Taylorova polynomu

    x = (1.0 - p);
    x = (x - 2 * x);
    while(e <= n){

      sum += x / e;
      x = x * (1 - p);
      e++;
    }
  }
  else if(p >= 1.0){ //pro druhý případ Taylorova polynomu

    for(e = 1.0; e <= n; e++){ // e se musí rovnat 1

      x *= (p - 1.0)/p;
      result = x / e;
      sum += result;      //pricteme vysledek k sum
    }
  }
  else if(p == 0){
    sum = -INFINITY;
  }
  else if(p < 0){
    sum = NAN;
  }
  return sum;
}


double cfrac_log(double p, unsigned int n){
  double cf = (2 * n - 1);
  double z;
  int n0 = n - 1;
  z = (p - 1.0)/(p + 1.0);
  double x = z * z;
  if(p > 0){
    for (; n0 >= 1; n0--){ //lomená funkce se provádí zenvnitř

      cf = (2 * n0 - 1) - (n0*n0*x / cf); //vždy se dělí předchozím výsledkem, než se n0 dostane do 1
    }
    cf = 2.0 * z / cf;
  }
  else if(p == 0){
    cf = -INFINITY;
  }
  else if(p < 0){
    cf = NAN;
  }
  return cf;
}

double power(double x, double y, unsigned int n){ //funkce ktera vypocita exponencialni fci
  double f = 1.0;
  double e = 1.0;
  double sum = 1.0;
  double pomoc = y;
  double pomoc1 = x;
  while(e < n){

    f *= e;
    sum += (y * x)/f;
    y = y * pomoc;
    x = pomoc1 * x;
    e++;
  }
  return sum;
}

double taylor_pow(double p, double y, unsigned int n){
  double sum = 1.0;
  double x = taylor_log(p, n);      //použijeme už vytvořenou funkci taylor log, kterou nahrajeme do proměnné x; do promenne f se nahrava faktorial
  sum = power(x, y, n);
  if(p == 0){
    sum = 0;
  }
  return sum;
}

double taylorcf_pow(double p, double y, unsigned int n){
  double sum = 1.0;
  double x = cfrac_log(p, n);
  sum = power(x, y, n);
  if(p == 0){
    sum = 0;
  }
  return sum;
}

double mylog(double p){
  double sum = 0.0;
  double x = 1.0;
  double controllsum = 0.0;
  unsigned int eT = 1;     //pocet iteraci provedeny v Taylorove polynomu
  unsigned int eCF = 1;    // pocet iteraci provedeny v zretezenem zlomku


  sum = taylor_log(p, eT);
  do{

    controllsum = sum;
    eT++;
    sum = taylor_log(p, eT);
  } while (fabs( sum - controllsum ) > EPS);

  x = cfrac_log(p, eCF);
  do{

    controllsum = x;
    eCF++;
    x = cfrac_log(p, eCF);
  }while(fabs( x - controllsum) > EPS);

  //nasledujici podminky uz pouze porovnavaji, ktery pocet iteraci je mensi
  if(eCF < eT){

    printf("Zvoleným způsobem je způsob zřetězených zlomků\nPotřebný počet iterací je: %d\n", eCF);
    return x;
  }
  else if(eT < eCF){

    printf("Zvoleným způsobem je způsob Taylorova polynomu\nPotřebný počet iterací je: %d\n", eT);
    return sum;
  }
  else{
    printf("Nejmenší potřebný počet iterací je pro obě možnosti stejný.\nPotřebný počet iterací je: %d\n", eT);
    return sum;
  }
}


double mypow(double p, double y){
  double sum = 0.0;
  double x = 1.0;
  double controllsum = 0.0;
  unsigned int eT = 1;
  unsigned int eCF = 1;

  sum = taylor_pow(p, y, eT);
  do {

    controllsum = sum;
    eT++;
    sum = taylor_pow(p, y, eT);
  } while(fabs( sum - controllsum ) > EPS);


  x = taylorcf_pow(p, y, eCF);
  do{

    controllsum = x;
    eCF++;
    x = taylorcf_pow(p, y, eCF);
  }while(fabs( x - controllsum) > EPS);

  if(eCF < eT){

    printf("Zvoleným způsobem je způsob zřetězených zlomků\nPotřebný počet iterací je: %d\n", eCF);
    return x;
  }
  else if (eT < eCF){

    printf("Zvoleným způsobem je způsob Taylorova polynomu\nPotřebný počet iterací je: %d\n", eT);
    return sum;
  }
  else{
    printf("Nejmenší potřebný počet iterací je pro obě možnosti stejný.\nPotřebný počet iterací je: %d\n", eT);
    return sum;
  }
}

int isItInfOrNan(char *arg, char *arg1){
  double c = 0.0;
  int n = 0;
  n = atoi(arg1);
  c = atof(arg);
  double x = 0.0;
  int help = 0;
  int len = strlen(arg1);
  int nDigits = floor(log10(abs(n))) + 1;

  if(strcmp(arg, "-nan") == 0){ //kontrola pomocí strcmp, protože strtod a atof vzdy prevede -nan na nan
    c = -NAN;
    x = -NAN;
    help = 1;
  }
  else if(isinf(c) == -1){
    x = NAN;
    help = 1;
  }
  else if(strcmp(arg, "nan") == 0 || isinf(c) == 1){
    x = c;
    help = 1;
  }
  if(help == 1 && nDigits == len){
    printf("       log(%g) = %.12g\n", c, x);
    printf(" cfrac_log(%g) = %.12g\n", c, x);
    printf("taylor_log(%g) = %.12g\n", c, x);
    return 1;
  }
  else{
    return 0;
  }
}

int isItInfOrNanPow(char *arg, char *arg1, char *arg2){
  double c = 0.0;
  double d = 0.0;
  double x = 0.0;
  int n = 0;
  int help = 0;
  n = atoi(arg2);
  c = atof(arg);
  d = atof(arg1);
  int len = strlen(arg2);
  int nDigits = floor(log10(abs(n))) + 1;
  if(strcmp(arg, "-nan") == 0){
    c = -NAN;
    x = -NAN;
    help = 1;
  }
  else if(strcmp(arg, "nan") == 0 || isinf(c) == 1 || isinf(c) == -1){
    x = fabs(c);
    help = 1;
  }
  if(help == 1 && nDigits == len){
    printf("         pow(%g,%g) = %.12g\n", c, d, x);
    printf("  taylor_pow(%g,%g) = %.12g\n", c, d, x);
    printf("taylorcf_pow(%g,%g) = %.12g\n", c, d, x);
    return 1;
  }
  else{
    return 0;
  }
}

void printMeLog(double c, unsigned int n){
  printf("       log(%g) = %.12g\n", c, log(c));
  printf(" cfrac_log(%g) = %.12g\n",c , cfrac_log(c, n));
  printf("taylor_log(%g) = %.12g\n",c , taylor_log(c, n));
}

void printMePow(double c, double d, unsigned int n){
  printf("         pow(%g,%g) = %.12g\n", c, d, pow(c, d));
  printf("  taylor_pow(%g,%g) = %.12g\n", c, d, taylor_pow(c, d, n));
  printf("taylorcf_pow(%g,%g) = %.12g\n", c, d, taylorcf_pow(c, d, n));
}

int main(int argc, char *argv[]) {
  char a[] = "--log";
  char b[] = "--pow";
  char k[] = "--mylog";
  char l[] = "--mypow";
  double c;
  double n;
  double d;
  int nDigits = 0;
  int len = 0;
  char *help;
  char *help1;


  if(argc >= 4){
    c = strtod(argv[2], &help);
    if(argc >= 3 && strcmp(argv[1], a) == 0 && argc == 4 && isItInfOrNan(argv[2], argv[3]) == 1){

      return 0;
    }
    else if(argc >= 3 && strcmp(argv[1], b) == 0 && argc == 5 && isItInfOrNanPow(argv[2], argv[3], argv[4]) == 1){

      return 0;
    }
    else if(strcmp(help, "\0") == 0){

      if (strcmp(argv[1], a) == 0 && argc == 4){

        n = atoi(argv[3]);
        len = strlen(argv[3]);
        nDigits = floor(log10(abs(n))) + 1;
        if (n <= 0 || len != nDigits){

          fprintf(stderr, "ERROR: a problem occurred!\n");
          return 0;
        }

        printMeLog(c, n);
      }
      else if (strcmp(argv[1], b) == 0 && argc == 5){

        d = strtod(argv[3], &help1);
        n = atoi(argv[4]);
        len = strlen(argv[4]);
        nDigits = floor(log10(abs(n))) + 1;
        if (n <= 0 || len != nDigits || strcmp(help1, "\0") != 0){

          fprintf(stderr, "ERROR: a problem occurred!\n");
          return 1;
        }

        printMePow(c, d, n);
      }
      else if(strcmp(argv[1], k) == 0 && argc == 3){

        printf("mylog(%g) = %.7e\n", c, mylog(c));
      }
      else if(strcmp(argv[1], l) == 0 && argc == 4 && c > 0){

        d = atof(argv[3]);
        printf("mypow(%g,%g) = %.7e\n", c, d, mypow(c, d));
      }
      else{

        fprintf(stderr, "ERROR: a problem occurred!\n");
      }
    }
    else{
      fprintf(stderr, "ERROR: a problem occurred!\n");
      return 1;
    }
  }
  else{
    fprintf(stderr, "ERROR: a problem occurred!\n");
    return 1;
  }
  return 0;
}
