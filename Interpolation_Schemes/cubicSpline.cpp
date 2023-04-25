#include <cstdio>
#include <vector>
#include "spline.h"
#include <iostream>
int main(int, char**) {

   std::vector<double> X ,Y;

   int N = 16;
   double a = 0;
   double b = 2;
   double dx = (b-a)/N;

   for (int i=0; i<N; i++){
     X.push_back(i*dx);
     Y.push_back(sin(2. * 3.14 * X[i]));
     std::cout << X[i] << " " << Y[i] << "\n";
   }


   tk::spline s(X,Y);
   double x=0.3, y=s(x);
   

   
   std::cout << x << " " << y << " " << sin(2.*3.14 * x)-y;

//   printf("spline at %f is %f with derivative %f\n", x, y, deriv);
}
