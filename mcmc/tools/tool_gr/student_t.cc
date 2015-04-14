#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

// Definitions

typedef double (*func_v)(vector<double>&, double);

inline double tsign(double x, double y)
{
  if (y>=0) 
    return  fabs(x);
  else 
    return -fabs(x);
}


//
// Brent's method
//
double zbrent(func_v func, vector<double>& v, double x1, double x2, double tol)
{
  const int ITMAX = 500;
  const double EPS=3.0e-8;

  double a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm;

  a = x1;
  b = x2;
  fa = (*func)(v,a);
  fb = (*func)(v,b);
  e = d = b - a;

  if (fa*fb > 0.0) {
    cerr << "root must be bracketed for zbrent" << endl;
    return FP_NAN;
  }

  c  = b;
  fc = fb;
  for (int iter=0; iter<ITMAX; iter++) {
    if (fb*fc > 0.0) {
      c  = a;
      fc = fa;
      d  = b-a;
      e  = d;
    }
    if(fabs(fc) < fabs(fb)) {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
    xm = 0.5*(c-b);

    if(fabs(xm) <= tol1 || fb == 0.0) return b;

    if(fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa;

      if(a == c) {
	p = 2.0*xm*s;
	q = 1.0-s;
      } else {
	q = fa/fc;
	r = fb/fc;
	p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q = (q-1.0)*(r-1.0)*(s-1.0);
      }

      if(p > 0.0) q=-q;
      p = fabs(p);
      if(2.0*p < min<double>(3.0*xm*q-abs(tol1*q),abs(e*q))) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a = b;
    fa = fb;
    if(fabs(d) > tol1) 
      b = b + d;
    else
      b = b + tsign(tol1, xm);

    fb = (*func)(v,b);
  }

  cerr << "zbrent exceeding maximum iterations" << endl;

  return b;
}      


//  Continued fraction evaluation.
//  Used by routine betai.
//  Numerical Recipes, chapter 6.4.
double betacf(double a, double b, double x)
{
  const int MAXIT = 100;
  const double EPS=3.0e-7, FPMIN=1.0e-30;

  int m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c   = 1.0;
  d   = 1.0 - qab*x/qap;
  if(fabs(d) < FPMIN) d=FPMIN;
  d = 1.0/d;
  h = d;
  for (int m=1; m<=MAXIT; m++) {
    m2 = 2*m;
    aa = m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1.0 + aa*d;
    if(fabs(d) < FPMIN)d=FPMIN;
    c = 1.0 + aa/c;
    if(fabs(c) < FPMIN)c=FPMIN;
    d = 1.0/d;
    h = h*d*c;
    aa =-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d;
    if(fabs(d) < FPMIN)d=FPMIN;
    c = 1.0+aa/c;
    if(abs(c) < FPMIN)c=FPMIN;
    d = 1.0/d;
    del = d*c;
    h = h*del;
    if(fabs(del-1.0) < EPS) return h;
  }
   
  cerr << "a or b too big, or MAXIT too small in betacf" << endl;
  return h;
}      

// Incomplete beta function.
double betai(double a, double b, double x)
{
  double bt;
  if(x < 0.0 || x > 1.0)
    cerr << "bad argument x in betai" << endl;
  if(x == 0.0 || x == 1.0) 
    bt=0.0;
  else
    bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
  
  if(x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}


double tfct(vector<double>& p, double x)
{
  return betai(p[0], p[1], x) - 1.0 + p[2];
}

/*
  The two-sided inverse Student's T distribution. This is solved
  using Brent's method by solving the relation:

  alpha = 1 - IncompleteBetaFunction( x, a, b ) 

  x = df / ( df + student_t^2 )
  a = df / 2
  b = 1/2

*/
double inv_student_t_2sided(double alpha, double df)
{
  double x1=0.0, x2=1.0, eps=1.0e-7;
  vector<double> p(3);

  p[0] = 0.5 * df;
  p[1] = 0.5;
  p[2] = fabs(alpha);

  double student_t = zbrent(tfct, p, x1, x2, eps );
  student_t = sqrt( df/student_t  - df );

  return student_t;
}

/*
  The one-sided version of the inverse Student's T distribution.
*/
double inv_student_t_1sided(double alpha, double df)
{
  double x1=0.0, x2=1.0, eps=1.0e-7;

  vector<double> p(3);
  p[0] = 0.5*df;
  p[1] = 0.5;
  p[2] = fabs(alpha-0.5)*2.0;

  double student_t = zbrent(tfct, p, x1, x2, eps );
  student_t = sqrt( df/student_t  - df );

  return tsign(student_t, alpha);
}


#ifdef MAIN
int main()
{
  cout << "alpha, DF? ";
  double df, alpha;
  cin >> alpha;
  cin >> df;

  cout << "Value=" << inv_student_t_1sided(alpha, df) << endl;
  return 0;
}

#endif
