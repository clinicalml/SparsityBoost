# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pwl_interp_2d.h"
//# include "r8lib.hpp"

int r8vec_bracket5 ( int nd, double xd[], double xi )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
//
//  Discussion:
//
//    We assume XD is sorted.
//
//    If XI is contained in the interval [XD(1),XD(N)], then the returned 
//    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
//
//    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
//
//    This code implements a version of binary search which is perhaps more
//    understandable than the usual ones.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data values.
//
//    Input, double XD[N], the sorted data.
//
//    Input, double XD, the query value.
//
//    Output, int R8VEC_BRACKET5, the bracket information.
//
{
  int b;
  int l;
  int m;
  int r;

  if ( xi < xd[0] || xd[nd-1] < xi )
  {
    b = -1;
  }
  else
  {
    l = 0;
    r = nd - 1;

    while ( l + 1 < r )
    {
      m = ( l + r ) / 2;
      if ( xi < xd[m] )
      {
        r = m;
      }
      else
      {
        l = m;
      }
    }
    b = l;
  }

  return b;
}
//****************************************************************************80

double *pwl_interp_2d ( int nxd, int nyd, double xd[], double yd[], double zd[], 
  int ni, double xi[], double yi[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NXD, NYD, the number of X and Y data values.
//
//    Input, double XD[NXD], YD[NYD], the sorted X and Y data.
//
//    Input, double ZD[NXD*NYD}, the Z data.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], YI[NI], the coordinates of the
//    interpolation points.
//
//    Output, double PWL_INTERP_2D[NI], the value of the interpolant.
//
{
  double alpha;
  double beta;
  double det;
  double dxa;
  double dxb;
  double dxi;
  double dya;
  double dyb;
  double dyi;
  double gamma;
  int i;
  int j;
  int k;
  double *zi;

  zi = new double[ni];

  for ( k = 0; k < ni; k++ )
  {
    i = r8vec_bracket5 ( nxd, xd, xi[k] );
    if ( i == -1 )
    {
      if (xi[k] < xd[0]) {
		  i = 0;
	  }
	  if (xi[k] > xd[nxd-1]) {
		  i = nxd - 2;
	  }
		//zi[k] = r8_huge ( );
      //continue;
    }

    j = r8vec_bracket5 ( nyd, yd, yi[k] );
    if ( j == -1 )
    {
      if (yi[k] < yd[0]) {
		  j = 0;
	  }
	  if (yi[k] > yd[nyd-1]) {
		  j = nyd - 2;
	  }
      //zi[k] = r8_huge ( );
      //continue;
    }

    if ( yi[k] < yd[j+1] + ( yd[j] - yd[j+1] ) * ( xi[i] - xd[i] ) / ( xd[i+1] - xd[i] ) )
    {
      dxa = xd[i+1] - xd[i];
      dya = yd[j]   - yd[j];

      dxb = xd[i]   - xd[i];
      dyb = yd[j+1] - yd[j];

      dxi = xi[k]   - xd[i];
      dyi = yi[k]   - yd[j];

      det = dxa * dyb - dya * dxb;

      alpha = ( dxi * dyb - dyi * dxb ) / det;
      beta =  ( dxa * dyi - dya * dxi ) / det;
      gamma = 1.0 - alpha - beta;

	  //cout << "alpha: " << alpha << "beta: " << beta << "gamma: " << gamma << "det: " << endl;
      zi[k] = alpha * zd[i+1+j*nxd] + beta * zd[i+(j+1)*nxd] + gamma * zd[i+j*nxd];
    }
    else
    {
      dxa = xd[i]   - xd[i+1];
      dya = yd[j+1] - yd[j+1];

      dxb = xd[i+1] - xd[i+1];
      dyb = yd[j]   - yd[j+1];

      dxi = xi[k]   - xd[i+1];
      dyi = yi[k]   - yd[j+1];

      det = dxa * dyb - dya * dxb;

      alpha = ( dxi * dyb - dyi * dxb ) / det;
      beta =  ( dxa * dyi - dya * dxi ) / det;
      gamma = 1.0 - alpha - beta;

	  //cout << "alpha: " << alpha << "beta: " << beta << "gamma: " << gamma << "det: " << endl;
      zi[k] = alpha * zd[i+(j+1)*nxd] + beta * zd[i+1+j*nxd] + gamma * zd[i+1+(j+1)*nxd];
    }
  }
  //cout << "zi[0]: " << zi[0] << endl;
  return zi;
}
