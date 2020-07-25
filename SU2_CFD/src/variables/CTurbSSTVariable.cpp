/*!
 * \file CTurbSSTVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */


#include "../../include/variables/CTurbSSTVariable.hpp"


CTurbSSTVariable::CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint, unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = kine;
    Solution(iPoint,1) = omega;
  }

  Solution_Old = Solution;

  sigma_om2 = constants[3];
  beta_star = constants[6];

  F1.resize(nPoint) = su2double(1.0);
  F2.resize(nPoint) = su2double(0.0);
  CDkw.resize(nPoint) = su2double(0.0);

  muT.resize(nPoint) = mut;

  if (config->GetUsing_SDD()) {
    Aij_ML.resize(nPoint,3,3);
    Aij_BL.resize(nPoint,3,3);
  }  

}

void CTurbSSTVariable::SetBlendingFunc(unsigned long iPoint, su2double val_viscosity,
                                       su2double val_dist, su2double val_density) {
  su2double arg2, arg2A, arg2B, arg1;

  AD::StartPreacc();
  AD::SetPreaccIn(val_viscosity);  AD::SetPreaccIn(val_dist);
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(Solution[iPoint], nVar);
  AD::SetPreaccIn(Gradient[iPoint], nVar, nDim);

  /*--- Cross diffusion ---*/

  CDkw(iPoint) = 0.0;
  for (unsigned long iDim = 0; iDim < nDim; iDim++)
    CDkw(iPoint) += Gradient(iPoint,0,iDim)*Gradient(iPoint,1,iDim);
  CDkw(iPoint) *= 2.0*val_density*sigma_om2/Solution(iPoint,1);
  CDkw(iPoint) = max(CDkw(iPoint), pow(10.0, -20.0));

  /*--- F1 ---*/

  arg2A = sqrt(Solution(iPoint,0))/(beta_star*Solution(iPoint,1)*val_dist+EPS*EPS);
  arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Solution(iPoint,1)+EPS*EPS);
  arg2 = max(arg2A, arg2B);
  arg1 = min(arg2, 4.0*val_density*sigma_om2*Solution(iPoint,0) / (CDkw(iPoint)*val_dist*val_dist+EPS*EPS));
  F1(iPoint) = tanh(pow(arg1, 4.0));

  /*--- F2 ---*/

  arg2 = max(2.0*arg2A, arg2B);
  F2(iPoint) = tanh(pow(arg2, 2.0));

  AD::SetPreaccOut(F1(iPoint)); AD::SetPreaccOut(F2(iPoint)); AD::SetPreaccOut(CDkw(iPoint));
  AD::EndPreacc();

}

void CTurbSSTVariable::InitAijML(unsigned long iPoint, su2double muT, su2double turb_ke, su2double rho, su2double **PrimGrad, su2double *delta_sdd) {

  unsigned short iDim, jDim;
  su2double **S_ij          = new su2double* [3];
  su2double **A_ij          = new su2double* [3];
  su2double **newA_ij       = new su2double* [3];
  su2double **delta3        = new su2double* [3];
  su2double **Eig_Vec       = new su2double* [3];
  su2double **Corners       = new su2double* [3];
  su2double *Eig_Val        = new su2double  [3];
  su2double *Bary_Coord     = new su2double  [2];
  su2double *New_Bary_Coord = new su2double  [2];
  for (iDim = 0; iDim < 3; iDim++){
    A_ij[iDim]         = new su2double [3];
    newA_ij[iDim]      = new su2double [3];
    S_ij[iDim]         = new su2double [3];
    delta3[iDim]       = new su2double [3] ();
    delta3[iDim][iDim] = 1.0;
    Eig_Vec[iDim]      = new su2double [3];
    Eig_Val[iDim]      = 0.0;
    Corners[iDim]      = new su2double [2];
  }
  su2double divVel = 0;

 /* --- Calculate rate of strain tensor --- */
  for (iDim = 0; iDim < 3; iDim++){
    divVel += PrimGrad[iDim+1][iDim];
  }

  for (iDim = 0; iDim < 3; iDim++) {
    for (jDim = 0 ; jDim < 3; jDim++) {
      S_ij[iDim][jDim] = 0.5*(PrimGrad[iDim+1][jDim] + PrimGrad[jDim+1][iDim]) - divVel*delta3[iDim][jDim]/3.0;
    }
  }

  /* --- Calculate anisotropic part of Reynolds Stress tensor --- */
  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      A_ij[iDim][jDim] = - muT * S_ij[iDim][jDim] / (rho*turb_ke);
      Aij_BL(iPoint,iDim,jDim) = A_ij[iDim][jDim];
//      cout << iDim << jDim << A_ij[iDim][jDim] << endl;
    }
  }

  /* --- Get ordered eigenvectors and eigenvalues of A_ij --- */
  EigenDecomposition(A_ij, Eig_Vec, Eig_Val, 3);

  /* compute convex combination coefficients */
  su2double c1c = Eig_Val[2] - Eig_Val[1];
  su2double c2c = 2.0 * (Eig_Val[1] - Eig_Val[0]);
  su2double c3c = 3.0 * Eig_Val[0] + 1.0;
  //cout << "old" << " " << c1c << " " << c2c << " " << c3c << endl;

  /* define barycentric traingle corner points */
  Corners[0][0] = 1.0;
  Corners[0][1] = 0.0;
  Corners[1][0] = 0.0;
  Corners[1][1] = 0.0;
  Corners[2][0] = 0.5;
  Corners[2][1] = 0.866025;

  /* define baseline barycentric coordinates */
  Bary_Coord[0] = Corners[0][0] * c1c + Corners[1][0] * c2c + Corners[2][0] * c3c;
  Bary_Coord[1] = Corners[0][1] * c1c + Corners[1][1] * c2c + Corners[2][1] * c3c;

  /* Add ML derived delta to barycentric coordinates */
  New_Bary_Coord[0] = Bary_Coord[0] + delta_sdd[0];
  New_Bary_Coord[1] = Bary_Coord[1] + delta_sdd[1];

  /* rebuild c1c,c2c,c3c based on perturbed barycentric coordinates */
  c3c = New_Bary_Coord[1] / Corners[2][1];
  c1c = New_Bary_Coord[0] - Corners[2][0] * c3c;
  c2c = 1 - c1c - c3c;
  //cout << "new" << " " << c1c << " " << c2c << " " << c3c << endl;

  //TODO - limit c1c, c2c, c3c to be realizable here?

  /* build new anisotropy eigenvalues */
  Eig_Val[0] = (c3c - 1) / 3.0;
  Eig_Val[1] = 0.5 *c2c + Eig_Val[0];
  Eig_Val[2] = c1c + Eig_Val[1];

  // TODO - ML derived perturbations to eigenvectors would go here
  EigenRecomposition(newA_ij, Eig_Vec, Eig_Val, 3);

  // Save newA_ij in field variable AijML and original A_ij in AijBL
  for (iDim = 0; iDim < 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      Aij_ML(iPoint,iDim,jDim) = newA_ij[iDim][jDim];
    }	    
  }

  /* --- Deallocate memory --- */
    for (iDim = 0; iDim < 3; iDim++){
      delete [] S_ij[iDim];
      delete [] A_ij[iDim];
      delete [] newA_ij[iDim];
      delete [] delta3[iDim];
      delete [] Eig_Vec[iDim];
      delete [] Corners[iDim];
    }
    delete [] S_ij;
    delete [] A_ij;
    delete [] newA_ij;
    delete [] delta3;
    delete [] Eig_Vec;
    delete [] Corners;
    delete [] Eig_Val;
    delete [] Bary_Coord;
    delete [] New_Bary_Coord;

}

void CTurbSSTVariable::EigenDecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val, unsigned short n){
  int iDim,jDim;
  su2double *e = new su2double [n];
  for (iDim= 0; iDim< n; iDim++){
    e[iDim] = 0;
    for (jDim = 0; jDim < n; jDim++){
      Eig_Vec[iDim][jDim] = A_ij[iDim][jDim];
    }
  }
  tred2(Eig_Vec, Eig_Val, e, n);
  tql2(Eig_Vec, Eig_Val, e, n);

  delete [] e;

}

void CTurbSSTVariable::EigenRecomposition(su2double **A_ij, su2double **Eig_Vec, const su2double *Eig_Val, unsigned short n){
  unsigned short i,j,k;
  su2double **tmp = new su2double* [n];
  su2double **deltaN = new su2double* [n];

  for (i= 0; i< n; i++){
    tmp[i] = new su2double [n];
    deltaN[i] = new su2double [n];
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) deltaN[i][j] = 1.0;
      else deltaN[i][j]=0.0;
    }
  }

  for (i= 0; i< n; i++){
    for (j = 0; j < n; j++){
      tmp[i][j] = 0.0;
      for (k = 0; k < n; k++){
        tmp[i][j] += Eig_Vec[i][k] * Eig_Val[k] * deltaN[k][j];
      }
    }
  }

  for (i= 0; i< n; i++){
    for (j = 0; j < n; j++){
      A_ij[i][j] = 0.0;
      for (k = 0; k < n; k++){
        A_ij[i][j] += tmp[i][k] * Eig_Vec[j][k];
      }
    }
  }

  for (i = 0; i < n; i++){
    delete [] tmp[i];
    delete [] deltaN[i];
  }
  delete [] tmp;
  delete [] deltaN;
}

void CTurbSSTVariable::tred2(su2double **V, su2double *d, su2double *e, unsigned short n) {
/* Author:

 * Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
 * Klema, Moler.
 * C++ version by Aashwin Mishra and Jayant Mukhopadhaya.

 * Reference:

 * Martin, Reinsch, Wilkinson,
 * TRED2,
 * Numerische Mathematik,
 * Volume 11, pages 181-195, 1968.

 * James Wilkinson, Christian Reinsch,
 * Handbook for Automatic Computation,
 * Volume II, Linear Algebra, Part 2,
 * Springer, 1971,
 * ISBN: 0387054146,
 * LC: QA251.W67.

 * Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
 * Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
 * Matrix Eigensystem Routines, EISPACK Guide,
 * Lecture Notes in Computer Science, Volume 6,
 * Springer Verlag, 1976,
 * ISBN13: 978-3540075462,
 * LC: QA193.M37

*/

  unsigned short i,j,k;

  for (j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  /* Householder reduction to tridiagonal form. */

  for (i = n-1; i > 0; i--) {

    /* Scale to avoid under/overflow. */

    su2double scale = 0.0;
    su2double h = 0.0;
    for (k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    }
    else {

      /* Generate Householder vector. */

      for (k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      su2double f = d[i-1];
      su2double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      /* Apply similarity transformation to remaining columns. */

      for (j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      su2double hh = f / (h + h);
      for (j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (k = j; k <= i-1; k++) {
            V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  /* Accumulate transformations. */

  for (i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    su2double h = d[i+1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (j = 0; j <= i; j++) {
        su2double g = 0.0;
        for (k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
}

void CTurbSSTVariable::tql2(su2double **V, su2double *d, su2double *e, unsigned short n) {

/* Author:

 * Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
 * Klema, Moler.
 * C++ version by Aashwin Mishra and Jayant Mukhopadhaya.

 * Reference:

 * Bowdler, Martin, Reinsch, Wilkinson,
 * TQL2,
 * Numerische Mathematik,
 * Volume 11, pages 293-306, 1968.

 * James Wilkinson, Christian Reinsch,
 * Handbook for Automatic Computation,
 * Volume II, Linear Algebra, Part 2,
 * Springer, 1971,
 * ISBN: 0387054146,
 * LC: QA251.W67.

 * Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
 * Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
 * Matrix Eigensystem Routines, EISPACK Guide,
 * Lecture Notes in Computer Science, Volume 6,
 * Springer Verlag, 1976,
 * ISBN13: 978-3540075462,
 * LC: QA193.M37

*/

  int i,j,k,l;
  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  su2double f = 0.0;
  su2double tst1 = 0.0;
  su2double eps = pow(2.0,-52.0);
  for (l = 0; l < n; l++) {

    /* Find small subdiagonal element */

    tst1 = max(tst1,(fabs(d[l]) + fabs(e[l])));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    /* If m == l, d[l] is an eigenvalue, */
    /* otherwise, iterate.               */

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  /* (Could check iteration count here.) */

        /* Compute implicit shift */

        su2double g = d[l];
        su2double p = (d[l+1] - g) / (2.0 * e[l]);
        su2double r = sqrt(p*p+1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        su2double dl1 = d[l+1];
        su2double h = g - d[l];
        for (i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        /* Implicit QL transformation. */

        p = d[m];
        su2double c = 1.0;
        su2double c2 = c;
        su2double c3 = c;
        su2double el1 = e[l+1];
        su2double s = 0.0;
        su2double s2 = 0.0;
        for (i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = sqrt(p*p+e[i]*e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          /* Accumulate transformation. */

          for (k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        /* Check for convergence. */

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  /* Sort eigenvalues and corresponding vectors. */

  for (i = 0; i < n-1; i++) {
    k = i;
    su2double p = d[i];
    for (j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < n; j++) {
          p = V[j][i];
          V[j][i] = V[j][k];
          V[j][k] = p;
      }
    }
  }
}

