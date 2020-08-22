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
    TKE_ML.resize(nPoint);
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

void CTurbSSTVariable::InitSDD(unsigned long iPoint, su2double muT, su2double turb_ke, su2double rho, su2double **PrimGrad, su2double *delta_sdd, su2double dist, su2double *coords) {

  unsigned short iDim, jDim;
  unsigned int *idx = new unsigned int[3];
  su2double **S_ij          = new su2double* [3];
  su2double **A_ij          = new su2double* [3];
  su2double **newA_ij       = new su2double* [3];
  su2double **delta3        = new su2double* [3];
  su2double **Eig_Vec       = new su2double* [3];
  su2double **Corners       = new su2double* [3];
  su2double *Eig_Val        = new su2double  [3];
  su2double *Bary_Coord     = new su2double  [2];
  su2double *New_Bary_Coord = new su2double  [2];
  su2double *h_r            = new su2double  [4];
  su2double *h_orig         = new su2double  [4];
  su2double *h_new          = new su2double  [4];
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
  su2double r3o2  = sqrt(3.0)/2.0;
  su2double r3i   = 1.0/sqrt(3.0);

 /* --- Calculate rate of strain tensor --- */
//  for (iDim = 0; iDim < 3; iDim++){
//    divVel += PrimGrad[iDim+1][iDim];
//  }
//  divVel = 0.0; // Hard code divVel = 0 as incompressible flow

  for (iDim = 0; iDim < 3; iDim++) {
    for (jDim = 0 ; jDim < 3; jDim++) {
      S_ij[iDim][jDim] = 0.5*(PrimGrad[iDim+1][jDim] + PrimGrad[jDim+1][iDim]) - divVel*delta3[iDim][jDim]/3.0;
    }
  }

  /* --- Calculate anisotropic part of Reynolds Stress tensor --- */
  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      A_ij[iDim][jDim] = - muT * S_ij[iDim][jDim] / (rho*turb_ke);
    }
  }
  A_ij[2][0] = 0.0; 
  A_ij[2][1] = 0.0; 
  A_ij[0][2] = 0.0; 
  A_ij[1][2] = 0.0; 

  /* --- Get ordered eigenvectors and eigenvalues of A_ij --- */
  EigenDecomposition(A_ij, Eig_Vec, Eig_Val, 3);
//  if (coords[0] > 4.64 && coords[0] < 4.65 && coords[1] > 0.61 && coords[1]<0.62) {
  //if (coords[0] > 5.92 && coords[0] < 5.94 && coords[1] > 0.65 && coords[1]<0.66) {
//  if (coords[0] > 5.96 && coords[0] < 5.97 && coords[1] > 0.64 && coords[1]<0.65) {

  cout << " " << endl;
  cout << "coord" << " " << coords[0] << ", " << coords[1]  << endl;
  cout << "OLD" << endl;
  cout << "EigVal: " << Eig_Val[0] << " " << Eig_Val[1] << " " << Eig_Val[2] << endl;
  cout << "V1: " << Eig_Vec[0][0] << " " << Eig_Vec[0][1] << " " << Eig_Vec[0][2] << endl;
  cout << "V2: " << Eig_Vec[1][0] << " " << Eig_Vec[1][1] << " " << Eig_Vec[1][2] << endl;
  cout << "V3: " << Eig_Vec[2][0] << " " << Eig_Vec[2][1] << " " << Eig_Vec[2][2] << endl;

  /* Find index to sort Eigenvalues. Eigenvalues and Eigenvectors themselves are not 
   * sorted here because we need the unsorted Eigenvectors in q_from_mat(). We only need sorted 
   * ones fro c1c,c2c and c3c calc. (and rebuilding after c1c perturb etc) */
  /* TODO Sort Eigenvalues only here */
  idx = EigenSort(Eig_Val);
  cout << "idx" << " " << idx[0] << " " << idx[1] << " " << idx[2] << endl;

  /* compute convex combination coefficients */
  su2double c1c = Eig_Val[idx[2]] - Eig_Val[idx[1]];
  su2double c2c = 2.0 * (Eig_Val[idx[1]] - Eig_Val[idx[0]]);
  su2double c3c = 3.0 * Eig_Val[idx[0]] + 1.0;
  cout << "old" << " " << c1c << " " << c2c << " " << c3c << endl;

  /* define barycentric traingle corner points */
  Corners[0][0] = 1.0;
  Corners[0][1] = 0.0;
  Corners[1][0] = 0.0;
  Corners[1][1] = 0.0;
  Corners[2][0] = 0.5;
  Corners[2][1] = 0.866025;

//  delta_sdd[0] = 0.403;
//  delta_sdd[1] = -0.38;
//  delta_sdd[2] = 0;
//  delta_sdd[3] = 0;
//  delta_sdd[4] = -0.9;

  /* define baseline barycentric coordinates */
  Bary_Coord[0] = Corners[0][0] * c1c + Corners[1][0] * c2c + Corners[2][0] * c3c;
  Bary_Coord[1] = Corners[0][1] * c1c + Corners[1][1] * c2c + Corners[2][1] * c3c;

  /* Add ML derived delta to barycentric coordinates */
  New_Bary_Coord[0] = Bary_Coord[0] + delta_sdd[0];
  New_Bary_Coord[1] = Bary_Coord[1] + delta_sdd[1];
  cout << "delta x,y" << " " << delta_sdd[0] << ", " << delta_sdd[1]  << endl;

  /* limit perturbation to be inside barycentric triagle for realizability */
  New_Bary_Coord[1] = max(0.0,min(r3o2,New_Bary_Coord[1]));  // eta constraint: 0 <= eta <= sqrt(3)/2
  New_Bary_Coord[0] = max(r3i*New_Bary_Coord[1],min(1.0-r3i*New_Bary_Coord[1],New_Bary_Coord[0])); // zeta constraint: eta/sqrt(3) <= zeta <= 1-(eta/sqrt(3)

//  /* rebuild c1c,c2c,c3c based on perturbed barycentric coordinates */
  c3c = New_Bary_Coord[1] / Corners[2][1];
  c1c = New_Bary_Coord[0] - Corners[2][0] * c3c;
  c2c = 1 - c1c - c3c;
  cout << "new" << " " << c1c << " " << c2c << " " << c3c << endl;

  /* build new anisotropy eigenvalues */
//  Eig_Val[0] = (c3c - 1) / 3.0;
//  Eig_Val[1] = 0.5 *c2c + Eig_Val[0];
//  Eig_Val[2] = c1c + Eig_Val[1];
  Eig_Val[1] = (c3c - 1) / 3.0;
  Eig_Val[2] = 0.5 *c2c + Eig_Val[1];
  Eig_Val[0] = c1c + Eig_Val[2];

  /* Get rotation quaternion */
  h_r[1] = delta_sdd[2];
  h_r[2] = delta_sdd[3];
  h_r[3] = delta_sdd[4];
  h_r[0] = sqrt(1.0 - (h_r[1]*h_r[1] + h_r[2]*h_r[2] + h_r[3]*h_r[3])); // assuming unit quarternion

  /* Get original quaternion */
  q_from_mat(Eig_Vec,h_orig);
  q_norm(h_orig); //h_r should already be unit length so only norm h_orig

  cout << "h_orig: " << " " << h_orig[0] << " " << h_orig[1] << " " << h_orig[2] << " " << h_orig[3] << endl;
  cout << "h_r:    " << " " << h_r[0] << " " << h_r[1] << " " << h_r[2] << " " << h_r[3] << endl;

  /* Apply rotation quaternion to original to get new rotated one */
  q_prod(h_r,h_orig,h_new);
  cout << "h_new:    " << " " << h_new[0] << " " << h_new[1] << " " << h_new[2] << " " << h_new[3] << endl;

  /* Convert new quaternion back to eigenvector */ 
  mat_from_q(h_new,Eig_Vec);
  cout << "NEW" << endl;
  cout << "EigVal: " << Eig_Val[0] << " " << Eig_Val[1] << " " << Eig_Val[2] << endl;
  cout << "V1: " << Eig_Vec[0][0] << " " << Eig_Vec[0][1] << " " << Eig_Vec[0][2] << endl;
  cout << "V2: " << Eig_Vec[1][0] << " " << Eig_Vec[1][1] << " " << Eig_Vec[1][2] << endl;
  cout << "V3: " << Eig_Vec[2][0] << " " << Eig_Vec[2][1] << " " << Eig_Vec[2][2] << endl;

  /* Rebuild new Aij tensor */
//  Eig_Val[0] = 0.54;
//  Eig_Val[1] = -0.3;
//  Eig_Val[2] = -0.238;
//  Eig_Vec[0][0] = 0.97; Eig_Vec[0][1] = 0.214; Eig_Vec[0][2] = 0.0;
//  Eig_Vec[1][0] = -0.21; Eig_Vec[1][1] = 0.976; Eig_Vec[1][2] = 0.0;
//  Eig_Vec[2][0] = 0.0; Eig_Vec[2][1] = 0.0; Eig_Vec[2][2] = 1.0;

  EigenRecomposition(newA_ij, Eig_Vec, Eig_Val, 3);

  // Save newA_ij in field variable AijML
  for (iDim = 0; iDim < 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      Aij_ML(iPoint,iDim,jDim) = newA_ij[iDim][jDim];
    }	    
  }
  cout << "A1j_old: " << A_ij[0][0] << " " << A_ij[0][1] << " " << A_ij[0][2] << endl;
  cout << "A2j_old: " << A_ij[1][0] << " " << A_ij[1][1] << " " << A_ij[1][2] << endl;
  cout << "A3j_old: " << A_ij[2][0] << " " << A_ij[2][1] << " " << A_ij[2][2] << endl;
  cout << "A1j_new: " << newA_ij[0][0] << " " << newA_ij[0][1] << " " << newA_ij[0][2] << endl;
  cout << "A2j_new: " << newA_ij[1][0] << " " << newA_ij[1][1] << " " << newA_ij[1][2] << endl;
  cout << "A3j_new: " << newA_ij[2][0] << " " << newA_ij[2][1] << " " << newA_ij[2][2] << endl;

  //Aij_ML(iPoint,0,0) = h_orig[0]; 
  //Aij_ML(iPoint,0,1) = h_orig[3]; 
  //Aij_ML(iPoint,0,2) = h_r[0]; 
  //Aij_ML(iPoint,1,1) = h_r[3]; 
  //Aij_ML(iPoint,1,2) = h_new[0]; 
  //Aij_ML(iPoint,2,2) = h_new[3]; 
  //exit(1);
  //}

  su2double delta_k = delta_sdd[5];

  if (dist<1e-10) {
    TKE_ML(iPoint) = 0.0;
  }
  else {
    TKE_ML(iPoint) = pow(10.0,delta_k)*turb_ke;
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
    delete [] h_r;
    delete [] h_orig;
    delete [] h_new;
}

void CTurbSSTVariable::EigenDecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val, unsigned short n){
  int iDim,jDim;
  su2double tmp, detV;
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
  su2double tmp;
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

}

void CTurbSSTVariable::q_from_mat(su2double **Vorig, su2double *q)  {
  unsigned short i;
  su2double t, tmp;

  su2double **V = new su2double* [3];
   for (i = 0; i < 3; i++){
    V[i] = new su2double [3];
    V[i][0] = +Vorig[i][0];
    V[i][1] = -Vorig[i][1];
    V[i][2] = Vorig[i][2];
  }

  if (V[2][2] < 0){
    if (V[0][0] > V[1][1]){
      t = 1.0 + V[0][0] - V[1][1] - V[2][2];
      q[0] = V[1][2]-V[2][1];
      q[1] = t;
      q[2] = V[0][1]+V[1][0];
      q[3] = V[2][0]+V[0][2];
    }
    else {
      t = 1.0 - V[0][0] + V[1][1] - V[2][2];
      q[0] = V[2][0]-V[0][2];
      q[1] = V[0][1]+V[1][0];
      q[2] = t; 
      q[3] = V[1][2]+V[2][1];
    }
  }
  else {
    if (V[0][0] < -V[1][1]){
      t = 1.0 - V[0][0] - V[1][1] + V[2][2];
      q[0] = V[0][1]-V[1][0];
      q[1] = V[2][0]+V[0][2];
      q[2] = V[1][2]+V[2][1];
      q[3] = t;
    }
    else {
      t = 1.0 + V[0][0] + V[1][1] + V[2][2];
      q[0] = t;
      q[1] = V[1][2]-V[2][1];
      q[2] = V[2][0]-V[0][2];
      q[3] = V[0][1]-V[1][0];
    }
  }
  for (unsigned short i = 0; i < 4; i++) { q[i] = 0.5*q[i] / sqrt(t); }
}

void CTurbSSTVariable::q_prod(su2double *q1, su2double *q2, su2double *qp)  {
  qp[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
  qp[1] = q1[0]*q2[1] + q2[0]*q1[1] - q1[2]*q2[3] + q2[2]*q1[3];
  qp[2] = q1[0]*q2[2] + q2[0]*q1[2] + q1[1]*q2[3] - q2[1]*q1[3];
  qp[3] = q1[0]*q2[3] + q2[0]*q1[3] - q1[1]*q2[2] + q2[1]*q1[2];
}

void CTurbSSTVariable::mat_from_q(su2double *q, su2double **V)  {
  su2double tmp1, tmp2;

  V[0][0] =  pow(q[1],2.0) - pow(q[2],2.0) - pow(q[3],2.0) + pow(q[0],2.0);
  V[1][1] = -pow(q[1],2.0) + pow(q[2],2.0) - pow(q[3],2.0) + pow(q[0],2.0);
  V[2][2] = -pow(q[1],2.0) - pow(q[2],2.0) + pow(q[3],2.0) + pow(q[0],2.0);
  
  tmp1 = q[1]*q[2];
  tmp2 = q[3]*q[0];
  V[0][1] = 2.0 * (tmp1 + tmp2);
  V[1][0] = 2.0 * (tmp1 - tmp2);
  tmp1 = q[1]*q[3];
  tmp2 = q[2]*q[0];
  V[0][2] = 2.0 * (tmp1 - tmp2);
  V[2][0] = 2.0 * (tmp1 + tmp2);
  tmp1 = q[2]*q[3];
  tmp2 = q[1]*q[0];
  V[1][2] = 2.0 * (tmp1 + tmp2);
  V[2][1] = 2.0 * (tmp1 - tmp2);
//  tmp1 = q[1]*q[2];
//  tmp2 = q[3]*q[0];
//  V[0][1] = 2.0 * (tmp1 - tmp2);
//  V[1][0] = 2.0 * (tmp1 + tmp2);
//  tmp1 = q[1]*q[3];
//  tmp2 = q[2]*q[0];
//  V[0][2] = 2.0 * (tmp1 + tmp2);
//  V[2][0] = 2.0 * (tmp1 - tmp2);
//  tmp1 = q[2]*q[3];
//  tmp2 = q[1]*q[0];
//  V[1][2] = 2.0 * (tmp1 - tmp2);
//  V[2][1] = 2.0 * (tmp1 + tmp2);

}

void CTurbSSTVariable::q_norm(su2double *q)  {
  su2double qnorm = sqrt(pow(q[0],2.0) + pow(q[1],2.0) + pow(q[2],2.0) + pow(q[3],2.0));
  for (unsigned short i = 0; i < 4; i++) { q[i] = q[i] / qnorm; }
  if (q[0] < 0) {
    for (unsigned short i = 0; i < 4; i++) { q[i] = -q[i]; }
  }
}

void CTurbSSTVariable::EigenSolve(su2double **A, su2double **Eig_Vec, su2double *Eig_Val){
  unsigned short i,j;

  su2double **As = new su2double* [3];
   for (i = 0; i < 3; i++){
    As[i] = new su2double [3];
     for (j = 0; j < 3; j++){
       As[i][j] = A[i][j];
     }
  }

  // Precondition the matrix by factoring out the maximum absolute value
  su2double max0 = max ( fabs(As[0][0]), fabs(As[0][1]) );
  su2double max1 = max ( fabs(As[0][2]), fabs(As[1][1]) );
  su2double max2 = max ( fabs(As[1][2]), fabs(As[2][2]) );
  su2double maxAbsElement = max(max(max0, max1), max2 );
  if ( maxAbsElement == 0.0 ) {
  // As is a zero matrix
  Eig_Val[0] = 0.0;
  Eig_Val[1] = 0.0;
  Eig_Val[2] = 0.0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      Eig_Vec[i][j] = 0.0;
      if (i==j){Eig_Vec[i][j] = 1.0;}
    }
  }
  return;
  }

  su2double invMaxAbsElement = 1.0 / maxAbsElement;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      As[i][j] *=  invMaxAbsElement;
    }
  }
  
  su2double norm = As[0][1]*As[0][1] + As[0][2] * As[0][2] + As[1][2] * As[1][2];
  
  if ( norm > 0.0 ){
    // Compute the eigenvalues of A.
    su2double q = ( As[0][0] + As[1][1] + As[2][2] ) / 3.0;
    su2double b00 = As[0][0] - q;
    su2double b11 = As[1][1] - q;
    su2double b22 = As[2][2] - q;

    su2double p = sqrt( (b00*b00 + b11*b11 + b22*b22 + norm*2.0) / 6.0 );
    su2double c00 = b11 * b22 - A[1][2] * A[1][2];
    su2double c01 = As[0][1] * b22 - As[1][2] * As[0][2];
    su2double c02 = As[0][1] * As[1][2] - b11 * As[0][2];
    su2double det = ( b00 * c00 - As[0][1] * c01 + As[0][2] * c02 ) / pow(p,3.0);

    su2double halfDet = det/2.0;
    halfDet = min(max(halfDet,-1.0),1.0);

    su2double angle = acos(halfDet) / 3.0;
    su2double twoThirdsPi = 2.09439510239319549;
    su2double beta2 = cos(angle)*2;
    su2double beta0 = cos(angle + twoThirdsPi)*2.0;
    su2double beta1 = -(beta0 + beta2);

    Eig_Val[0] = q + p * beta0;
    Eig_Val[1] = q + p * beta1;
    Eig_Val[2] = q + p * beta2;

    if (halfDet >= 0.0 ) {
      ComputeEigenvector0(As, Eig_Val[2], Eig_Vec[2] );
      ComputeEigenvector1(As, Eig_Vec[2], Eig_Val[1], Eig_Vec[1]);
      Eig_Vec[0] = Cross(Eig_Vec[1], Eig_Vec[2]);
      cout << "opt 1" << endl; 
//      for (j = 0; j < 3; j++) {
//        Eig_Vec[j][1] *= -1.0;
//        Eig_Vec[j][2] *= -1.0;
//      }
    }
    else {
      ComputeEigenvector0(As, Eig_Val[0], Eig_Vec[0] );
      ComputeEigenvector1(As, Eig_Vec[0], Eig_Val[1], Eig_Vec[1]);
      Eig_Vec[2] = Cross(Eig_Vec[0], Eig_Vec[1]);
//      for (j = 0; j < 3; j++) {
//        Eig_Vec[j][1] *= -1.0;
//        Eig_Vec[j][2] *= -1.0;
//      }
    }
  }
  else {
  // The matrix is diagonal.
    Eig_Val[0] = As[0][0];
    Eig_Val[1] = As[1][1];
    Eig_Val[2] = As[2][2];
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        Eig_Vec[i][j] = 0.0;
        if (i==j){Eig_Vec[i][j] = 1.0;}
      }
    }
  }

  // Revert the scaling.
  Eig_Val[0] *= maxAbsElement;
  Eig_Val[1] *= maxAbsElement;
  Eig_Val[2] *= maxAbsElement;
}

su2double* CTurbSSTVariable::Cross(su2double U[3], su2double V[3]) {
  su2double *C  = new su2double [3];
  C[0] = U[1]*V[2] - U[2]*V[1];
  C[1] = U[2]*V[0] - U[0]*V[2];
  C[2] = U[0]*V[1] - U[1]*V[0];
  return C;
}

void CTurbSSTVariable::ComputeOrthogonalComplement(su2double *W, su2double *U, su2double *V) {

  // Robustly compute a right−handed orthonormal set{U,V,W}. 
  su2double invLength;
  if (fabs(W[0]) > fabs(W[1])){
    invLength = 1.0 / sqrt(W[0]*W[0] + W[2]*W[2]);
    U[0] = -W[2]*invLength;
    U[1] = 0.0;
    U[2] = W[0]*invLength;
  }
  else {
  invLength = 1.0 / sqrt(W[1]*W[1] + W[2]*W[2]);
  U[0] = 0.0;
  U[1] = W[2]*invLength;
  U[2] = -W[1]*invLength;
  }
  V = Cross(W,U);
}

void CTurbSSTVariable::ComputeEigenvector0(su2double **A,  su2double eval0, su2double *evec0) {
  su2double row0[3] = { A[0][0] - eval0, A[0][1] , A[0][2] };
  su2double row1[3] = { A[0][1], A[1][1] - eval0 , A[1][2] };
  su2double row2[3] = { A[0][2], A[1][2] , A[2][2] - eval0 };
  su2double *r0xr1 = Cross(row0, row1);
  su2double *r0xr2 = Cross(row0, row2);
  su2double *r1xr2 = Cross(row1, row2);
  su2double d0 = Dot(r0xr1, r0xr1);
  su2double d1 = Dot(r0xr2, r0xr2);
  su2double d2 = Dot(r1xr2, r1xr2);
  su2double dmax = d0;
  
  unsigned int imax = 0;
  if (d1>dmax) {
    dmax = d1;
    imax = 1;
  }
  if (d2>dmax) {
    imax = 2;
  }
  if (imax == 0) {
    evec0 = Divide(r0xr1, sqrt(d0));
  }
  else if (imax == 1) {
    evec0 = Divide(r0xr2, sqrt(d1));
  }
  else {
    evec0 = Divide(r1xr2, sqrt(d2));
  }
}

void CTurbSSTVariable::ComputeEigenvector1(su2double **A, su2double *evec0, su2double eval1, su2double *evec1) {
  su2double *U  = new su2double [3];
  su2double *V  = new su2double [3];
  su2double *AU = new su2double [3];
  su2double *AV = new su2double [3];

  ComputeOrthogonalComplement(evec0, U, V);

  AU[0] = A[0][0]*U[0] + A[0][1]*U[1] + A[0][2]*U[2];
  AU[1] = A[0][1]*U[0] + A[1][1]*U[1] + A[1][2]*U[2];
  AU[2] = A[0][2]*U[0] + A[1][2]*U[1] + A[2][2]*U[2];

  AV[0] = A[0][0]*V[0] + A[0][1]*V[1] + A[0][2]*V[2];
  AV[1] = A[0][1]*V[0] + A[1][1]*V[1] + A[1][2]*V[2];
  AV[2] = A[0][2]*V[0] + A[1][2]*V[1] + A[2][2]*V[2];

  su2double m00 = U[0]*AU[0] + U[1]*AU[1] + U[2]*AU[2] - eval1;
  su2double m01 = U[0]*AV[0] + U[1]*AV[1] + U[2]*AV[2];
  su2double m11 = V[0]*AV[0] + V[1]*AV[1] + V[2]*AV[2] - eval1;

  su2double absM00 = fabs(m00);
  su2double absM01 = fabs(m01);
  su2double absM11 = fabs(m11);
  su2double maxAbsComp;
  if ( absM00 >= absM11 ) {
    maxAbsComp = max(absM00, absM01);
    if (maxAbsComp > 0.0 ) {
      if ( absM00 >= absM01 ) {
        m01 /= m00;
        m00 = 1.0 / sqrt(1.0 + m01*m01);
        m01 *= m00;
      }
      else {
        m00 /= m01;
        m01 = 1.0 / sqrt(1.0 + m00*m00);
        m00 *= m01;
      }
      evec1 = Subtract(Multiply(m01 , U), Multiply(m00, V));
    }
    else {
      evec1 = U;
    }
  }

  else {
    maxAbsComp = max(absM11, absM01);
    if (maxAbsComp > 0.0) {
      if (absM11 >= absM01) {
        m01 /= m11;
        m11 = 1.0/sqrt(1.0 + m01*m01);
        m01 *= m11;
      }
      else {
        m11 /= m01;
        m01 = 1.0 / sqrt(1.0 + m11*m11);
        m11 *= m01;
      }
      evec1 = Subtract(Multiply(m11 , U), Multiply(m01, V));
    }
    else {
      evec1 = U;
    }
  }
}

su2double* CTurbSSTVariable::Multiply(su2double s, su2double *U) {
  su2double *product  = new su2double [3];
  for (unsigned int i = 0; i < 3; i++) {product[i] = s*U[i];}
  return product;
}

su2double* CTurbSSTVariable::Subtract(su2double *U, su2double *V) {
  su2double *diff  = new su2double [3];
  for (unsigned int i = 0; i < 3; i++) {diff[i] = U[i] - V[i];}
  return diff;
}

su2double* CTurbSSTVariable::Divide(su2double *U, su2double s) {
  su2double *div  = new su2double [3];
  su2double invS = 1.0/s;
  for (unsigned int i = 0; i < 3; i++) {div[i] = invS*U[i];}
  return div;
}

su2double CTurbSSTVariable::Dot(su2double *U, su2double *V) {
  su2double dot = U[0]*V[0] + U[1]*V[1] + U[2]*V[2];
  return dot;
}

unsigned int* CTurbSSTVariable::EigenSort(su2double *v) {
  unsigned int i, max_i, min_i;
  unsigned int *idx = new unsigned int[3];
  su2double min_v, max_v;
  if (v[0]!=v[1] && v[0]!=v[2]) {
    max_v = -99.0;
    min_v = 99.0;
    for (unsigned int i = 0; i < 3; i++) {
      if (v[i] > max_v) {
        max_v = v[i];
        max_i = i;
      }
      if (v[i] < min_v) {
        min_v = v[i];
        min_i = i;
      }
    }
    idx[0] = min_i;
    idx[2] = max_i;
    idx[1] = 3 - (min_i + max_i);  
  }
  else {
    idx[0] = 0;
    idx[1] = 1;
    idx[2] = 2;
  }
  return idx;
}
