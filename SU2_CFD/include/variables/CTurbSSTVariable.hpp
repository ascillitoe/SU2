/*!
 * \file CTurbSSTVariable.hpp
 * \brief Declaration of the variables of the SST turbulence model.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "CTurbVariable.hpp"

/*!
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSSTVariable final : public CTurbVariable {
protected:
  su2double sigma_om2;
  su2double beta_star;
  VectorType F1;
  VectorType F2;    /*!< \brief Menter blending function for blending of k-w and k-eps. */
  VectorType CDkw;  /*!< \brief Cross-diffusion. */
  CVectorOfMatrix Aij_ML;
  VectorType TKE_ML;


public:
  /*!
   * \brief Constructor of the class.
   * \param[in] kine - Turbulence kinetic energy (k) (initialization value).
   * \param[in] omega - Turbulent variable value (initialization value).
   * \param[in] mut - Eddy viscosity (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants -
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint,
                   unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSSTVariable() override = default;

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   */
  void SetBlendingFunc(unsigned long iPoint, su2double val_viscosity, su2double val_dist, su2double val_density) override;

  void InitSDD(unsigned long iPoint, su2double muT, su2double turb_ke, su2double rho, su2double **PrimGrad, su2double *delta_sdd, su2double dist) override;

  /*!
   * \brief Get the first blending function.
   */
  inline su2double GetF1blending(unsigned long iPoint) const override { return F1(iPoint); }

  /*!
   * \brief Get the second blending function.
   */
  inline su2double GetF2blending(unsigned long iPoint) const override { return F2(iPoint); }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline su2double GetCrossDiff(unsigned long iPoint) const override { return CDkw(iPoint); }

  /*!
   * \brief Get the ML derived anisotropy tensor at iPoint
   * \return val_Aij_ML - Value of Aij_ML
   */
  inline su2double **GetAijML(unsigned long iPoint) override  { return Aij_ML[iPoint]; }

  /*!
   * \brief Get the ML derived turbulent kinetic energy at iPoint
   * \return val_Aij_ML - Value of TKE_ML
   */
  inline su2double GetTKEML(unsigned long iPoint) override  { return TKE_ML(iPoint); }

  /*!
   * \brief Decomposes the symmetric matrix A_ij, into eigenvectors and eigenvalues
   * \param[in] A_i: symmetric matrix to be decomposed
   * \param[in] Eig_Vec: strores the eigenvectors
   * \param[in] Eig_Val: stores the eigenvalues
   * \param[in] n: order of matrix A_ij
   */
  static void EigenDecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val, unsigned short n);

  /*!
   * \brief Recomposes the eigenvectors and eigenvalues into a matrix
   * \param[in] A_ij: recomposed matrix
   * \param[in] Eig_Vec: eigenvectors
   * \param[in] Eig_Val: eigenvalues
   * \param[in] n: order of matrix A_ij
   */
  static void EigenRecomposition(su2double **A_ij, su2double **Eig_Vec, const su2double *Eig_Val, unsigned short n);

  /*!
   * \brief tred2
   * \param[in] V: matrix that needs to be decomposed
   * \param[in] d: holds eigenvalues
   * \param[in] e: supplemental data structure
   * \param[in] n: order of matrix V
   */
  static void tred2(su2double **V, su2double *d, su2double *e, unsigned short n);

  /*!
   * \brief tql2
   * \param[in] V: matrix that will hold the eigenvectors
   * \param[in] d: array that will hold the ordered eigenvalues
   * \param[in] e: supplemental data structure
   * \param[in] n: order of matrix V
   */
  static void tql2(su2double **V, su2double *d, su2double *e, unsigned short n);

  /*!
   * \brief Calculate the quarternion from a rotation matrix
   * \param[in] V: matrix that will hold the rotation matrix
   * \param[in] q: array that will hold the calcuated quaternions
   */
  static void q_from_mat(su2double **V, su2double *q);

  /*!
   * \brief Calculates the product of two quarternions
   * \param[in] q1: The first quarternion
   * \param[in] q2: The second quarternion
   * \param[in] qp: The productg of q1 and q2
   */
  static void q_prod(su2double *q1, su2double *q2, su2double *qp);

  /*!
   * \brief Calculates the rotation matrix from a quarternion
   * \param[in] q: The quarternion
   * \param[in] V: The rotation matrix
   */
  static void mat_from_q(su2double *q, su2double **V);

  /*!
   * \brief Normalises a quarternions to the unit quarternion
   * \param[in] q: The quarternion
   */
  static void q_norm(su2double *q);

  static unsigned int *EigenSort(su2double *Eig_Val);

};
