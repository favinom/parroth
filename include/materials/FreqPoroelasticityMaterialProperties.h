//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * This material uses a LinearInterpolation object to define the dependence
 * of the material's value on a variable.
 */
class FreqPoroelasticityMaterialProperties : public Material
{
public:
  static InputParameters validParams();

  FreqPoroelasticityMaterialProperties(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  int const _dim;

  VariableGradient const & _grad_disp_r_x;
  VariableGradient const & _grad_disp_i_x;
  VariableGradient const & _grad_disp_r_y;
  VariableGradient const & _grad_disp_i_y;
  VariableGradient const & _grad_disp_r_z;
  VariableGradient const & _grad_disp_i_z;
  VariableValue    const & _p_r;
  VariableValue    const & _p_i;

  Real const _mu;
  Real const _lambda;
  Real const _alpha;
  Real const _m;
  Real const _k;

  MaterialProperty<Real> & _mu_prop;
  MaterialProperty<Real> & _lambda_prop;
  MaterialProperty<Real> & _alpha_prop;
  MaterialProperty<Real> & _m_prop;
  MaterialProperty<Real> & _k_prop;

  MaterialProperty<RealTensorValue> & _eps_r_prop;
  MaterialProperty<RealTensorValue> & _eps_i_prop;
  MaterialProperty<RealTensorValue> & _sigma_r_prop;
  MaterialProperty<RealTensorValue> & _sigma_i_prop;

  RealTensorValue _U_r;
  RealTensorValue _U_i;
  RealTensorValue _I;

};
