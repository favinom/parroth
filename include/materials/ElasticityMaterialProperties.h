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
class ElasticityMaterialProperties : public Material
{
public:
  static InputParameters validParams();

  ElasticityMaterialProperties(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  int const _dim;

  VariableGradient const & _grad_disp_x;
  VariableGradient const & _grad_disp_y;
  VariableGradient const & _grad_disp_z;
  Real const _mu;
  Real const _lambda;

  MaterialProperty<Real> & _mu_prop;
  MaterialProperty<Real> & _lambda_prop;
  MaterialProperty<RealTensorValue> & _sigma_prop;

  RealTensorValue _U;

};
