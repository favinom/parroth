//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"
#include "Function.h"

// Forward Declarations
class ElasticityInterfaceKernel;

template <>
InputParameters validParams<ElasticityInterfaceKernel>();

/**
 * Implements a reaction to establish ReactionRate=k_f*u-k_b*v
 * at interface.
 */
class ElasticityInterfaceKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  ElasticityInterfaceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  int const _dim;
  int const _component;
  bool const _isJump;
  /// Forward reaction rate coefficient
  Real _mu;
  Real _lambda;

  /// Backward reaction rate coefficient
  //Real _eps;
  Function const & _eps;


  VariableValue const & _disp_x;
  VariableValue const & _disp_y;
  VariableValue const & _disp_z;
  VariableValue const & _disp_x_n;
  VariableValue const & _disp_y_n;
  VariableValue const & _disp_z_n;

  VariableGradient const & _grad_disp_x;
  VariableGradient const & _grad_disp_y;
  VariableGradient const & _grad_disp_z;
  VariableGradient const & _grad_disp_x_n;
  VariableGradient const & _grad_disp_y_n;
  VariableGradient const & _grad_disp_z_n;

};