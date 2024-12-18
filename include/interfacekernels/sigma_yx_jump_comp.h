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
class sigma_yx_jump_comp;

template <>
InputParameters validParams<sigma_yx_jump_comp>();

/**
 * Implements a reaction to establish ReactionRate=k_f*u-k_b*v
 * at interface.
 */
class sigma_yx_jump_comp : public InterfaceKernel
{
public:
  static InputParameters validParams();

  sigma_yx_jump_comp(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Forward reaction rate coefficient
  Real _mu;
  Real _lambda;

  /// Backward reaction rate coefficient
  //Real _eps;
  Function const & _eps;

  VariableValue const & _disp_x;
  VariableValue const & _disp_x_n;
};
