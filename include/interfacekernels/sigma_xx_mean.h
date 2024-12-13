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
class sigma_xx_mean;

template <>
InputParameters validParams<sigma_xx_mean>();

/**
 * Implements a reaction to establish ReactionRate=k_f*u-k_b*v
 * at interface.
 */
class sigma_xx_mean : public InterfaceKernel
{
public:
  static InputParameters validParams();

  sigma_xx_mean(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Forward reaction rate coefficient
  Real _kf;

  /// Backward reaction rate coefficient
  //Real _eps;
  Function const & _eps;
};
