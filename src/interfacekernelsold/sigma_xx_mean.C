//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "sigma_xx_mean.h"

registerMooseObject("MooseApp", sigma_xx_mean);

InputParameters
sigma_xx_mean::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("mu", "Forward reaction rate coefficient.");
  //params.addRequiredParam<Real>("eps", "Backward reaction rate coefficient.");
  params.addRequiredParam<FunctionName>("eps", "A function that describes the body force");
  return params;
}

sigma_xx_mean::sigma_xx_mean(const InputParameters & parameters) :
InterfaceKernel(parameters),
_kf(getParam<Real>("mu")),
//_eps(getParam<Real>("eps"))
_eps(getFunction("eps"))
{
}

Real
sigma_xx_mean::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);

  Real j_u=_neighbor_value[_qp]-_u[_qp];
  Real j_test;

  switch (type)
  {
    case Moose::Element:
      j_test = -_test[_i][_qp];
      break;
    case Moose::Neighbor:
      j_test = _test_neighbor[_i][_qp] ;
      break;
  }
  r = 1.0/eps*j_u*j_test;
  r = 2.0*_kf*r;
  return r;
}

Real
sigma_xx_mean::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);
  switch (type)
  {
    case Moose::ElementElement:
      jac = _kf/eps*  _phi[_j][_qp]  * _test[_i][_qp];
      break;
    case Moose::NeighborNeighbor:
      jac = _kf/eps* _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp] ;
      break;
    case Moose::NeighborElement:
      jac = -_kf/eps* _phi[_j][_qp] * _test_neighbor[_i][_qp] ;
      break;
    case Moose::ElementNeighbor:
      jac = -_kf/eps* _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;
  }
  return jac;
}
