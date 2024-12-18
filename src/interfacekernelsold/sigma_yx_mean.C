//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "sigma_yx_mean.h"

registerMooseObject("MooseApp", sigma_yx_mean);

InputParameters
sigma_yx_mean::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("mu", "Forward reaction rate coefficient.");
  //params.addRequiredParam<Real>("eps", "Backward reaction rate coefficient.");
  params.addRequiredParam<FunctionName>("eps", "A function that describes the body force");
  params.addRequiredCoupledVar("disp_x", "Displacement along y");
  return params;
}

sigma_yx_mean::sigma_yx_mean(const InputParameters & parameters) :
InterfaceKernel(parameters),
_kf(getParam<Real>("mu")),
//_eps(getParam<Real>("eps"))
_eps(getFunction("eps")),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_x_n(coupledNeighborGradient("disp_x"))
{
}

Real
sigma_yx_mean::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);

  Real j_u=_neighbor_value[_qp]-_u[_qp];
  Real j_test;

  RealVectorValue m_grad_disp_x=0.5*(_grad_disp_x[_qp]+_grad_disp_x_n[_qp]);

  switch (type)
  {
    case Moose::Element:
      j_test = -_test[_i][_qp];
      break;
    case Moose::Neighbor:
      j_test = _test_neighbor[_i][_qp] ;
      break;
  }
  r = (0.5*m_grad_disp_x(1)+0.5/eps*j_u)*j_test;
  r = 2.0*_kf*r;
  return r;
}

Real
sigma_yx_mean::computeQpJacobian(Moose::DGJacobianType type)
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
