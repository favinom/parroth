//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "sigma_xx_mean_comp.h"

registerMooseObject("MooseApp", sigma_xx_mean_comp);

InputParameters
sigma_xx_mean_comp::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("mu", "Forward reaction rate coefficient.");
  params.addRequiredParam<Real>("lambda", "Lame first parameter.");  //cambio
  //params.addRequiredParam<Real>("eps", "Backward reaction rate coefficient.");
  params.addRequiredParam<FunctionName>("eps", "A function that describes the body force");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  return params;
}

sigma_xx_mean_comp::sigma_xx_mean_comp(const InputParameters & parameters) :
InterfaceKernel(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda")),  //cambio
//_eps(getParam<Real>("eps"))
_eps(getFunction("eps")),
_grad_disp_y(coupledGradient("disp_y")),
_grad_disp_y_n(coupledNeighborGradient("disp_y"))
{
}

Real
sigma_xx_mean_comp::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);

  Real j_u=_neighbor_value[_qp]-_u[_qp];
  Real j_test;

  RealVectorValue m_grad_disp_y=0.5*(_grad_disp_y[_qp]+_grad_disp_y_n[_qp]);
  
  switch (type)
  {
    case Moose::Element:
      j_test = -_test[_i][_qp];
      break;
    case Moose::Neighbor:
      j_test = _test_neighbor[_i][_qp] ;
      break;
  }
  r = (1.0/eps*(2.0*_mu+_lambda)*j_u+_lambda*m_grad_disp_y(1))*j_test;
  return r;
}

Real
sigma_xx_mean_comp::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);
  switch (type)
  {
    case Moose::ElementElement:
      jac = (2*_mu+_lambda)/(2.0*eps)*  _phi[_j][_qp]  * _test[_i][_qp];
      break;
    case Moose::NeighborNeighbor:
      jac = (2*_mu+_lambda)/(2.0*eps)* _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp] ;
      break;
    case Moose::NeighborElement:
      jac = -(2*_mu+_lambda)/(2.0*eps)* _phi[_j][_qp] * _test_neighbor[_i][_qp] ;
      break;
    case Moose::ElementNeighbor:
      jac = -(2*_mu+_lambda)/(2.0*eps)* _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;
  }
  return jac;
}
