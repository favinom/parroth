//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "sigma_yx_jump_comp.h"

registerMooseObject("MooseApp", sigma_yx_jump_comp);

//defineLegacyParams(MyInterfaceDiffusion);

InputParameters
sigma_yx_jump_comp::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("mu", "Forward reaction rate coefficient.");
  params.addRequiredParam<Real>("lambda", "Lame first parameter.");  //cambio
  //params.addRequiredParam<Real>("eps", "Backward reaction rate coefficient.");
  params.addRequiredParam<FunctionName>("eps", "A function that describes the body force");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  return params;
}

sigma_yx_jump_comp::sigma_yx_jump_comp(const InputParameters & parameters) :
InterfaceKernel(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda")),  //cambio
//_eps(getParam<Real>("eps"))
_eps(getFunction("eps")),
_disp_x(coupledValue("disp_x")),
_disp_x_n(coupledNeighborValue("disp_x"))
{
}

Real
sigma_yx_jump_comp::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);

  RealTensorValue P(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
  RealTensorValue R;
  for (int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      R(i,j)=_normals[_qp](i)*_normals[_qp](j);

  P=P-R;
  //std::cout<<P<<std::endl;

  RealVectorValue grad_u=P*_grad_u[_qp];
  RealVectorValue grad_u_n=P*_grad_neighbor_value[_qp];

  RealVectorValue grad_test=P*_grad_test[_i][_qp];
  RealVectorValue grad_test_n=P*_grad_test_neighbor[_i][_qp];

  RealVectorValue m_grad_u=0.5*(grad_u+grad_u_n);
  RealVectorValue m_grad_test;
  Real j_disp_x=_disp_x_n[_qp]-_disp_x[_qp];

  switch (type)
  {
  case Moose::Element:
    m_grad_test=0.5*grad_test;
    break;
  case Moose::Neighbor:
    m_grad_test=0.5*grad_test_n;
    break;
  }
  r = (eps*(2.0*_mu+_lambda)*m_grad_u(1)+_lambda*j_disp_x)*m_grad_test(1);
  return r;
}

Real
sigma_yx_jump_comp::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  Real eps=_eps.value(_t,_q_point[_qp]);

  RealTensorValue P(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
  RealTensorValue R;
  for (int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      R(i,j)=_normals[_qp](i)*_normals[_qp](j);

  P=P-R;
  //std::cout<<P<<std::endl;

  RealVectorValue grad_phi=P*_grad_phi[_j][_qp];
  RealVectorValue grad_phi_n=P*_grad_phi_neighbor[_j][_qp];

  RealVectorValue grad_test=P*_grad_test[_i][_qp];
  RealVectorValue grad_test_n=P*_grad_test_neighbor[_i][_qp];


  switch (type)
  {
    case Moose::ElementElement:
      jac = (2.0*_mu+_lambda)*eps/2.0*(  grad_phi ) * grad_test;
      break;
    case Moose::NeighborNeighbor:
      jac = (2.0*_mu+_lambda)*eps/2.0*(  grad_phi_n  ) * grad_test_n;
      break;
    case Moose::NeighborElement:
      jac = (2.0*_mu+_lambda)*eps/2.0*(    grad_phi ) * grad_test_n;
      break;
    case Moose::ElementNeighbor:
      jac = (2.0*_mu+_lambda)*eps/2.0*(  grad_phi_n ) * grad_test;
      break;
  }
  return jac;
}
