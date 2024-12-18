//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElasticityInterfaceKernel.h"

#include "myHelpers.h"

registerMooseObject("MooseApp", ElasticityInterfaceKernel);


InputParameters
ElasticityInterfaceKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("mu", "Forward reaction rate coefficient.");
  params.addRequiredParam<Real>("lambda", "Lame first parameter.");  //cambio
  params.addRequiredParam<int>("component", "Lame first parameter.");  //cambio
  params.addRequiredParam<bool>("jump", "Lame first parameter.");  //cambio
  params.addRequiredParam<FunctionName>("eps", "A function that describes the body force");
  //params.addRequiredParam<Real>("eps", "Backward reaction rate coefficient.");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z", "Displacement along z");
  return params;
}

ElasticityInterfaceKernel::ElasticityInterfaceKernel(const InputParameters & parameters) :
InterfaceKernel(parameters),
_dim(_mesh.dimension()),
_component(getParam<int>("component")),
_isJump(getParam<bool>("jump")),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda")),   //cambio
//_eps(getParam<Real>("eps"))
_eps(getFunction("eps")),
_disp_x(coupledValue("disp_x")),
_disp_y(coupledValue("disp_y")),
_disp_z(_dim>2 ? coupledValue("disp_z") : _zero),
_disp_x_n(coupledNeighborValue("disp_x")),
_disp_y_n(coupledNeighborValue("disp_y")),
_disp_z_n(_dim>2 ? coupledNeighborValue("disp_z") : _zero),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(coupledGradient("disp_y")),
_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
_grad_disp_x_n(coupledNeighborGradient("disp_x")),
_grad_disp_y_n(coupledNeighborGradient("disp_y")),
_grad_disp_z_n(_dim>2 ? coupledNeighborGradient("disp_z") : _grad_zero)
{
}

Real
ElasticityInterfaceKernel::computeQpResidual(Moose::DGResidualType type)
{
  RealTensorValue U;
  RealTensorValue U_n;
  RealTensorValue U_a;
  RealVectorValue u;
  RealVectorValue u_n;
  RealVectorValue u_j;

  RealTensorValue U_i;
  RealTensorValue I(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);

  U.zero();
  U_n.zero();
  U_a.zero();
  u.zero();
  u_n.zero();
  u_j.zero();

  U_i.zero();

  fillRowOfTensor(_grad_disp_x[_qp],0,U);
  fillRowOfTensor(_grad_disp_y[_qp],1,U);
  fillRowOfTensor(_grad_disp_z[_qp],2,U);

  fillRowOfTensor(_grad_disp_x_n[_qp],0,U_n);
  fillRowOfTensor(_grad_disp_y_n[_qp],1,U_n);
  fillRowOfTensor(_grad_disp_z_n[_qp],2,U_n);

  u(0)=_disp_x[_qp];
  u(1)=_disp_y[_qp];
  u(2)=_disp_z[_qp];

  u_n(0)=_disp_x_n[_qp];
  u_n(1)=_disp_y_n[_qp];
  u_n(2)=_disp_z_n[_qp];

  Real eps=_eps.value(_t,_q_point[_qp]);

  U_a=eps*0.5*(U+U_n);
  u_j=u_n-u;

  RealTensorValue P(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
  RealTensorValue R;

  for (int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j)
    {
      R(i,j)=_normals[_qp](i)*_normals[_qp](j);
    }
  }
  P=P-R;

  RealTensorValue uj_o_n;

  ///// Controlar la normal no olvida

  for (int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j)
    {
      uj_o_n(i,j)=u_j(i)*_normals[_qp](j);
    }
  }

  U_i=U_a*P+uj_o_n;

  RealTensorValue Ei=0.5*(U_i+U_i.transpose());
  RealTensorValue Si=2.0*_mu*Ei+_lambda*Ei.tr()*I;

  RealTensorValue ST=Si*P;
  RealVectorValue Sn=Si*_normals[_qp];

  if (_isJump)
  {
    RealVectorValue v;
    RealTensorValue V;

    switch (type)
    {
    case Moose::Element:
      v = 0.5*_grad_test[_i][_qp];
      break;
    case Moose::Neighbor:
      v = 0.5*_grad_test_neighbor[_i][_qp];
      break;
    }

    fillRowOfTensor(v,_component,V);
    V=V*P;
    return ST.contract(V);
  }
  else
  {
    Real v;
    switch (type)
    {
    case Moose::Element:
      v = -_test[_i][_qp];
      break;
    case Moose::Neighbor:
      v = _test_neighbor[_i][_qp];
      break;
    }
    return Sn(_component)*v/eps;
  }

// RealVectorValue grad_u=P*_grad_u[_qp];
// RealVectorValue grad_u_n=P*_grad_neighbor_value[_qp];

// RealVectorValue grad_test=P*_grad_test[_i][_qp];
// RealVectorValue grad_test_n=P*_grad_test_neighbor[_i][_qp];

// RealVectorValue m_grad_u=0.5*(grad_u+grad_u_n);
// Real j_disp_y=_disp_y_n[_qp]-_disp_y[_qp];

//r=eps*m_grad_u*m_grad_test+j_disp_y*m_grad_test(1);
//r=_mu*r;

//return r;
}

Real
ElasticityInterfaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  return 0.0;
  // Real jac = 0;
  // Real eps=_eps.value(_t,_q_point[_qp]);

  // RealTensorValue P(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
  // RealTensorValue R;
  // for (int i=0; i<3; ++i)
  //   for(int j=0; j<3; ++j)
  //     R(i,j)=_normals[_qp](i)*_normals[_qp](j);

  // P=P-R;
  // //std::cout<<P<<std::endl;

  // RealVectorValue grad_phi=P*_grad_phi[_j][_qp];
  // RealVectorValue grad_phi_n=P*_grad_phi_neighbor[_j][_qp];

  // RealVectorValue grad_test=P*_grad_test[_i][_qp];
  // RealVectorValue grad_test_n=P*_grad_test_neighbor[_i][_qp];


  // switch (type)
  // {
  //   case Moose::ElementElement:
  //     jac = _mu*eps/4.0*(  grad_phi ) * grad_test;
  //     break;
  //   case Moose::NeighborNeighbor:
  //     jac = _mu*eps/4.0*(  grad_phi_n  ) * grad_test_n;
  //     break;
  //   case Moose::NeighborElement:
  //     jac = _mu*eps/4.0*(    grad_phi ) * grad_test_n;
  //     break;
  //   case Moose::ElementNeighbor:
  //     jac = _mu*eps/4.0*(  grad_phi_n ) * grad_test;
  //     break;
  // }
  // return jac;
}
