#include"myHelpers.h"

void fillRowOfTensor(RealVectorValue const & r,int const comp,RealTensorValue & V)
{
	for (int i=0; i<3; ++i)
		V(comp,i)=r(i);
}


void computeStress(RealTensorValue const & U,double const mu,double const lambda,RealTensorValue & S)
{
	RealTensorValue E=0.5*(U+U.transpose());
	double trE=E.tr();
	RealTensorValue I(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	S=2.0*mu*E+lambda*trE*I;
}
