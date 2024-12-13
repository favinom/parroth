#include"myHelpers.h"

void fillRowOfTensor(RealVectorValue const & r,int const comp,RealTensorValue & V)
{
	for (int i=0; i<3; ++i)
		V(comp,i)=r(i);
}
