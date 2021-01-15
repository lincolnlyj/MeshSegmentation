#include "StructsAndConsts.h"

Face::Face()
{
	State = UNDISCOVERED;
	Kind = A;
}

Vertice::Vertice()
{
	ColorData.R = 0;
	ColorData.G = 0;
	ColorData.B = 0;
}

bool operator < (const Edge & A, const Edge & B)
{
	return A.Distance > B.Distance;
}
