#include "StructsAndConsts.h"

bool operator<(const Edge & A, const Edge & B)
{
	return A.Distance > B.Distance;
}

Face::Face()
{
	State = UNDISCOVERED;
	Kind = DETERMINED;
}

Color::Color(int r, int g, int b)
{
	R = r;
	G = g;
	B = b;
}

