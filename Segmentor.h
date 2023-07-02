#pragma once

#define HAVE_INT8_T
#include "Mesh.h"

class Segmentor
{
public:
	void assignLengthValuesOfVertices(Mesh* mesh);
	void setColorValuesToVertices(Mesh* mesh);
};
