#pragma once

#define HAVE_INT8_T
#include "Mesh.h"
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoCone.h>

class Painter
{
public:
	SoSeparator* getShapeSep(Mesh* mesh);
	SoSeparator* getThickLineSep(Mesh* mesh, int edgeIndex);
	SoSeparator* getSpheresSep(Mesh* mesh, float deltaX, float deltaY, float scale);
};
