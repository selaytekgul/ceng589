#define HAVE_INT8_T
#include "Painter.h"
#include "Segmentor.h"
#include "KMeans.h"
#include "GraphOperations.h"
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer* viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator* root = new SoSeparator;
	root->ref();

	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();

	//mesh->loadOff("tr_reg_000.off");
	//mesh->loadOff("for timing/centaur.off");
	//mesh->loadOff("for timing/man.off");
	//mesh->loadOff("for timing/weirdSphere.off");
	//mesh->loadOff("for fprinting/horse0.off");
	//mesh->loadOff("for fprinting/man0.off");
	//mesh->loadOff("0.off");
	//mesh->loadOff("car.off");
	//mesh->loadOff("coffeecup.off");
	//mesh->loadOff("bunny.off");
	//mesh->loadOff("cube3.off");
	mesh->createOpenCube(20.0f);
	GraphOperations::generateDiskParameterization(mesh, GraphOperations::ParameterizationMethod::HARMONIC);
	//Segmentor segmentor(mesh);
	KMeans kmeans(mesh);
	//segmentor.assignDiameterValuesOfVertices();
	//segmentor.setColorValuesToVertices();
	//kmeans.assignClusterIdsOfVertices();
	kmeans.setColorValuesToVertices();

	root->addChild( painter->getShapeSep(mesh) );


	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
