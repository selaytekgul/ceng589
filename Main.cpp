#define HAVE_INT8_T
#include "Painter.h"
#include "Segmentor.h"
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer* viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator* root = new SoSeparator;
	root->ref();

	//stuff to be drawn on screen must be added to the root
//	SoCone * cone = new SoCone;
//	root->addChild(cone);

	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();
	Segmentor* segmentor = new Segmentor();

	mesh->loadOff("0.off");
	//mesh->loadOff("car.off");
	//mesh->loadOff("coffeecup.off");
	//mesh->loadOff("bunny.off");
	//mesh->loadOff("cube3.off");
//	mesh->createCube(20.0f);

	segmentor->assignLengthValuesOfVertices(mesh);
	segmentor->setColorValuesToVertices(mesh);

	//cout << "my (verts[4]) 1-ring neighborhood is: \n";
	//for (int nv = 0; nv < mesh->verts[4]->vertList.size(); nv++)
	//	cout << mesh->verts[4]->vertList[nv] << " neighbb\n";


	//cout << "my (verts[4]) 1-ring neighborhood is: \n";
	//for (int ne = 0; ne < mesh->verts[4]->edgeList.size(); ne++)
	//	if (mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v1i == 4)
	//		cout << mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v2i << " nnnnnnnnn\n";
	//	else
	//		cout << mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v1i << " nnnnnnnnn\n";

//		cout << mesh->verts[4]->vertList[nv] << " neighbb\n";

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
