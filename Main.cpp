#define HAVE_INT8_T
#include "Painter.h"
#include "Segmentor.h"
#include "KMeans.h"
#include "GraphOperations.h"
#include "DijstraImplementations.h"
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <random>

//#include <Eigen/Dense>

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
	mesh->loadOff("for fprinting/man0.off");
	//mesh->loadOff("0.off");
	//mesh->loadOff("car.off");
	//mesh->loadOff("coffeecup.off");
	//mesh->loadOff("bunny.off");
	//mesh->loadOff("cube3.off");
	//mesh->loadOff("faces/face.off");
	//mesh->loadOff("faces/face-low.off");
	//mesh->loadOff("faces/facem.off");
	//mesh->loadOff("faces/facem-low.off");

	//mesh->loadOff("doubleOpenCube3.off");
	//mesh->createCube(20.0f);
	//mesh->createOpenCube(20.0f);
	//mesh->createDoubleOpenCube(20.0f);

	//GraphOperations::generateDiskParameterization(mesh, GraphOperations::ParameterizationMethod::UNIFORM);
	//Segmentor segmentor(mesh);
	KMeans kmeans(mesh);
	//segmentor.assignDiameterValuesOfVertices();
	//segmentor.setColorValuesToVertices();
	//kmeans.assignClusterIdsOfVertices();
	kmeans.setColorValuesToVertices();


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
	std::srand(std::time(NULL)); //seeding for the first time only!
	const int min = 0;
	const int max = mesh->verts.size();
	const int range = max - min + 1;
	const int num1 = std::rand() % range + min;
	const int num2 = std::rand() % range + min;
	const int source_vert_Id = 213;
	const int dest_vert_Id = 451;
	mesh->samples = { mesh->verts[source_vert_Id]->idx , mesh->verts[dest_vert_Id]->idx };
	Graph gmesh = Dijkstra::meshToGraph(mesh);
	Dijkstra::timing(gmesh, source_vert_Id, dest_vert_Id, Dijkstra::ARRAY);
	//Dijkstra::fprinting(gmesh, Dijkstra::ARRAY);
	//Dijkstra::fprintingOnce(gmesh,source_vert_Id, dest_vert_Id, Dijkstra::ARRAY);
	Dijkstra::pathDrawing(mesh, gmesh, source_vert_Id, dest_vert_Id, Dijkstra::MIN_HEAP);

	for (size_t i = 0; i < mesh->edges.size(); i++)
	{
		//if (mesh->edges[i]->isItBoundary)
		//if (mesh->edges[i]->isItInLongestBoundary)
		//if (mesh->edges[i]->isItPathPart)
		//if (mesh->edges[i]->isItTraversed)
		if (mesh->edges[i]->isInShortestPath)
			root->addChild(painter->getThickLineSep(mesh, i));
	}

	root->addChild(painter->getSpheresSep(mesh, 0.f, 0.f, 1.f));

	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
