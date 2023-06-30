#define HAVE_INT8_T
#include "Mesh.h"
#include "Painter.h"
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

/* a = b - c */
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];

#define innerProduct(v,q) \
((v)[0] * (q)[0] + \
	(v)[1] * (q)[1] + \
	(v)[2] * (q)[2])

#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

void crossProductFunction(int v_A[], int v_B[], int c_P[]) {
	c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
	c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

float rayIntersectsTriangle(float* p, float* d,
	float* v0, float* v1, float* v2) {

	float e1[3], e2[3], h[3], s[3], q[3];
	float a, f, u, v;
	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct(h, d, e2);
	a = innerProduct(e1, h);

	if (a > -0.00001 && a < 0.00001)
		//return(false);
		return(-1);

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));

	if (u < 0.0 || u > 1.0)
		//return(false);
		return(-1);

	crossProduct(q, s, e1);
	v = f * innerProduct(d, q);

	if (v < 0.0 || u + v > 1.0)
		//return(false);
		return(-1);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	float t = f * innerProduct(e2, q);

	if (t > 0.00001) // ray intersection
		//return(true);
		return(t);

	else // this means that there is a line intersection
		 // but not a ray intersection
		//return (false);
		return (-2);
}


float calculateLength(int v[]) {
	int square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	float root = sqrt(square);
	return root;
}



void normalize(int v[], float nor[]) {
	float len = calculateLength(v);
	nor[0] = v[0] / len;
	nor[1] = v[1] / len;
	nor[2] = v[2] / len;
}


void triangleNormalVectors(Mesh* mesh)
{
	for (size_t i = 0; i < mesh->tris.size(); i++)
	{
		int vertex1ID = mesh->tris[i]->v1i;
		int vertex2ID = mesh->tris[i]->v2i;
		int vertex3ID = mesh->tris[i]->v3i;
		mesh->verts[vertex1ID]->length = -5.0;
		mesh->verts[vertex2ID]->length = -5.0;
		mesh->verts[vertex3ID]->length = -5.0;
	}
	for (size_t i = 0; i < mesh->tris.size(); i++)
	{
		int vertex1ID = mesh->tris[i]->v1i;
		int vertex2ID = mesh->tris[i]->v2i;
		int vertex3ID = mesh->tris[i]->v3i;

		int x1_A = mesh->verts[vertex1ID]->coords[0];
		int y1_A = mesh->verts[vertex1ID]->coords[1];
		int z1_A = mesh->verts[vertex1ID]->coords[2];

		int x2_B = mesh->verts[vertex2ID]->coords[0];
		int y2_B = mesh->verts[vertex2ID]->coords[1];
		int z2_B = mesh->verts[vertex2ID]->coords[2];

		int x3_C = mesh->verts[vertex3ID]->coords[0];
		int y3_C = mesh->verts[vertex3ID]->coords[1];
		int z3_C = mesh->verts[vertex3ID]->coords[2];

		int x_B_A = x2_B - x1_A;
		int y_B_A = y2_B - y1_A;
		int z_B_A = z2_B - z1_A;

		int x_C_A = x3_C - x1_A;
		int y_C_A = y3_C - y1_A;
		int z_C_A = z3_C - z1_A;

		int v_B_A[3] = { x_B_A, y_B_A, z_B_A };
		int v_C_A[3] = { x_C_A, y_C_A, z_C_A };
		int c_P[3];
		crossProductFunction(v_B_A, v_C_A, c_P);
		printf("(B-A): x=%d y=%d z=%d\n", x_B_A, y_B_A, z_B_A);
		printf("(C-A): x=%d y=%d z=%d\n", x_C_A, y_C_A, z_C_A);
		printf("CrossProduct x=%d y=%d z=%d\n", c_P[0], c_P[1], c_P[2]);

		float len = calculateLength(c_P);
		float len_B_A = calculateLength(v_B_A);
		float len_C_A = calculateLength(v_C_A);
		float sinAlpha = len / (len_B_A * len_C_A);
		float arcSin = asin(sinAlpha);
		float norm[3];
		normalize(c_P, norm);
		float norm_final[3];
		norm_final[0] = norm[0] * arcSin;
		norm_final[1] = norm[1] * arcSin;
		norm_final[2] = norm[2] * arcSin;
		printf("Normal x=%f y=%f z=%f\n\n", norm_final[0], norm_final[1], norm_final[2]);

		//int x1_A = mesh->verts[vertex1ID]->coords[0];
		//int y1_A = mesh->verts[vertex1ID]->coords[1];
		//int z1_A = mesh->verts[vertex1ID]->coords[2];

		float p[3];
		p[0] = x1_A;
		p[1] = y1_A;
		p[2] = z1_A;

		float d[3];
		/*d[0] = norm_final[0];
		d[1] = norm_final[1];
		d[2] = norm_final[2];*/
		d[0] = c_P[0];
		d[1] = c_P[1];
		d[2] = c_P[2];
		
		for (size_t j = 0; j < mesh->tris.size(); j++)
		{
			if (i != j) {
				int idxofVertex0inTriangle = mesh->tris[j]->v1i;
				int idxofVertex1inTriangle = mesh->tris[j]->v2i;
				int idxofVertex2inTriangle = mesh->tris[j]->v3i;
				float v0[3];
				float v1[3];
				float v2[3];

				float x0 = mesh->verts[idxofVertex0inTriangle]->coords[0];
				float x1 = mesh->verts[idxofVertex1inTriangle]->coords[0];
				float x2 = mesh->verts[idxofVertex2inTriangle]->coords[0];
				v0[0] = x0;
				v1[0] = x1;
				v2[0] = x2;


				float y0 = mesh->verts[idxofVertex0inTriangle]->coords[1];
				float y1 = mesh->verts[idxofVertex1inTriangle]->coords[1];
				float y2 = mesh->verts[idxofVertex2inTriangle]->coords[1];
				v0[1] = y0;
				v1[1] = y1;
				v2[1] = y2;
				
				float z0 = mesh->verts[idxofVertex0inTriangle]->coords[2];
				float z1 = mesh->verts[idxofVertex1inTriangle]->coords[2];
				float z2 = mesh->verts[idxofVertex2inTriangle]->coords[2];
				v0[2] = z0;
				v1[2] = z1;
				v2[2] = z2;
				float intersects = rayIntersectsTriangle(p, d, v0, v1, v2);
				if (intersects > 0) {
					printf("p =%f.2, %f.2, %f.2\n", p[0], p[1], p[2]);
					printf("d =%f.2, %f.2, %f.2\n", d[0], d[1], d[2]);
					printf("v0 =%f.2, %f.2, %f.2\n", v0[0], v0[1], v0[2]);
					printf("v1 =%f.2, %f.2, %f.2\n", v1[0], v1[1], v1[2]);
					printf("v2 =%f.2, %f.2, %f.2\n", v2[0], v2[1], v2[2]);
					printf("Intersects = %f.2 \n\n\n\n", intersects);
					if (mesh->verts[vertex1ID]->length < intersects) {
						mesh->verts[vertex1ID]->length = intersects;
						mesh->verts[vertex2ID]->length = intersects;
						mesh->verts[vertex3ID]->length = intersects;
					}
				}
			}
		}
	}
}

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

	mesh->loadOff("0.off");
//	mesh->createCube(20.0f);

	triangleNormalVectors(mesh);

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
