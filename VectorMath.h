#pragma once

#include <vector>
#include <array>
#include <cmath>

namespace VectorMath
{
	inline void vector(float v_A[3], const float v_B[3], const float v_C[3]);
	inline void crossProduct(float CP[3], const float v_A[3], const float v_B[3]);
	inline float calculateLengthOfVector(const float v[3]);
	inline float distanceBetweenVectors(const float v[3], const float y[3]);
	inline float innerProduct(const float v[3], const float q[3]);
	inline float rayTriangleIntersectLength(const float* p, const float* d, const float* v0, const float* v1, const float* v2);
	inline float calculateAngleBetweenVectors(const float* vertex1, const float* vertex2);
	inline void midpoint(float v_A[3], const float v_B[3], const float v_C[3]);


	//v_A = v_B - v_C
	void vector(float v_A[3], const float v_B[3], const float v_C[3])
	{
		v_A[0] = v_B[0] - v_C[0];
		v_A[1] = v_B[1] - v_C[1];
		v_A[2] = v_B[2] - v_C[2];
	}
	
	//CP = v_A x v_B
	void crossProduct(float CP[3], const float v_A[3], const float v_B[3])
	{
		CP[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
		CP[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
		CP[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
	}

	//v . q
	float innerProduct(const float v[3], const float q[3])
	{
		return v[0] * q[0] + v[1] * q[1] + v[2] * q[2];
	}

	//|v|
	float calculateLengthOfVector(const float v[3])
	{
		float square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		float root = sqrt(square);
		return root;
	}
	
	//v --- y
	float distanceBetweenVectors(const float v[3], const float y[3])
	{
		const float deltaX = y[0] - v[0];
		const float deltaY = y[1] - v[1];
		const float deltaZ = y[2] - v[2];

		//return fabs(sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ));
		return sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
	}

	float rayTriangleIntersectLength(const float* p, const float* d, const float* v0, const float* v1, const float* v2)
	{
		//https://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
		float e1[3], e2[3], h[3], s[3], q[3];
		float a, f, u, v;
		vector(e1, v1, v0);
		vector(e2, v2, v0);
		crossProduct(h, d, e2);
		a = innerProduct(e1, h);
		if (a > -0.00001 && a < 0.00001)
			return -1;

		f = 1 / a;
		vector(s, p, v0);
		u = f * (innerProduct(s, h));

		if (u < 0.0 || u > 1.0)
			return -1;

		crossProduct(q, s, e1);
		v = f * innerProduct(d, q);

		if (v < 0.0 || u + v > 1.0)
			return -1;

		// at this stage we can compute t to find out where
		// the intersection point is on the line
		float t = f * innerProduct(e2, q);

		if (t > 0.00001) // ray intersection
		{
			const float directionVector[3] = {d[0], d[1], d[2]};
			const float lenghtOfDirectionVector = calculateLengthOfVector(directionVector);
			return t * lenghtOfDirectionVector;
		}
		else // this means that there is a line intersection
			 // but not a ray intersection
			return -2;
	}

	float calculateAngleBetweenVectors(const float* vertex1, const float* vertex2)
	{
		// As described at https://www.jwwalker.com/pages/angle-between-vectors.html
		const float dot12 = innerProduct(vertex1, vertex2);
		const float dot11 = innerProduct(vertex1, vertex1);
		const float dot22 = innerProduct(vertex2, vertex2);
		return std::acos(dot12 / std::sqrt(dot11 * dot22));
	}

	//v_A : midpoint between v_B and v_C
	void midpoint(float v_A[3], const float v_B[3], const float v_C[3])
	{
		// Calculate the midpoint between v_B and v_C
		v_A[0] = (v_B[0] + v_C[0]) / 2.0f;
		v_A[1] = (v_B[1] + v_C[1]) / 2.0f;
		v_A[2] = (v_B[2] + v_C[2]) / 2.0f;
	}
}
