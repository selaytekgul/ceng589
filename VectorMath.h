#pragma once

#include <vector>
#include <array>
#include <cmath>
namespace VectorMath
{
	inline void vector(float v_A[], const float v_B[], const float v_C[]);
	inline void crossProduct(const float v_A[], const float v_B[], float CP[]);
	inline float calculateLengthOfVector(const float v[]);
	inline float innerProduct(const float v[], const float q[]);
	float rayTriangleIntersectLength(float* p, float* d, float* v0, float* v1, float* v2);
	//void normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr);

	/* a = b - c */
	void vector(float v_A[], const float v_B[], const float v_C[]) {
		v_A[0] = v_B[0] - v_C[0];
		v_A[1] = v_B[1] - v_C[1];
		v_A[2] = v_B[2] - v_C[2];
	}
	void crossProduct(const float v_A[], const float v_B[], float CP[]) {
		CP[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
		CP[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
		CP[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
	}

	float innerProduct(const float v[], const float q[]) {
		return v[0] * q[0] + v[1] * q[1] + v[2] * q[2];
	}

	float calculateLengthOfVector(const float v[]) {
		float square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		float root = sqrt(square);
		return root;
	}

	float rayTriangleIntersectLength(float* p, float* d, float* v0, float* v1, float* v2) {
		float e1[3], e2[3], h[3], s[3], q[3];
		float a, f, u, v;
		vector(e1, v1, v0);
		vector(e2, v2, v0);
		crossProduct(d, e2, h);
		a = innerProduct(e1, h);
		if (a > -0.00001 && a < 0.00001)
			return(-1);

		f = 1 / a;
		vector(s, p, v0);
		u = f * (innerProduct(s, h));

		if (u < 0.0 || u > 1.0)
			return(-1);

		crossProduct(s, e1, q);
		v = f * innerProduct(d, q);

		if (v < 0.0 || u + v > 1.0)
			return(-1);

		// at this stage we can compute t to find out where
		// the intersection point is on the line
		float t = f * innerProduct(e2, q);

		if (t > 0.00001) // ray intersection
		{
			const float directionVector[3] = { d[0], d[1], d[2] };
			const float lenghtOfDirectionVector = calculateLengthOfVector(directionVector);
			return(t * lenghtOfDirectionVector);
		}
		else // this means that there is a line intersection
			 // but not a ray intersection
			return (-2);
	}

	//void normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr) {
	//	float minValue = std::numeric_limits<float>::max();
	//	float maxValue = std::numeric_limits<float>::min();
	//
	//	// Find the minimum and maximum values in the input vector
	//	for (float value : inputArr) {
	//		if (value < minValue)
	//			minValue = value;
	//		if (value > maxValue)
	//			maxValue = value;
	//	}
	//	int i = 0;
	//	// Normalize the input vector values and store them in the output vector
	//	for (float value : inputArr) {
	//		float normalizedValue = (value - minValue) / (maxValue - minValue);
	//		outputArr[i] = normalizedValue;
	//		i++;
	//	}
	//}
}
