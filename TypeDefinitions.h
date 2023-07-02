#pragma once

#include <array>

typedef std::array<int, 3> triVertsIds;
typedef std::array<float, 3> vertCoords;
typedef std::array<std::array<float, 3>, 3> triVertsCoords;
typedef std::array<std::array<float, 3>, 2> triOtherVertsCoords;

namespace TD {
	inline void fillWith(float* retArr, vertCoords arr);

	void fillWith(float* retArr, const vertCoords arr) {
		int size = arr.size();
		for (size_t i = 0; i < size; i++)
		{
			retArr[i] = arr[i];
		}
	}
}