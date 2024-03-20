#pragma once

#include <array>

typedef std::array<int, 3> triVertsIds;
typedef std::array<int, 2> edgeVertsIds;
typedef std::array<float, 3> vertCoords;
typedef std::array<std::array<float, 3>, 3> triVertsCoords;
typedef std::array<std::array<float, 3>, 2> triOtherVertsCoords;

namespace TD
{
	template <typename T, typename U>
	inline void fillWith(T retArr, const U& arr, const int size);

	template <typename T, typename U>
	void fillWith(T retArr, const U& arr, const int size) {
		for (size_t i = 0; i < size; i++)
			retArr[i] = arr[i];
	}
}