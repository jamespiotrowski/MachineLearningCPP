#pragma once

#include "Array.h"
//#include <iostream>

using namespace std;

template<typename T>
void AllCombinations(const int &length, const Array<T> &vals, const Array<T> &currentCombination, Array<Array<T>> &combinations) {
	// Base Case
	if (currentCombination.GetSize() == length) {
		combinations.Add(currentCombination);
		return;
	}
	for (int i = 0; i < vals.GetSize(); i++) {
		Array<T> newArr = currentCombination;
		newArr.Add(vals.At(i));
		AllCombinations(length, vals, newArr, combinations);
	}
	return;
}

template<typename T>
Array<Array<T>> AllCombinations(const Array<T>& vals, const unsigned int& len) {
	Array<Array<T>> c;
	if (vals.GetSize() == 0) {
		return c;
	}
	Array<T> currentCombination;
	AllCombinations(len, vals, currentCombination, c);
	return c;
}

template<typename T>
void SomeCombinations(const int& length, const Array<T>& vals, const Array<T>& currentCombination, Array<Array<T>>& combinations, const int& numToKeep, int& numCombinations) {
	// Base Case
	if (currentCombination.GetSize() == length) {
		if ((numCombinations % 100) < numToKeep) {
			combinations.Add(currentCombination);
		}
		return;
	}
	for (int i = 0; i < vals.GetSize(); i++) {
		numCombinations += 1;
		Array<T> newArr = currentCombination;
		newArr.Add(vals.At(i));
		SomeCombinations(length, vals, newArr, combinations, numToKeep, numCombinations);
	}
	return;
}

template<typename T>
Array<Array<T>> SomeCombinations(const Array<T>& vals, const unsigned int& len, const double& percentage) {
	Array<Array<T>> c;
	if (vals.GetSize() == 0 || percentage > 1.0 || percentage <= 0.0) {
		return c;
	}

	int numToKeep = 100 * percentage;
	int numCombinations = 0;

	/*cout << "Value Set Contains: " << vals.GetSize() << endl
		<< "Generating Combinations of Length: " << len << endl
		<< "Only Keeping " << numToKeep << "% of Combinations generated" << endl;*/

	Array<T> currentCombination;
	SomeCombinations(len, vals, currentCombination, c, numToKeep, numCombinations);

	return c;
}

