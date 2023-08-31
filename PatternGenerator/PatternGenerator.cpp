#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "BitMap.h"
using namespace std;

#define M_PI 3.141592653589793;

#define SCALES_SQUARE 1
#define SCALES_DIAMOND 2
#define SCALES_CIRCLE 1
#define SCALES_TRIANGLE 2
#define SCALES_STAR 1
#define SCALES_PENTAGON 1
#define SCALES_STRIPE 1

/**********************************************************************************************
###############################################################################################
##### ARRAY CLASS
###############################################################################################
**********************************************************************************************/
template <class T> class Array {
private:

	/* Array Members */
	unsigned int arraySize = 0;
	unsigned int maxSize = 0;
	T* arr = nullptr;
	unsigned int growthFactor = 100;

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void resize() {
		maxSize += growthFactor;
		T* newArr = new T[maxSize];
		for (unsigned int i = 0; i < arraySize; i++) {
			newArr[i] = arr[i];
		}
		delete[] arr;
		arr = newArr;
		newArr = nullptr;
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void clear() {
		if (arr != nullptr) {
			delete[] arr;
		}
		maxSize = 0;
		arraySize = 0;
		arr = nullptr;
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	int partition(int start, int end)
	{
		T pivot = arr[start];
		int count = 0;
		for (int i = start + 1; i <= end; i++) {
			if (arr[i] <= pivot)
				count++;
		}

		// Giving pivot element its correct position
		int pivotIndex = start + count;
		swap(arr[pivotIndex], arr[start]);

		// Sorting left and right parts of the pivot element
		int i = start, j = end;

		while (i < pivotIndex && j > pivotIndex) {

			while (arr[i] <= pivot) {
				i++;
			}

			while (arr[j] > pivot) {
				j--;
			}

			if (i < pivotIndex && j > pivotIndex) {
				swap(arr[i++], arr[j--]);
			}
		}

		return pivotIndex;
	}

	void quickSort(int start, int end)
	{

		// base case
		if (start >= end)
			return;

		// partitioning the array
		int p = partition(start, end);

		// Sorting the left part
		quickSort(start, p - 1);

		// Sorting the right part
		quickSort(p + 1, end);
	}

public:

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	Array() : arraySize(0), maxSize(0), arr(nullptr), growthFactor(1) { }

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	~Array() { clear(); }

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void operator=(const Array& copy) {
		clear();
		arraySize = copy.arraySize;
		maxSize = copy.maxSize;
		arr = new T[maxSize];
		growthFactor = 100;
		for (unsigned int i = 0; i < arraySize; i++) {
			arr[i] = copy.arr[i];
		}
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	Array(const Array& copy) {
		(*this) = copy;
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void Add(T item) {
		if (arraySize == maxSize) {
			resize();
		}
		arr[arraySize] = item;
		arraySize += 1;
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	T& operator[](unsigned int i) {
		return arr[i];
	}

	T& At(unsigned int i) const {
		return arr[i];
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	unsigned int GetSize() const { return arraySize; }

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void reset() { arraySize = 0; }

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	bool exists(T item) { 
		for (unsigned int i = 0; i < arraySize; i++) {
			if (item == arr[i]) {
				return true;
			}
		}
		return false;
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void sort() {
		quickSort(0, arraySize - 1);
	}

	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void remove(unsigned int i) {
		if (i < 0 || i >= arraySize) {
			return;
		}
		for (unsigned int j = i; j < arraySize - 1; j++) {
			arr[j] = arr[j + 1];
		}
		arraySize = arraySize - 1;
	}


	/****************************************************************
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-

	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	****************************************************************/
	void removeDuplicates() { // WARNING THIS WILL SORT YOUR ARRAY
		sort();
		unsigned int i = 0;
		while (i < (arraySize - 1)) {
			if (arr[i] == arr[i + 1]) {
				remove(i + 1);
			}
			else {
				i += 1;
			}
		}
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
enum PatternType {
	SQUARE
	, CIRCLE 
};

int RoundDouble(double d) {
	int i = d;
	double decimal = d - i;
	return (decimal >= 0.5) ? i + 1 : i;
}

struct Coordinate {
	int y;
	int x;
	Coordinate() { }
	Coordinate(int y, int x) : x(x), y(y) { }
	string ToString() const { return ("{" + to_string(y) + "," + to_string(x) + "}"); }
	bool operator==(const Coordinate& c) const { return ((y == c.y) && (x == c.x)); }
	bool operator<=(const Coordinate& c) const {
		if (this->y < c.y) { return true; }
		if (this->y > c.y) { return false; }
		return (this->x <= c.x);
	}
	bool operator>(const Coordinate& c) const { return !((*this) <= c); }
};

class Angle {
public:
	Coordinate c;
	bool intersect = true;

	Angle() : c(0, 0), intersect(false){ }

	Angle(Coordinate c1, bool i) : c(c1){
		intersect = i;
	}

	void operator=(const Angle& copy) {
		c = copy.c;
		intersect = copy.intersect;
	}

	Angle(const Angle& copy) {
		(*this) = copy;
	}

	bool operator==(const Angle& a) { return this->c == a.c; }
	bool operator<=(const Angle& a) { return (this->c <= a.c); }
	bool operator>(const Angle& a) { return !((*this) <= a); }

	string ToString() const {
		string s = "[" + c.ToString() + ",";
		if (intersect) {
			s += "i]";
		}
		else {
			s += "t]";
		}
		return s;
	}
};

class Edge {
public:
	Coordinate c1, c2;
	double slope = 1;
	double b = 0;
	bool infSlope = false;

	Edge() {}

	Edge(Coordinate ca, Coordinate cb) {
		if (cb > ca) {
			c1 = ca;
			c2 = cb;
		}
		c1 = cb;
		c2 = ca;

		if ((c2.x - c1.x) == 0) {
			infSlope = true;
		}
		else {
			slope = (double)(c2.y - c1.y) / (double)(c2.x - c1.x);
			b = (double)c1.y - (slope * (double)c1.x);
		}
	}

	bool operator==(const Edge& e) {
		return (c1 == e.c1) && (c2 == e.c2);
	}

	bool operator<=(const Edge& e) {
		return (c1 <= e.c1) && (c2 <= e.c2);
	}

	bool operator>(const Edge& e) {
		return !((*this) <= e);
	}

	static bool SharesPoint(const Edge& edge1, const Edge& edge2) {
		return (edge1.c1 == edge2.c1) || (edge1.c1 == edge2.c2) || (edge1.c2 == edge2.c1) || (edge1.c2 == edge2.c2);
	}

	static Coordinate GetSharedPoint(const Edge& edge1, const Edge& edge2){
		if (edge1.c1 == edge2.c1) { return edge1.c1; }
		if (edge1.c1 == edge2.c2) { return edge1.c1; }
		if (edge1.c2 == edge2.c1) { return edge1.c2; }
		if (edge1.c2 == edge2.c2) { return edge1.c2; }
		cout << "ERROR: Edge::GetSharedPoint - Edges do not share a point, returning (0,0)" << endl;
		return Coordinate(0, 0);
	}

	string ToString() const {
		string s = "[" + c1.ToString() + " - " + c2.ToString() + "] :";
		if (!infSlope) {
			s += to_string(slope) + " ";
		}
		else {
			s += "(inf) ";
		}
		return s;
	}

	double getValueAtY(unsigned int y) const {
		return ((infSlope) ? c2.x : ((((double)y) - b) / slope));
	}
};

class Polygon {
private:
	Array<Coordinate> points;
	Array<Edge> edges;
	Array<Angle> angles;

	unsigned int maxY, minY, maxX, minX;

public:

	bool pointWithinPolygonRange(const Coordinate& c) const {
		return c.y >= minY && c.y <= maxY && c.x <= maxX && c.x >= minX;
	}

	static Array<Coordinate> ComputeStraitLine(Coordinate c1, Coordinate c2) {
		Array<Coordinate> line;
		if (c1.y == c2.y) {
			int greaterX = (c2.x > c1.x) ? c2.x : c1.x;
			int lesserX = (c2.x > c1.x) ? c1.x : c2.x;
			for (unsigned int i = lesserX + 1; i < greaterX; i++) {
				line.Add(Coordinate(c2.y, i));
			}
		}
		else {
			int numCoordinates = (c2.y > c1.y) ? c2.y - c1.y : c1.y - c2.y;
			double slope = 0;
			double b = 0;
			bool infSlope = false;

			if ((c2.x - c1.x) == 0) {
				infSlope = true;
			}
			else {
				slope = (double)(c2.y - c1.y) / (double)(c2.x - c1.x);
				b = c1.y - (slope * c1.x);
			}

			int startY = (c2.y > c1.y) ? c1.y : c2.y;
			Coordinate c;
			double d;
			for (int i = 0; i < numCoordinates; i += 1) {
				c.y = i + startY;
				c.x = (infSlope) ? c2.x : RoundDouble((((double)(c.y)) - b) / slope);
				line.Add(c);
			}
		}
		if (!line.exists(c1)) { line.Add(c1); }
		if (!line.exists(c2)) { line.Add(c2); }
		return line;
	}

	// Angle from the x = 0 line
	static Coordinate ComputePointGivenAngleAndDistance(double angle, double distance, Coordinate startingPoint) {
		double piCalc = angle * M_PI;
		piCalc = piCalc / 180.0;
		double x = startingPoint.x + (distance * sin(piCalc));
		double y = startingPoint.y + (distance * cos(piCalc));
		return Coordinate(RoundDouble(y), RoundDouble(x));
	}

	static Coordinate ComputeCentroid(Coordinate c1, Coordinate c2, Coordinate c3) {
		double new_x = ((double)(c1.x + c2.x + c3.x)) / 3.0;
		double new_y = ((double)(c1.y + c2.y + c3.y)) / 3.0;
		return Coordinate(RoundDouble(new_y), RoundDouble(new_x));
	}

	static bool inLineWithAngle(const Coordinate& c, const Angle& a) {
		return (c.y == a.c.y) && (c.x <= a.c.x);
	}

	static bool intersectsAngle(const Coordinate& c, const Angle& a) {
		if (c.y == a.c.y) {
			if (c.x <= a.c.x) {
				return a.intersect;
			}
		}
		return false;
	}

	static bool intersectsEdge(const Coordinate& c, const Edge& e) {
		if ((c.y <= e.c1.y && c.y >= e.c2.y) || (c.y <= e.c2.y && c.y >= e.c1.y)) {
			double x = e.getValueAtY(c.y);
			return x >= ((double)c.x);
		}
		return false;
	}

	bool isInsidePolygon(const Coordinate& c) const {

		if (pointWithinPolygonRange(c)) {
			unsigned int intersections = 0;
			for (unsigned int i = 0; i < edges.GetSize(); i++) {
				if (intersectsEdge(c, edges.At(i))) { // If ray would intersect an edge
					intersections += 1;
				}
			}

			for (unsigned int i = 0; i < angles.GetSize(); i++) {
				if (inLineWithAngle(c, angles.At(i))) { // If ray would intersect an edge
					if (intersectsAngle(c, angles.At(i))) {
						intersections -= 1; // would have double counted for two edges meeting, sub 1
					}
					else {
						intersections -= 2; // would have double counted for two edges meeting, sub 2 (because it doesnt intersect the angle)
					}
				}
			}

			return ((intersections % 2) != 0);
		}
		return false;
	}

	Polygon(Array<Edge> e) : edges(e) {
		// Remove dupes from list of edges
		edges.removeDuplicates();

		// Get all points
		for (unsigned int i = 0; i < edges.GetSize(); i++) {
			Array<Coordinate> line = ComputeStraitLine(edges[i].c1, edges[i].c2);
			for (unsigned int j = 0; j < line.GetSize(); j++) {
				points.Add(line[j]);
			}
		}
		// Remove dupes
		points.removeDuplicates();

		// Find angles
		for (unsigned int i = 0; i < edges.GetSize(); i++) {
			for (unsigned int j = 0; j < edges.GetSize(); j++) {
				if (i == j) {
					continue;
				}
				// If edges share a point
				if (Edge::SharesPoint(edges[i], edges[j])) {
					// Get the point
					Coordinate c = Edge::GetSharedPoint(edges[i], edges[j]);
					// Get the other 2 coordinates
					Coordinate other1 = (edges[i].c1 == c) ? edges[i].c2 : edges[i].c1;
					Coordinate other2 = (edges[j].c1 == c) ? edges[j].c2 : edges[j].c1;
					// Determine if a ray would intersect the two edges or just touch it
					bool intersection = ((other1.y > c.y) && (other2.y < c.y)) || ((other1.y < c.y) && (other2.y > c.y));
					Angle a(c, intersection);
					angles.Add(a);
				}
			}
		}
		// Remove dupes - not doing this because bug with copying the intersection bool between objects
		angles.removeDuplicates();	

		if (points.GetSize() > 0) {
			maxY = points[0].y;
			minY = points[0].y;
			maxX = points[0].x;
			minX = points[0].x;
			for (unsigned int i = 1; i < points.GetSize(); i++) {
				if (points[i].y > maxY) { maxY = points[i].y; }
				if (points[i].y < minY) { minY = points[i].y; }
				if (points[i].x > maxX) { maxX = points[i].x; }
				if (points[i].x < minX) { minX = points[i].x; }
			}
		}
		else {
			maxY = 0;
			minY = 0;
			maxX = 0;
			minX = 0;
		}
	}

	void plotPolygon(char** grid) {
		for (unsigned int i = 0; i < points.GetSize(); i++) {
			grid[points[i].y][points[i].x] = 1;
		}
	}

	void printPolygon() {
		string p = "Points : { ";
		for (unsigned int i = 0; i < points.GetSize(); i++) {
			p += points[i].ToString() + " ";
		}
		p += " }\n";

		string e = "Edges : { ";
		for (unsigned int i = 0; i < edges.GetSize(); i++) {
			e += edges[i].ToString() + " ";
		}
		e += " }\n";

		string a = "Angles : { ";
		for (unsigned int i = 0; i < angles.GetSize(); i++) {
			a += angles[i].ToString() + " ";
		}
		a += " }\n";

		string r = "Width Range: {" + to_string(minX) + " - " + to_string(maxX)
			+ "} | Height Range: {" + to_string(minY) + " - " + to_string(maxY) + "}\n";

		string t = p + e + a + r;
		cout << t << endl;
	}
	
	
	// way to determine if point intersects with an angle or if it just touches it
};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class UnitPattern {
protected:
	char** pattern = nullptr;
	int height = 0, width = 0, numberOfScales = 0;
	double* scales = nullptr;

private: 

	// possible optimization..
	void FillInUntilEdge(const int& y, int& w) {
		while (w < width && pattern[y][w] != 1) {
			pattern[y][w] = 1;
			w += 1;
		}
	}

protected:

	void GetCenter(const int& height, const int& width, int& centerHeight, int& centerWidth) {
		centerHeight = height / 2;
		centerWidth = width / 2;
	}

	void SetScale(int numberOfScales, double* scales) {
		this->scales = new double[numberOfScales];
		for (int i = 0; i < numberOfScales; i++) {
			this->scales[i] = scales[i];
		}
	}
		
	// For round.. use for circle.. maybe
	bool IsInsidePolygon(Coordinate c, const Array<Coordinate>& loneEdgePoints) {
		int edges = 0;
		int sp = width - 1;
		for (int w = c.x; w < width; w += 1) {
			if (w < sp) {
				if (pattern[c.y][w] == 1 && pattern[c.y][w + 1] != 1) {
					edges += 1;
				}
			}
			else {
				if (pattern[c.y][w] == 1 && pattern[c.y][w - 1] != 1) {
					edges += 1;
				}
			}
		}

		for (unsigned int i = 0; i < loneEdgePoints.GetSize(); i++) {
			if (loneEdgePoints.At(i).y == c.y && loneEdgePoints.At(i).x >= c.x) {
				edges -= 1;
			}
		}

		return (((edges % 2) == 1));
	}

	// For a true polygon
	void FillInPolygon(const Polygon& p) {
		Coordinate c;
		for (unsigned int h = 0; h < height; h++) {
			c.y = h;
			for (unsigned int w = 0; w < width; w++) {
				c.x = w;
				if (p.isInsidePolygon(c)) {
					pattern[h][w] = 1; // Could optimize by using fillInUntilEdge, but not general and may not be necesarry
				}
			}
		}
	}

	// for a not true polygon, though... a circle could just be a bunch of edges. we will see
	void FillInPolygon(const Array<Coordinate> &loneEdgePoints) {
		Coordinate c;
		for (unsigned int h = 0; h < height; h++) {
			c.y = h;
			for (unsigned int w = 0; w < width; w++) {
				c.x = w;
				if (IsInsidePolygon(c, loneEdgePoints)) {
					pattern[h][w] = 1; // Could optimize by using fillInUntilEdge, but not general and may not be necesarry
				}
			}
		}
	}

	/* this.... shouldnt be used. */
	void FillInBruteForce(const Coordinate& c) {
		if (c.x < 0 || c.x >= width || c.y < 0 || c.y >= height) { return; }
		if (pattern[c.y][c.x] != 0) { return; }
		pattern[c.y][c.x] = 1;
		/* Down */
		if ((c.y + 1) < height) {
			if (pattern[c.y + 1][c.x] != 1) {
				Coordinate down(c.y + 1, c.x);
				FillInBruteForce(down);
			}
		}
		/* Up */
		if ((c.y - 1) >= 0) {
			if (pattern[c.y - 1][c.x] != 1) {
				Coordinate up(c.y - 1, c.x);;
				FillInBruteForce(up);
			}
		}
		/* right */
		if ((c.x + 1) < width) {
			if (pattern[c.y][c.x + 1] != 1) {
				Coordinate right(c.y, c.x + 1);
				FillInBruteForce(right);
			}
		}
		/* left */
		if ((c.x - 1) >= 0) {
			if (pattern[c.y][c.x - 1] != 1) {
				Coordinate left(c.y, c.x - 1);
				FillInBruteForce(left);
			}
		}
	}


private:

	static void DetermineHeightAndWidthWithTrueCenter(int& height, int& width) {
		if ((height % 2) == 0) {
			height += 1;
		}

		if ((width % 2) == 0) {
			width += 1;
		}
	}

	void clear() {
		for (int i = 0; i < height; i++) {
			delete[] pattern[i];
		}
		delete[] pattern;
		delete[] scales;
		pattern = nullptr;
		height = 0;
		width = 0;
		numberOfScales = 0;
	}

	char ToString(char p) { return (p == 0) ? '0' : '1'; }
	string ToSymbol(char p) { return (p == 0) ? "." : "#"; }

public:

	UnitPattern(int height, int width) : height(height), width(width){
		DetermineHeightAndWidthWithTrueCenter(this->height, this->width);
		pattern = new char* [this->height];
		for (int i = 0; i < this->height; i++) {
			pattern[i] = new char[this->width];
			for (int j = 0; j < this->width; j++) {
				pattern[i][j] = 0;
			}
		}
	}

	~UnitPattern() {
		clear();
	}

	void operator=(const UnitPattern& copy) {
		clear();
		numberOfScales = copy.numberOfScales;
		height = copy.height;
		width = copy.width;
		pattern = new char* [height];
		for (int i = 0; i < height; i++) {
			pattern[i] = new char[width];
			for (int j = 0; j < width; j++) {
				pattern[i][j] = copy.pattern[i][j];
			}
		}

		this->scales = new double[numberOfScales];
		for (int i = 0; i < numberOfScales; i++) {
			scales[i] = copy.scales[i];
		}
	}

	UnitPattern(const UnitPattern& copy) {
		(*this) = copy;
	}

	void PrintPattern(ostream& out) {
		string s = "";
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				s += ToSymbol(pattern[i][j]);
			}
			s += "\n";
		}
		s += "\n";

		out << s;
	}

	int GetHeight() const { return height; }
	int GetWidth() const { return width; }

	char At(const int& h, const int& w) const {
		return pattern[h][w];
	}

	virtual void GenerateUnitPattern() = 0;

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class SquarePattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth);					 // Detmine center coordinates
		int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit
		/* Compute the pattern area */
		int startHeight = centerHeight - radius;
		int startWidth = centerWidth - radius;
		int endHeight = centerHeight + radius;
		int endWidth = centerWidth + radius;

		/* Compute the pattern area */
		Coordinate* c = new Coordinate[4];
		c[0] = Coordinate(startHeight, startWidth);
		c[1] = Coordinate(startHeight, endWidth);
		c[2] = Coordinate(endHeight, endWidth);
		c[3] = Coordinate(endHeight, startWidth);

		Array<Edge> edges;
		edges.Add(Edge(c[0], c[1]));
		edges.Add(Edge(c[1], c[2]));
		edges.Add(Edge(c[2], c[3]));
		edges.Add(Edge(c[3], c[0]));

		delete[] c;

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	SquarePattern(int height, int width, double scale) : UnitPattern(height, width) {
		SetScale(SCALES_SQUARE, &scale);
		GenerateUnitPattern();
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class DiamondPattern : public UnitPattern {
private:
	const int widthScale = 0;
	const int heightScale = 1;

	void GenerateUnitPattern() {
		int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth); // Detmine center coordinates
		int widthRadius = (scales[widthScale] * width) / 2;
		int heightRadius = (scales[heightScale] * height) / 2;

		/* Compute the pattern area */
		Coordinate* c = new Coordinate[4];
		c[0] = Coordinate(centerHeight + heightRadius, centerWidth);
		c[1] = Coordinate(centerHeight, centerWidth - widthRadius);
		c[2] = Coordinate(centerHeight - heightRadius, centerWidth);
		c[3] = Coordinate(centerHeight, centerWidth + widthRadius);

		Array<Edge> edges;
		edges.Add(Edge(c[0], c[1]));
		edges.Add(Edge(c[1], c[2]));
		edges.Add(Edge(c[2], c[3]));
		edges.Add(Edge(c[3], c[0]));

		delete[] c;

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	DiamondPattern(int height, int width, double widthScale, double heightScale) : UnitPattern(height, width) {
		double* scales = new double[SCALES_DIAMOND];
		scales[0] = widthScale;
		scales[1] = heightScale;
		SetScale(SCALES_DIAMOND, scales);
		GenerateUnitPattern();
		delete[] scales;
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class CirclePattern : public UnitPattern {
private:
	const int scale1 = 0;

	int PythagoreanTheorem_Edge(int a_or_b, double c) {
		return sqrt((c * c) - (a_or_b * a_or_b));
	}

	void GenerateUnitPattern() {
		int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth); // Detmine center coordinates
		double radius = (scales[scale1] * ((height < width) ? height : width)) / 2;
		/* Compute the pattern area */
		Array<Coordinate> circlePoints;
		Coordinate coords[8];

		for (int i = 1; i <= (int)radius; i++) {
			int yPoint = PythagoreanTheorem_Edge(i, radius);

			coords[0] = Coordinate(centerHeight - yPoint, centerWidth + i);
			coords[1] = Coordinate(centerHeight - yPoint, centerWidth - i);
			coords[2] = Coordinate(centerHeight + yPoint, centerWidth + i);
			coords[3] = Coordinate(centerHeight + yPoint, centerWidth - i);
			coords[4] = Coordinate(centerHeight - i, centerWidth + yPoint);
			coords[5] = Coordinate(centerHeight - i, centerWidth - yPoint);
			coords[6] = Coordinate(centerHeight + i, centerWidth + yPoint);
			coords[7] = Coordinate(centerHeight + i, centerWidth - yPoint);

			for (int j = 0; j < 8; j++) {
				if (!circlePoints.exists(coords[j])) {
					circlePoints.Add(coords[j]);
				}
			}
		}

		coords[0] = Coordinate(centerHeight, centerWidth + radius);
		coords[1] = Coordinate(centerHeight, centerWidth - radius);
		coords[2] = Coordinate(centerHeight + radius, centerWidth);
		coords[3] = Coordinate(centerHeight - radius, centerWidth);
		for (int j = 0; j < 4; j++) {
			if (!circlePoints.exists(coords[j])) {
				circlePoints.Add(coords[j]);
			}
		}

		Array<int> edgePointCounter;
		for (unsigned int i = 0; i < height; i++) {
			edgePointCounter.Add(0);
		}

		for (int i = 0; i < circlePoints.GetSize(); i++) {
			edgePointCounter[circlePoints[i].y] += 1;
			pattern[circlePoints[i].y][circlePoints[i].x] = 1;
		}

		int minHeight = 0, maxHeight = 0;
		for (unsigned int i = 0; i < edgePointCounter.GetSize(); i++) {
			if (edgePointCounter[i] >= 1) {
				minHeight = i;
				break;
			}
		}

		

		for (int i = ((int)edgePointCounter.GetSize()) - 1; i >= 0; i--) {
			if (edgePointCounter[i] >= 1) {
				maxHeight = i;
				break;
			}
		}

		bool foundMinCoord = false, foundMaxCoord = false;
		Coordinate minCoordinate;
		Coordinate maxCoordinate;
		for (unsigned int j = 0; j < circlePoints.GetSize(); j++) {

			if (minHeight == circlePoints[j].y) {
				if (foundMinCoord) {
					if (minCoordinate.x > circlePoints[j].x) {
						minCoordinate = circlePoints[j];
					}
				}
				else {
					minCoordinate = circlePoints[j];
					foundMinCoord = true;
				}
			}

			if (maxHeight == circlePoints[j].y) {
				if (foundMaxCoord) {
					if (maxCoordinate.x > circlePoints[j].x) {
						maxCoordinate = circlePoints[j];
					}
				}
				else {
					maxCoordinate = circlePoints[j];
					foundMaxCoord = true;
				}
			}
		
		}

		Array<Coordinate> loneEdgePoints;
		loneEdgePoints.Add(minCoordinate);
		loneEdgePoints.Add(maxCoordinate);

		/*
		for (unsigned int i = 0; i < edgePointCounter.GetSize(); i++) {
			if (edgePointCounter[i] == 1) {
				for (unsigned int j = 0; j < circlePoints.GetSize(); j++) {
					if (i == circlePoints[j].y) {
						loneEdgePoints.Add(circlePoints[j]);
						break;
					}
				}
			}
		}
		*/

		

		FillInPolygon(loneEdgePoints);
	}

public:

	CirclePattern(int height, int width, double scale) : UnitPattern(height, width) {
		double s = scale;
		SetScale(SCALES_CIRCLE, &s);
		GenerateUnitPattern();
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class TrianglePattern : public UnitPattern {
private:
	const int widthScale = 0;
	const int heightScale = 1;

	void GenerateUnitPattern() {
		int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth); // Detmine center coordinates
		int widthRadius = (scales[widthScale] * width) / 2;
		int heightRadius = (scales[heightScale] * height) / 2;
		/* Compute the pattern area */
		Coordinate highestPoint;
		Coordinate leftMostPoint;
		Coordinate rightMostPoint;
		/* Compute the pattern area */
		Coordinate* c = new Coordinate[3];
		c[0] = Coordinate(centerHeight - heightRadius, centerWidth);
		c[1] = Coordinate(centerHeight + heightRadius, centerWidth - widthRadius);
		c[2] = Coordinate(centerHeight + heightRadius, centerWidth + widthRadius);

		Array<Edge> edges;
		edges.Add(Edge(c[0], c[1]));
		edges.Add(Edge(c[1], c[2]));
		edges.Add(Edge(c[2], c[0]));;

		delete[] c;

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	TrianglePattern(int height, int width, double widthScale, double heightScale) : UnitPattern(height, width) {
		double* scales = new double[SCALES_TRIANGLE];
		scales[0] = widthScale;
		scales[1] = heightScale;
		SetScale(SCALES_TRIANGLE, scales);
		GenerateUnitPattern();
		delete[] scales;
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class PentagonPattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x);					 // Detmine center coordinates
		int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit

		double split = 360.0 / 5.0;

		/* Compute the pattern area */
		Coordinate* c = new Coordinate[5];
		c[0] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (0.0 * split), radius, centerCoord);
		c[1] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (1.0 * split), radius, centerCoord);
		c[2] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (2.0 * split), radius, centerCoord);
		c[3] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (3.0 * split), radius, centerCoord);
		c[4] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (4.0 * split), radius, centerCoord);

		Array<Edge> edges;
		edges.Add(Edge(c[0], c[1]));
		edges.Add(Edge(c[1], c[2]));
		edges.Add(Edge(c[2], c[3]));
		edges.Add(Edge(c[3], c[4]));
		edges.Add(Edge(c[4], c[0]));

		delete[] c;

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	PentagonPattern(int height, int width, double scale) : UnitPattern(height, width) {
		SetScale(SCALES_PENTAGON, &scale);
		GenerateUnitPattern();
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class StarPattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x);					 // Detmine center coordinates
		int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit

		double split = 360.0 / 5.0;

		/* Compute the pattern area */
		Coordinate* c = new Coordinate[10];
		c[0] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (0.0 * split), radius, centerCoord);
		c[2] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (1.0 * split), radius, centerCoord);
		c[4] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (2.0 * split), radius, centerCoord);
		c[6] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (3.0 * split), radius, centerCoord);
		c[8] = Polygon::ComputePointGivenAngleAndDistance(180.0 + (4.0 * split), radius, centerCoord);

		c[1] = Polygon::ComputeCentroid(c[0], c[2], centerCoord);
		c[3] = Polygon::ComputeCentroid(c[2], c[4], centerCoord);
		c[5] = Polygon::ComputeCentroid(c[4], c[6], centerCoord);
		c[7] = Polygon::ComputeCentroid(c[6], c[8], centerCoord);
		c[9] = Polygon::ComputeCentroid(c[8], c[0], centerCoord);


		Array<Edge> edges;
		edges.Add(Edge(c[0], c[1]));
		edges.Add(Edge(c[1], c[2]));
		edges.Add(Edge(c[2], c[3]));
		edges.Add(Edge(c[3], c[4]));
		edges.Add(Edge(c[4], c[5]));
		edges.Add(Edge(c[5], c[6]));
		edges.Add(Edge(c[6], c[7]));
		edges.Add(Edge(c[7], c[8]));
		edges.Add(Edge(c[8], c[9]));
		edges.Add(Edge(c[9], c[0]));

		delete[] c;

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);

		p.printPolygon();
	}

public:

	StarPattern(int height, int width, double scale) : UnitPattern(height, width) {
		SetScale(SCALES_STAR, &scale);
		GenerateUnitPattern();
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class HorizontalStripe : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x);					 // Detmine center coordinates
		int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit
	
		/* Compute the pattern area */
		Coordinate* c = new Coordinate[4];
		c[0] = Coordinate(centerCoord.y + radius, 0);
		c[1] = Coordinate(centerCoord.y - radius, 0);
		c[2] = Coordinate(centerCoord.y - radius, height - 1);
		c[3] = Coordinate(centerCoord.y + radius, height - 1);


		Array<Edge> edges;
		edges.Add(Edge(c[0], c[1]));
		edges.Add(Edge(c[1], c[2]));
		edges.Add(Edge(c[2], c[3]));
		edges.Add(Edge(c[3], c[0]));

		delete[] c;

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	HorizontalStripe(int height, int width, double scale) : UnitPattern(height, width) {
		SetScale(SCALES_STRIPE, &scale);
		GenerateUnitPattern();
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class Pattern {
private:
	int height = 0;
	int width = 0;
	int horizontalOffset = 0;
	int verticalOffset = 0;
	char** canvas = nullptr;

	void clear() {
		for (int i = 0; i < height; i++) {
			delete[] canvas[i];
		}
		delete[] canvas;
		canvas = nullptr;
		height = 0;
		width = 0;
		horizontalOffset = 0;
		verticalOffset = 0;
	}

	char ToString(char p) { return (p == 0) ? '0' : '1'; }
	string ToSymbol(char p) { return (p == 0) ? "." : "#"; }

public:
	Pattern() { }

	Pattern(int height, int width, int verticalOffset, int horizontalOffset, UnitPattern *unitPattern) : verticalOffset(verticalOffset), horizontalOffset(horizontalOffset), height(height), width(width) {
		canvas = new char*[height];
		for (unsigned h = 0; h < height; h++) {
			canvas[h] = new char[width];
			for (int w = 0; w < width; w++) {
				canvas[h][w] = 0;
			}
		}

		GeneratePattern(unitPattern);
	}

	~Pattern() {
		clear();
	}

	void operator=(const Pattern& copy) {
		clear();
		verticalOffset = copy.verticalOffset;
		horizontalOffset = copy.verticalOffset;
		height = copy.height;
		width = copy.width;
		canvas = new char* [height];
		for (int i = 0; i < height; i++) {
			canvas[i] = new char[width];
			for (int j = 0; j < width; j++) {
				canvas[i][j] = copy.canvas[i][j];
			}
		}
	}

	Pattern(const Pattern& copy) {
		(*this) = copy;
	}

	void PrintPattern(ostream& out) {
		string s = "";
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				s += ToSymbol(canvas[i][j]);
			}
			s += "\n";
		}
		s += "\n";

		out << s;
	}

	void GeneratePattern(UnitPattern* unitPattern) {
		if (unitPattern == nullptr) {
			return;
		}

		int heightStep = unitPattern->GetHeight() + verticalOffset;	// The amount of height steps to take whilst pasting the unit pattern
		int widthStep = unitPattern->GetWidth() + horizontalOffset;	// The amount of width steps to take whilst pasting the unit pattern
		int innerHeightLimit = 0;	// Used inside the loop to limit the inner for loops
		int innerWidthLimit = 0;	// Used inside the loop to limit the inner for loops

		for (int oH = 0; oH < height; oH += heightStep) {
			innerHeightLimit = oH + unitPattern->GetHeight();
			innerHeightLimit = (innerHeightLimit > height) ? height : innerHeightLimit;
			for (int oW = 0; oW < width; oW += widthStep) {
				innerWidthLimit = oW + unitPattern->GetWidth();
				innerWidthLimit = (innerWidthLimit > width) ? width : innerWidthLimit;
				for (int iH = oH; iH < innerHeightLimit; iH += 1) {
					for (int iW = oW; iW < innerWidthLimit; iW += 1) {
						canvas[iH][iW] = unitPattern->At(iH - oH, iW - oW);
					}
				}
			}
		}
	}

	void SavePatternToBmp(string fileName) {
		PixelMatrix pm;
		int colorVal;
		for (int h = 0; h < height; h++) {
			pm.push_back(vector<Pixel>());
			for (int w = 0; w < width; w++) {
				colorVal = (canvas[h][w] == 0) ? 255 : 0;
				Pixel p;
				p.red = colorVal;
				p.blue = colorVal;
				p.green = colorVal;
				pm[h].push_back(p);
			}
		}

		Bitmap bm;
		bm.fromPixelMatrix(pm);

		bm.save(fileName);
	}

};

int main()
{		
	
	string fName = "C:\\Users\\james\\Code\\CPP\\MachineLearningCPP\\MachineLearningCPP\\PatternGenerator\\Output\\Outfile.txt";

	Pattern P;

	ofstream outFile;
	outFile.open(fName);

	UnitPattern* sq = //new HorizontalStripe(50, 50, 0.8);
		//new PentagonPattern(50, 50, 0.7);
		new StarPattern(150, 150, 0.95);
		//new TrianglePattern(100, 100, 0.50, 0.50);
		//new CirclePattern(100, 100, 0.99);
		//new DiamondPattern(100, 100, 0.5, 0.5);  
		//new SquarePattern(100, 100, 0.95);
	sq->PrintPattern(outFile);
	P = Pattern(1000, 2000 , 0, 0, sq);
	delete sq;

	P.SavePatternToBmp("C:\\Users\\james\\Code\\CPP\\MachineLearningCPP\\MachineLearningCPP\\PatternGenerator\\Output\\Outfile.bmp");


	outFile.close();

	return 0;
}
