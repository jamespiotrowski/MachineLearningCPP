#include <iostream>
#include <fstream>
#include <string>
#include "BitMap.h"
using namespace std;

#define SCALES_SQUARE 1
#define SCALES_DIAMOND 2


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
	bool operator==(const Coordinate& c) { return ((y == c.y) && (x == c.x)); }
};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class Line {
private: 
	Coordinate* line = nullptr;
	int numCoordinates = 0;

	void clear() {
		delete[] line;
		line = nullptr;
		numCoordinates = 0;
	}

public:

	Line() { }

	Line(Coordinate c1, Coordinate c2) {
		ComputeStraitLine(c1, c2);
	}
	
	~Line() { clear(); }

	void operator=(const Line& copy) {
		clear();
		this->numCoordinates = copy.numCoordinates;
		this->line = new Coordinate[numCoordinates];
		for (int i = 0; i < numCoordinates; i++) {
			line[i] = copy.line[i];
		}
	}

	Line(const Line& copy) {
		(*this) = copy;
	}

	void ComputeStraitLine(Coordinate c1, Coordinate c2) {
		numCoordinates = (c2.y > c1.y) ? c2.y - c1.y : c1.y - c2.y;
		line = new Coordinate[numCoordinates];
		double slope = (double)(c2.y - c1.y) / (double)(c2.x - c1.x);
		double b = c1.y - (slope * c1.x);
		
		int startY = (c2.y > c1.y) ? c1.y : c2.y;
		Coordinate c;
		double d;
		for (int i = 0; i < numCoordinates; i += 1) {
			c.y = i + startY;
			c.x = RoundDouble((((double)(c.y)) - b) / slope);
			line[i] = c;
		}
	}

	Coordinate operator[](int i) {
		return line[i];
	}

	int GetNumCoordinates() const {
		return numCoordinates;
	}

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

	void PlotLines(Array<Line> lines) {
		for (int i = 0; i < lines.GetSize(); i++) {
			for (int j = 0; j < lines[i].GetNumCoordinates(); j++) {
				pattern[lines[i][j].y][lines[i][j].x] = 1;
			}
		}
	}

	bool IsInsidePolygon(Coordinate c) {
		int edges = 0;
		for (int w = c.x; w < width; w += 1) {
			if (pattern[c.y][w] == 1) {
				edges += 1;
			}
		}
		return (((edges % 2) == 1));
	}

	bool IsEdge(Array<Line>& lines, Coordinate c) {
		for (int i = 0; i < lines.GetSize(); i++) {
			for (int j = 0; j < lines[i].GetNumCoordinates(); j++) {
				if (lines[i][j] == c) {
					return true;
				}
			}
		}
		return false;
	}

	void FillInPolygon(Array<Line>& lines) {
		Coordinate c;
		Array<Coordinate> edgePoints;
		for (int i = 0; i < lines.GetSize(); i++) {
			for (int j = 0; j < lines[i].GetNumCoordinates(); j++) {
				if (!edgePoints.exists(lines[i][j])) {
					edgePoints.Add(lines[i][j]);
				}
			}
		}

		int* edgeCount = new int[height];
		for (int i = 0; i < height; i++) {
			edgeCount[i] = 0;
		}

		for (int i = 0; i < edgePoints.GetSize(); i++) {
			edgeCount[edgePoints[i].y] += 1;
		}

		for (int h = 0; h < height; h++) {
			if ((edgeCount[h] > 1)) {
				for (int w = 0; w < width; w++) {
					c.y = h;
					c.x = w;
					if (IsInsidePolygon(c)) {
						pattern[h][w] = 1;
					}
				}
			}
		}

		delete[] edgeCount;
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
		/* Plot pattern */
		for (int i = startHeight; i <= endHeight; i++) {
			for (int j = startWidth; j <= endWidth; j++) {
				pattern[i][j] = 1;
			}
		}
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
		Coordinate highestPoint(centerHeight + heightRadius, centerWidth);
		Coordinate lowestPoint(centerHeight - heightRadius, centerWidth);
		Coordinate leftMostPoint(centerHeight, centerWidth - widthRadius);
		Coordinate rightMostPoint(centerHeight, centerWidth + widthRadius);

		Array<Line> arr;
		arr.Add(Line(highestPoint, leftMostPoint));
		arr.Add(Line(highestPoint, rightMostPoint));
		arr.Add(Line(lowestPoint, rightMostPoint));
		arr.Add(Line(lowestPoint, leftMostPoint));

		PlotLines(arr);
		FillInPolygon(arr);
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
{	string fName = "C:\\Users\\james\\Code\\CPP\\MachineLearningCPP\\MachineLearningCPP\\PatternGenerator\\Output\\Outfile.txt";

	Pattern P;

	ofstream outFile;
	outFile.open(fName);

	UnitPattern *sq = new DiamondPattern(100, 100, 0.15, 0.85);
	sq->PrintPattern(outFile);
	P = Pattern(950, 950 , 50, 0, sq);
	//P.PrintPattern(outFile);
	delete sq;

	P.SavePatternToBmp("C:\\Users\\james\\Code\\CPP\\MachineLearningCPP\\MachineLearningCPP\\PatternGenerator\\Output\\Outfile.bmp");


	outFile.close();


	return 0;
}
