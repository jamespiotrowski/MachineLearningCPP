#pragma once

#include "Polygon.h"
#include "BitMap.h"

using namespace std;

enum PatternType
{
	SQUARE
	, RECTANGLE
	, DIAMOND
	, TRIANGLE
	, HORIZONTAL_STRIPES
	, VERTICAL_STRIPES
	, CIRCLE
	, HEXAGON
	, PENTAGON
	, HEPTAGON
	, STAR
	, OCTAGON
};

int GetScalesForPattern(const PatternType &pt) {
	switch (pt) {
	case(SQUARE): case(HORIZONTAL_STRIPES): case(VERTICAL_STRIPES): case(CIRCLE): case(STAR) : { return 1;}
	case(RECTANGLE): case(DIAMOND): { return 2; }
	case(TRIANGLE): { return 3; }
	case(PENTAGON): { return 5; }
	case(HEXAGON): { return 6; }
	case(HEPTAGON): { return 7; }
	case(OCTAGON): { return 8; }
	default: { return 0; }
	}
	return 0;
}

string GetNameForPattern(const PatternType& pt) {
	switch (pt) {
	case(SQUARE): { return "Square"; }
	case(HORIZONTAL_STRIPES): { return "HorizontalStripe"; }
	case(VERTICAL_STRIPES): { return "VerticalStripe"; }
	case(CIRCLE): { return "Circle"; }
	case(STAR): { return "Star"; }
	case(RECTANGLE): { return "Rectangle"; }
	case(DIAMOND): { return "Diamond"; }
	case(TRIANGLE): { return "Triangle"; }
	case(PENTAGON): { return "Pentagon"; }
	case(HEXAGON): { return "Hexagon"; }
	case(HEPTAGON): { return "Heptagon"; }
	case(OCTAGON): { return "Octogon"; }
	default: { return ""; }
	}
	return "";
}

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class UnitPattern {
public:

	static void DetermineHeightAndWidthWithTrueCenter(int& height, int& width) {
		if ((height % 2) == 0) {
			height += 1;
		}

		if ((width % 2) == 0) {
			width += 1;
		}
	}

	static void GetCenter(const int& height, const int& width, int& centerHeight, int& centerWidth) {
		centerHeight = height / 2;
		centerWidth = width / 2;
	}

protected:
	char** pattern = nullptr;
	int height = 0, width = 0, numberOfScales = 0;
	double* scales = nullptr;
	bool verticalOffsetAllowed = true;
	bool horizontalOffsetAllowed = true;
	PatternType patternType;

private:

	// possible optimization..
	void FillInUntilEdge(const int& y, int& w) {
		while (w < width && pattern[y][w] != 1) {
			pattern[y][w] = 1;
			w += 1;
		}
	}

protected:

	void SetScale(int numberOfScales, const double* scales) {
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

		for (int i = 0; i < loneEdgePoints.GetSize(); i++) {
			if (loneEdgePoints.At(i).y == c.y && loneEdgePoints.At(i).x >= c.x) {
				edges -= 1;
			}
		}

		return (((edges % 2) == 1));
	}

	// For a true polygon
	void FillInPolygon(const Polygon& p) {
		Coordinate c;
		for (int h = 0; h < height; h++) {
			c.y = h;
			for (int w = 0; w < width; w++) {
				c.x = w;
				if (p.isInsidePolygon(c)) {
					pattern[h][w] = 1; // Could optimize by using fillInUntilEdge, but not general and may not be necesarry
				}
			}
		}
	}

	// for a not true polygon, though... a circle could just be a bunch of edges. we will see
	void FillInPolygon(const Array<Coordinate>& loneEdgePoints) {
		Coordinate c;
		for (int h = 0; h < height; h++) {
			c.y = h;
			for (int w = 0; w < width; w++) {
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

	UnitPattern(int height, int width, PatternType patternType) : height(height), width(width), patternType(patternType) {
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
	int GetScales() const { return numberOfScales; }
	PatternType GetPatternType() const { return patternType; }
	bool allowsVerticalOffset() const { return verticalOffsetAllowed; }
	bool allowsHorizontalOffset() const { return horizontalOffsetAllowed; }

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

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

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

	SquarePattern(int height, int width, const double* s) : UnitPattern(height, width, SQUARE) {
		SetScale(GetScalesForPattern(patternType), s);
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

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

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

	DiamondPattern(int height, int width, const double* s) : UnitPattern(height, width, DIAMOND) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
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

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

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
		for (int i = 0; i < height; i++) {
			edgePointCounter.Add(0);
		}

		for (int i = 0; i < circlePoints.GetSize(); i++) {
			edgePointCounter[circlePoints[i].y] += 1;
			pattern[circlePoints[i].y][circlePoints[i].x] = 1;
		}

		int minHeight = 0, maxHeight = 0;
		for (int i = 0; i < edgePointCounter.GetSize(); i++) {
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
		for (int j = 0; j < circlePoints.GetSize(); j++) {

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

		FillInPolygon(loneEdgePoints);
	}

public:

	CirclePattern(int height, int width, const double* s) : UnitPattern(height, width, CIRCLE) {
		SetScale(GetScalesForPattern(patternType), s);
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
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x); // Detmine center coordinates
		
		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

		double t = ((height < width) ? height : width) / 2;
		//t = t * t;
		//t = t + t;
		//t = sqrt(t);

		int sides = 3;

		double split = 360.0 / sides;

		//cout << t << endl;
		//cout << scales[0] * t << endl;
		//cout << scales[0] << endl;

		Array<Coordinate> c;
		for (int i = 0; i < sides; i++) {
			c.Add(Polygon::ComputePointGivenAngleAndDistance(180.0 + (i * split), scales[i] * t, centerCoord));
		}
		
		Array<Edge> edges;
		for (int i = 0; i < sides; i++) {
			edges.Add(Edge(c[i], c[(i+1) % sides]));
		}

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	TrianglePattern(int height, int width, const double* s) : UnitPattern(height, width, TRIANGLE) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
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
		GetCenter(height, width, centerCoord.y, centerCoord.x); // Detmine center coordinates

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

		double t = ((height < width) ? height : width) / 2;
		//t = t * t;
		//t = t + t;
		//t = sqrt(t);

		int sides = 5;

		double split = 360.0 / sides;

		Array<Coordinate> c;
		for (int i = 0; i < sides; i++) {
			c.Add(Polygon::ComputePointGivenAngleAndDistance(180.0 + (i * split), scales[i] * t, centerCoord));
		}

		Array<Edge> edges;
		for (int i = 0; i < sides; i++) {
			edges.Add(Edge(c[i], c[(i + 1) % sides]));
		}

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	PentagonPattern(int height, int width, const double* s) : UnitPattern(height, width, PENTAGON) {
		SetScale(GetScalesForPattern(patternType), s);
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

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;


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

	}

public:

	StarPattern(int height, int width, const double* s) : UnitPattern(height, width, STAR) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class HorizontalStripePattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x);					 // Detmine center coordinates
		int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = false;


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

	HorizontalStripePattern(int height, int width, const double* s) : UnitPattern(height, width, HORIZONTAL_STRIPES) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class VerticalStripePattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x);					 // Detmine center coordinates
		int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = false;


		/* Compute the pattern area */
		Coordinate* c = new Coordinate[4];
		c[0] = Coordinate(0, centerCoord.x + radius);
		c[1] = Coordinate(0, centerCoord.x - radius);
		c[2] = Coordinate(width - 1, centerCoord.x - radius);
		c[3] = Coordinate(width - 1, centerCoord.x + radius);


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

	VerticalStripePattern(int height, int width, const double* s) : UnitPattern(height, width, VERTICAL_STRIPES) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
	}

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class RectanglePattern : public UnitPattern {
private:
	const int widthScale = 0;
	const int heightScale = 1;

	void GenerateUnitPattern() {
		int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth); // Detmine center coordinates
		int widthRadius = (scales[widthScale] * width) / 2;
		int heightRadius = (scales[heightScale] * height) / 2;

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

		/* Compute the pattern area */
		int startHeight = centerHeight - heightRadius;
		int startWidth = centerWidth - widthRadius;
		int endHeight = centerHeight + heightRadius;
		int endWidth = centerWidth + widthRadius;

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

	RectanglePattern(int height, int width, const double* s) : UnitPattern(height, width, RECTANGLE) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class HexagonPattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x); // Detmine center coordinates

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

		double t = ((height < width) ? height : width) / 2;
		//t = t * t;
		//t = t + t;
		//t = sqrt(t);

		int sides = 6;

		double split = 360.0 / sides;

		Array<Coordinate> c;
		for (int i = 0; i < sides; i++) {
			c.Add(Polygon::ComputePointGivenAngleAndDistance(180.0 + (i * split), scales[i] * t, centerCoord));
		}

		Array<Edge> edges;
		for (int i = 0; i < sides; i++) {
			edges.Add(Edge(c[i], c[(i + 1) % sides]));
		}

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	HexagonPattern(int height, int width, const double* s) : UnitPattern(height, width, HEXAGON) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class HeptagonPattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x); // Detmine center coordinates

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

		double t = ((height < width) ? height : width) / 2;
		//t = t * t;
		//t = t + t;
		//t = sqrt(t);

		int sides = 7;

		double split = 360.0 / sides;

		Array<Coordinate> c;
		for (int i = 0; i < sides; i++) {
			c.Add(Polygon::ComputePointGivenAngleAndDistance(180.0 + (i * split), scales[i] * t, centerCoord));
		}

		Array<Edge> edges;
		for (int i = 0; i < sides; i++) {
			edges.Add(Edge(c[i], c[(i + 1) % sides]));
		}

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	HeptagonPattern(int height, int width, const double* s) : UnitPattern(height, width, HEPTAGON) {
		SetScale(GetScalesForPattern(patternType), s);
		GenerateUnitPattern();
	}

};


/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class OctogonPattern : public UnitPattern {
private:
	const int scale1 = 0; // Scale that square uses

	void GenerateUnitPattern() {
		Coordinate centerCoord;
		GetCenter(height, width, centerCoord.y, centerCoord.x); // Detmine center coordinates

		verticalOffsetAllowed = true;
		horizontalOffsetAllowed = true;

		double t = ((height < width) ? height : width) / 2;
		//t = t * t;
		//t = t + t;
		//t = sqrt(t);

		int sides = 8;

		double split = 360.0 / sides;

		Array<Coordinate> c;
		for (int i = 0; i < sides; i++) {
			c.Add(Polygon::ComputePointGivenAngleAndDistance(180.0 + (i * split), scales[i] * t, centerCoord));
		}

		Array<Edge> edges;
		for (int i = 0; i < sides; i++) {
			edges.Add(Edge(c[i], c[(i + 1) % sides]));
		}

		Polygon p(edges);
		p.plotPolygon(pattern);
		FillInPolygon(p);
	}

public:

	OctogonPattern(int height, int width, const double* s) : UnitPattern(height, width, OCTAGON) {
		SetScale(GetScalesForPattern(patternType), s);
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
	bool clipping = false;
	bool centerPattern = true;

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

	Pattern(int height, int width, int verticalOffset, int horizontalOffset, bool clipping, bool center, const Array<UnitPattern*>& unitPatterns, const Array<int> &patternSet)
		: verticalOffset(verticalOffset), horizontalOffset(horizontalOffset)
		, height(height), width(width)
		, centerPattern(center), clipping(clipping)
	{
		UnitPattern::DetermineHeightAndWidthWithTrueCenter(height, width);

		canvas = new char* [height];
		for (unsigned h = 0; h < height; h++) {
			canvas[h] = new char[width];
			for (int w = 0; w < width; w++) {
				canvas[h][w] = 0;
			}
		}
		GeneratePattern(unitPatterns, patternSet);
	}

	Pattern(int height, int width, int verticalOffset, int horizontalOffset, bool clipping, bool center, UnitPattern* unitPattern) 
		: verticalOffset(verticalOffset), horizontalOffset(horizontalOffset)
		, height(height), width(width) 
		, centerPattern(center), clipping(clipping)
	{
		UnitPattern::DetermineHeightAndWidthWithTrueCenter(height, width);
		canvas = new char* [height];
		for (unsigned h = 0; h < height; h++) {
			canvas[h] = new char[width];
			for (int w = 0; w < width; w++) {
				canvas[h][w] = 0;
			}
		}

		Array<int> pattern;
		pattern.Add(0);
		
		Array<UnitPattern*> up;
		up.Add(unitPattern);

		GeneratePattern(up, pattern);
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
		clipping = copy.clipping;
		centerPattern = copy.centerPattern;
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

	void GeneratePattern(const Array<UnitPattern*>& unitPatterns, const Array<int>& patternSet) {
		if (unitPatterns.GetSize() == 0) {
			return;
		}

		// Unit Pattern Dimensions
		int unitHeight = unitPatterns.At(0)->GetHeight();
		int unitWidth = unitPatterns.At(0)->GetWidth();

		// Total space taken up by a unit and it's offset
		double tHeightSpace = (unitHeight + verticalOffset);
		double tWidthSpace = (unitWidth + horizontalOffset);

		// The number of total (unit + offset) that can fit in the pattern space
		double tHeightUnits;
		double tWidthUnits;

		if (clipping) {
			tHeightUnits = ceil((double)height / tHeightSpace);
			tWidthUnits = ceil((double)width / tWidthSpace);
		}
		else {
			tHeightUnits = floor((double)height / tHeightSpace);
			tWidthUnits = floor((double)width / tWidthSpace);
		}

		// The space needed to FULLY fit everything
		double virtualHeight = tHeightUnits * tHeightSpace;
		double virtualWidth = tWidthUnits * tWidthSpace;

		// Where the pattern begins
		int startHeight = 0;
		int startWidth = 0;

		// Where the pattern ends
		int endHeight = tHeightUnits * tHeightSpace;
		int endWidth = tWidthUnits * tWidthSpace;

		if (centerPattern) {

			double heightDiff = virtualHeight - height;
			double widthDiff = virtualWidth - width;
	
			startHeight = (((-1.0 * heightDiff) / 2.0) + ((double)verticalOffset / 2.0));
			startWidth = (((-1.0 * widthDiff) / 2.0) + ((double)horizontalOffset / 2.0));

			endHeight = height + ((heightDiff / 2.0) - ((double)verticalOffset / 2.0));
			endWidth = width + ((widthDiff / 2.0) - ((double)horizontalOffset / 2.0));

			//cout << startHeight << "," << endHeight << " | " << startWidth << "," << endWidth << endl;

		}

		int heightStep = unitHeight + verticalOffset;	// The amount of height steps to take whilst pasting the unit pattern
		int widthStep = unitWidth + horizontalOffset;	// The amount of width steps to take whilst pasting the unit pattern
		int innerHeightLimit = 0;	// Used inside the loop to limit the inner for loops
		int innerWidthLimit = 0;	// Used inside the loop to limit the inner for loops

		int slider = 0;

		// Adjust this only to plot within realm
		for (int oH = startHeight; oH < endHeight; oH += heightStep) {
			innerHeightLimit = oH + unitHeight;
			innerHeightLimit = (innerHeightLimit > height) ? height : innerHeightLimit;
			for (int oW = startWidth; oW < endWidth; oW += widthStep) {
				innerWidthLimit = oW + unitWidth;
				innerWidthLimit = (innerWidthLimit > width) ? width : innerWidthLimit;
				for (int iH = ((oH >= 0) ? oH : 0); iH < innerHeightLimit; iH += 1) {
					for (int iW = ((oW >= 0) ? oW : 0); iW < innerWidthLimit; iW += 1) {
						canvas[iH][iW] = unitPatterns.At(patternSet.At(slider))->At(iH - oH, iW - oW);
					}
				}
				slider = ((slider + 1) % (unitPatterns.GetSize()));
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