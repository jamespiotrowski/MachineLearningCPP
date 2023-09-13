#pragma once

#include "Combination.h"
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

#define M_PI 3.141592653589793;

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/

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

	Angle() : c(0, 0), intersect(false) { }

	Angle(Coordinate c1, bool i) : c(c1) {
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

	static Coordinate GetSharedPoint(const Edge& edge1, const Edge& edge2) {
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

	double getValueAtY(int y) const {
		return ((infSlope) ? c2.x : ((((double)y) - b) / slope));
	}
};

class Polygon {
private:
	Array<Coordinate> points;
	Array<Edge> edges;
	Array<Angle> angles;

	int maxY, minY, maxX, minX;

public:

	bool pointWithinPolygonRange(const Coordinate& c) const {
		return c.y >= minY && c.y <= maxY && c.x <= maxX && c.x >= minX;
	}

	static Array<Coordinate> ComputeStraitLine(Coordinate c1, Coordinate c2) {
		Array<Coordinate> line;
		if (c1.y == c2.y) {
			int greaterX = (c2.x > c1.x) ? c2.x : c1.x;
			int lesserX = (c2.x > c1.x) ? c1.x : c2.x;
			for (int i = lesserX + 1; i < greaterX; i++) {
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
			int intersections = 0;
			for (int i = 0; i < edges.GetSize(); i++) {
				if (intersectsEdge(c, edges.At(i))) { // If ray would intersect an edge
					intersections += 1;
				}
			}

			for (int i = 0; i < angles.GetSize(); i++) {
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
		for (int i = 0; i < edges.GetSize(); i++) {
			Array<Coordinate> line = ComputeStraitLine(edges[i].c1, edges[i].c2);
			for (int j = 0; j < line.GetSize(); j++) {
				points.Add(line[j]);
			}
		}
		// Remove dupes
		points.removeDuplicates();

		// Find angles
		for (int i = 0; i < edges.GetSize(); i++) {
			for (int j = 0; j < edges.GetSize(); j++) {
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
			for (int i = 1; i < points.GetSize(); i++) {
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
		for (int i = 0; i < points.GetSize(); i++) {
			grid[points[i].y][points[i].x] = 1;
		}
	}

	void printPolygon() {
		string p = "Points : { ";
		for (int i = 0; i < points.GetSize(); i++) {
			p += points[i].ToString() + " ";
		}
		p += " }\n";

		string e = "Edges : { ";
		for (int i = 0; i < edges.GetSize(); i++) {
			e += edges[i].ToString() + " ";
		}
		e += " }\n";

		string a = "Angles : { ";
		for (int i = 0; i < angles.GetSize(); i++) {
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
