#include <iostream>

using namespace std;

#define SCALES_SQUARE 1

enum PatternType {
	SQUARE
	, CIRCLE 
};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class UnitPattern {
protected:

	char** pattern = nullptr;
	unsigned int height = 0, width = 0, numberOfScales = 0;
	double* scales = nullptr;

	void GetCenter(const unsigned int& height, const unsigned int& width, unsigned int& centerHeight, unsigned int& centerWidth) {
		centerHeight = height / 2;
		centerWidth = width / 2;
	}

	void SetScale(unsigned int numberOfScales, double* scales) {
		this->scales = new double[numberOfScales];
		for (unsigned int i = 0; i < numberOfScales; i++) {
			this->scales[i] = scales[i];
		}
	}

private:

	static void DetermineHeightAndWidthWithTrueCenter(unsigned int& height, unsigned int& width) {
		if ((height % 2) == 0) {
			height += 1;
		}

		if ((width % 2) == 0) {
			width += 1;
		}
	}

	void clear() {
		for (unsigned int i = 0; i < height; i++) {
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

	UnitPattern(unsigned int height, unsigned int width) : height(height), width(width){
		DetermineHeightAndWidthWithTrueCenter(this->height, this->width);
		pattern = new char* [this->height];
		for (unsigned int i = 0; i < this->height; i++) {
			pattern[i] = new char[this->width];
			for (unsigned int j = 0; j < this->width; j++) {
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
		for (unsigned int i = 0; i < height; i++) {
			pattern[i] = new char[width];
			for (unsigned int j = 0; j < width; j++) {
				pattern[i][j] = copy.pattern[i][j];
			}
		}

		this->scales = new double[numberOfScales];
		for (unsigned int i = 0; i < numberOfScales; i++) {
			scales[i] = copy.scales[i];
		}
	}

	UnitPattern(const UnitPattern& copy) {
		(*this) = copy;
	}

	void PrintPattern(ostream& out) {
		string s = "";
		for (unsigned int i = 0; i < height; i++) {
			for (unsigned int j = 0; j < width; j++) {
				s += ToSymbol(pattern[i][j]);
			}
			s += "\n";
		}
		s += "\n";

		out << s;
	}

	virtual void GeneratePattern() = 0;

};

/**********************************************************************************************
###############################################################################################
#####
###############################################################################################
**********************************************************************************************/
class SquarePattern : public UnitPattern {
private:
	const unsigned int scale1 = 0; // Scale that square uses

public:

	SquarePattern(unsigned int height, unsigned int width, double scale) : UnitPattern(height, width) {
		SetScale(SCALES_SQUARE, &scale);
		GeneratePattern();
	}

	void GeneratePattern() {
		unsigned int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth);					 // Detmine center coordinates
		unsigned int radius = (scales[scale1] * ((height < width) ? height : width)) / 2; // Determine squares "radius" using minimum dimention of unit
		/* Compute the pattern area */
		unsigned int startHeight = centerHeight - radius;
		unsigned int startWidth = centerWidth - radius;
		unsigned int endHeight = centerHeight + radius;
		unsigned int endWidth = centerWidth + radius;
		/* Plot pattern */
		for (unsigned int i = startHeight; i <= endHeight; i++) {
			for (unsigned int j = startWidth; j <= endWidth; j++) {
				pattern[i][j] = 1;
			}
		}
	}
};


int main()
{
	UnitPattern* sq = new SquarePattern(28, 28, 0.5);
	sq->PrintPattern(cout);
	delete sq;

	sq = new SquarePattern(28, 28, 0.25);
	sq->PrintPattern(cout);
	delete sq;

	sq = new SquarePattern(28, 28, 0.75);
	sq->PrintPattern(cout);
	delete sq;

	sq = new SquarePattern(28, 28, 1.0);
	sq->PrintPattern(cout);
	delete sq;

	return 0;
}
