#include <iostream>

using namespace std;

enum PatternType {
	SQUARE
	, CIRCLE 
};

class UnitPattern {
private:

	static void DetermineHeightAndWidthWithTrueCenter(unsigned int& height, unsigned int& width) {
		if ((height % 2) == 0) {
			height += 1;
		}

		if ((width % 2) == 0) {
			width += 1;
		}
	}

protected:

	char** pattern = nullptr;
	unsigned int height = 0, width = 0;


	void GetCenter(const unsigned int& height, const unsigned int& width, unsigned int& centerHeight, unsigned int& centerWidth) {
		centerHeight = height / 2;
		centerWidth = width / 2;
	}

public:

	UnitPattern(unsigned int height, unsigned int width) : height(height), width(width){
		DetermineHeightAndWidthWithTrueCenter(this->height, this->width);
		pattern = new char* [height];
		for (unsigned int i = 0; i < height; i++) {
			pattern[i] = new char[width];
		}
	}

	virtual void GeneratePattern() = 0;

};

class SquarePattern : public UnitPattern {
private:

	double scale = 0.5;

public:
	SquarePattern(unsigned int height, unsigned int width, double scale) : UnitPattern(height, width) {
		GeneratePattern();
	}

	void GeneratePattern() {
		unsigned int centerHeight = 0, centerWidth = 0;
		GetCenter(height, width, centerHeight, centerWidth);
	}

};


int main()
{

	return 0;
}
