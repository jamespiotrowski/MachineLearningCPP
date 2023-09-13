#include <fstream>
#include "Pattern.h"

using namespace std;

/*

	- SOME UNIT PATTERNS HAVE MULTIPLE SCALES
		- To generate every pattern with every possible combination of scales is unreasonable. 

	- After further thought, scales add a lot of complexity... here are future considerations:
		- Rotations for certain shapes
		- scale options (all different possible scales)

	- What I think we will do:
		- All shapes have 1 or 2 scales
		- We will generate images with different versions of the same unit pattern in the same picture. 

*/


class PatternGenerator {

private:
	Array<PatternType> patternList;
	int totalImages = 0;

	int unitPatternWidth = 50;
	int unitPatternHeight = 50;

	int patternWidth = 50;
	int patternHeight = 50;

	double minScale = 0;
	double scaleStep = 0;
	double maxScale = 0;

	int allowedNumberOfScales = 0;

	const int minPixelsAllowed = 20;
	const int maximumPossibleScales = 8;

	bool clipping = false;
	bool center = true;

	double percentageOfPatternsToKeep = 0.01;

	string outputDirectory = "";

	Array<Array<UnitPattern*>> unitPatterns;
	Array<Array<int>> unitPatternIndexes;

	void deallocateAllUnitPattens() {
		for (int i = 0; i < unitPatterns.GetSize(); i++) {
			for (int j = 0; j < unitPatterns[i].GetSize(); j++) {
				if (unitPatterns[i][j] != nullptr) {
					delete unitPatterns[i][j];
					unitPatterns[i][j] = nullptr;
				}
			}
		}
		unitPatterns.reset();
		unitPatternIndexes.reset();
	}
	
	void cleanAndStandardizeMembers(bool smartScaleDetection) {


		/* The unit pattern must always be smaller than the big pattern */
		if (unitPatternWidth > patternWidth) { patternWidth = unitPatternWidth; }
		if (unitPatternHeight > patternHeight) { patternHeight = unitPatternHeight; }
		/* No dupes */
		patternList.removeDuplicates();

		/* Scales */
		double smallestDimension = (unitPatternWidth < unitPatternHeight) ? unitPatternWidth : unitPatternHeight;
		if (smartScaleDetection) {
			scaleStep = 1.0 / smallestDimension;
			minScale = 0.2;
			maxScale = 0.9;
		}

		/*
			Scale is between 0 and 1
			When the scale changes, there needs to be a visible change in the image

			for example:
				- A scale step of 0.01 is acceptable for 100 pixels because ((0.01) * (100)) = 1 pixel
				- A scale step of 0.01 is NOT acceptable for 50 pixels because ((0.01) * (50)) = 0 pixels

			So, we take 1 / (minimum dimension) and if the scale step is smaller than that, we correct. 
		*/
		if ((1.0 / smallestDimension) > scaleStep) {
			scaleStep = 1.0 / smallestDimension;
		}

		/* we should not make any shapes where the lesser dimension of the shape will be just a few pixels */
		// Anything dimension that is less than 20 pixels might be too "grainy" to recognize
		if ((minScale * smallestDimension) < minPixelsAllowed) {
			minScale = (double)minPixelsAllowed / smallestDimension;
		}

		/* we should allow at least 10 pixels to border the picture */
		int allowedBorder = (minPixelsAllowed / 2);
		double border = (smallestDimension)-(smallestDimension * maxScale);
		if (border < allowedBorder) {
			maxScale = (smallestDimension - allowedBorder) / smallestDimension;
		}

	}

	UnitPattern* GetUnitPattern(PatternType p, double* scaleSet) {
		switch (p) {
		case(SQUARE): { return new SquarePattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(RECTANGLE): { return new RectanglePattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(DIAMOND): { return new DiamondPattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(TRIANGLE): { return new TrianglePattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(HORIZONTAL_STRIPES): { return new HorizontalStripePattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(VERTICAL_STRIPES): { return new VerticalStripePattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(CIRCLE): { return new CirclePattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(HEXAGON): { return new HexagonPattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(PENTAGON): { return new PentagonPattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(HEPTAGON): { return new HeptagonPattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(STAR): { return new StarPattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		case(OCTAGON): { return new OctogonPattern(unitPatternHeight, unitPatternWidth, scaleSet); }
		}
		return nullptr;
	}

	int GetNumberOfUnitPatternsPerPattern(const int &verticalOffset, const int &horizontalOffset) const {

		double tUnitHeight = unitPatternHeight + verticalOffset;
		double tUnitWidth = unitPatternWidth + horizontalOffset;

		double tHeightUnits = (double)patternHeight / tUnitHeight;
		double tWidthUnits = (double)patternWidth / tUnitWidth;

		if (clipping) {
			tHeightUnits = ceil(tHeightUnits);
			tWidthUnits = ceil(tWidthUnits);
		}
		else {
			tHeightUnits = floor(tHeightUnits);
			tWidthUnits = floor(tWidthUnits);
		}
		
		return (int)(tHeightUnits * tWidthUnits);
	}

	Array<Array<int>> GetPatternCombinations(int pattern, const int& verticalOffset, const int& horizontalOffset) const {
		int totalUnitsPerPattern = GetNumberOfUnitPatternsPerPattern(verticalOffset, horizontalOffset);
		return SomeCombinations(unitPatternIndexes.At(pattern), totalUnitsPerPattern, percentageOfPatternsToKeep);
	}

	Pattern GetPattern(int pattern, const int& verticalOffset, const int& horizontalOffset, const Array<int>& combination) const {
		//string s = "";
		//for (int j = 0; j < combination.GetSize(); j++) {
		//	s = s + to_string(combination.At(j)) + " ";
		//}
		//cout << s << endl;
		return Pattern(patternHeight, patternWidth, verticalOffset, horizontalOffset, clipping, center, unitPatterns.At(pattern), combination);
	}

	/*******************************************

	*******************************************/

	void GenerateAllUnitPatterns() {

		deallocateAllUnitPattens();

		Array<double> scales;
		for (double s = minScale; s < maxScale; s += scaleStep) {
			scales.Add(s);
		}

		for (int p = 0; p < patternList.GetSize(); p++) {
			Array<UnitPattern*> unitPatternSet;
			unitPatterns.Add(unitPatternSet);
		}

		double* scaleForPattern = new double[maximumPossibleScales];
		for (int i = 0; i < scales.GetSize(); i++) {
			// All scales are the same here - for now
			for (int j = 0; j < maximumPossibleScales; j++) {
				scaleForPattern[j] = scales[i];
			}

			for (int p = 0; p < patternList.GetSize(); p++) {
				unitPatterns[p].Add(GetUnitPattern(patternList[p], scaleForPattern));
			}
		}

		for (int i = 0; i < unitPatterns.GetSize(); i++) {
			unitPatternIndexes.Add(Array<int>());
			for (int j = 0; j < unitPatterns[i].GetSize(); j++) {
				unitPatternIndexes[i].Add(j);
			}
		}

		delete[] scaleForPattern;
	}

public:

	PatternGenerator(
		Array<PatternType> patternList
		, int unitPatternWidth
		, int unitPatternHeight
		, int patternWidth
		, int patternHeight
		, double minScale
		, double scaleStep
		, double maxScale
		, int allowedNumberOfScales
		, string outputDirectory
		, bool clipping = false
		, bool center = true
		, double percentageOfPatternsToKeep = 0.01
	) : patternList(patternList), unitPatternWidth(unitPatternWidth), unitPatternHeight(unitPatternHeight)
		, patternWidth(patternWidth), patternHeight(patternHeight), minScale(minScale), scaleStep(scaleStep)
		, maxScale(maxScale), outputDirectory(outputDirectory), allowedNumberOfScales(allowedNumberOfScales)
		, clipping(clipping), center(center), percentageOfPatternsToKeep(percentageOfPatternsToKeep)
	{
		cleanAndStandardizeMembers(false);
		allowedNumberOfScales = 1; // Not going to incorporate multiple scales just yet.
		GenerateAllUnitPatterns();
	}

	~PatternGenerator() {
		deallocateAllUnitPattens();
	}

	void MakePatterns(/* Will add stuff in here to control the offsets */) const {
		int verticalOffset = 25;
		int horizontalOffset = 25;
		for (int currentPattern = 0; currentPattern < patternList.GetSize(); currentPattern++) {
			Array<Array<int>> allCombinations = GetPatternCombinations(currentPattern, verticalOffset, horizontalOffset);
			//cout << allCombinations.GetSize() << endl;
			int tImgs = 0;
			string outputFile;
			string currentPatternString = GetNameForPattern(patternList.At(currentPattern));
			for (int i = 0; i < allCombinations.GetSize(); i++) {
				/*string s = "";
				for (int j = 0; j < allCombinations[i].GetSize(); j++) {
					s = s + to_string(allCombinations[i][j]) + " ";
				}
				cout << s << endl;*/
				Pattern p = GetPattern(currentPattern, verticalOffset, horizontalOffset, allCombinations[i]);
				outputFile = outputDirectory + currentPatternString + "_" + to_string(tImgs) + ".bmp";
				p.SavePatternToBmp(outputFile);
				tImgs += 1;
			}
			//cout << "Done!" << endl;
		}
	}

	void SaveUnitPatternPNGs() {
		string outputFile;
		int tImgs = 0;
		for (int p = 0; p < patternList.GetSize(); p++) {
			string currentPattern = GetNameForPattern(patternList[p]);
			for (int i = 0; i < unitPatterns[p].GetSize(); i++) {
				UnitPattern* up = unitPatterns[p][i];
				outputFile = outputDirectory + currentPattern + "_" + to_string(tImgs) + ".bmp";
				Pattern p = Pattern(unitPatternHeight, unitPatternWidth, 0, 0, clipping, center, up); 
				p.SavePatternToBmp(outputFile);
				tImgs += 1;
			}
		}
	}
};


/* Possible options:

	- Size of input
		- Size of units
		- Size of larger pattern
		- (or) Smart Scale detection

	- Scales
		- Use all scales
		- Minimum Scale
		- Max Scale
		- Scale Step


	- Balance data set

	- Patterns to include

*/

int main()
{	
	/*
	Array<int> a;
	a.Add(1);
	a.Add(2);
	a.Add(3);
	a.Add(4);

	Array<Array<int>> c = AllCombinations(a, 3);

	for (int i = 0; i < c.GetSize(); i++) {
		string s = "";
		for (int j = 0; j < c[i].GetSize(); j++) {
			s += to_string(c[i][j]) + " ";
		}
		cout << s << endl;
	}
	*/
	Array<PatternType> patternList;
	patternList.Add(SQUARE);
	//patternList.Add(RECTANGLE);	
	patternList.Add(TRIANGLE);
	patternList.Add(PENTAGON);
	patternList.Add(STAR);
	//patternList.Add(VERTICAL_STRIPES);
	//patternList.Add(HORIZONTAL_STRIPES);
	patternList.Add(CIRCLE);
	patternList.Add(DIAMOND);
	patternList.Add(HEXAGON);
	patternList.Add(OCTAGON);
	patternList.Add(HEPTAGON);

	PatternGenerator pg(
		patternList
		, 100	// Unit Pattern height
		, 100	// Unit Pattern width
		, 300	// Total Pattern height
		, 300	// Total Pattern width
		, 0.2	// Starting scale
		, 0.1	// Scale step
		, 0.9	// Ending scale
		, 1		// Scales allowed
		, "C:\\Users\\james\\Code\\CPP\\MachineLearningCPP\\MachineLearningCPP\\PatternGenerator\\Output\\" // Output Directory
		, false	// clipping
		, true	// center
		, 0.01	// percentage of combinations to keep
	);

	pg.MakePatterns();

	//pg.SaveUnitPatternPNGs();

	/*
	double scaleSet[1] = { 0.8 };
	UnitPattern *up = new SquarePattern(100, 100, scaleSet);
	Pattern P(2000, 2000, 0, 0, true, true, up);
	P.SavePatternToBmp("C:\\Users\\james\\Code\\CPP\\MachineLearningCPP\\MachineLearningCPP\\PatternGenerator\\Output\\output.png");
	delete up;
	*/

	return 0;
}
