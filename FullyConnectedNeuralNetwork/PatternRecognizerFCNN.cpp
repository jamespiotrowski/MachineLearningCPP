#include "FCNN.h"
#include "Array.h"
#include <fstream>
using namespace std;

void PrintArray(ostream& out, const double arr[], const unsigned int& arrSize) {
    string s = "[";
    for (unsigned int i = 0; i < arrSize; i++) {
        s += to_string(arr[i]);
        if (i < (arrSize - 1)) {
            s += ",";
        }
    }
    s += "]";
    out << s << endl;
}

struct Pattern {
    static Array<string> Classes;
    int LableId = 0;
    string LableName = "";
    int height = 0, width = 0, totalInput = 0;
    double* data = nullptr;

    double operator[](const int& i) const { return (double)data[i]; }

    Pattern() { }
    ~Pattern() { if (data != nullptr) { delete[] data; } }

    Pattern(const string& s) {
        int currentIndex = 0;
        int commasFound = 0;
        string currentString = "";

        // Get data dimensions
        while (currentIndex < s.size() && commasFound < 3) {
            if (s[currentIndex] == ',') {
                commasFound += 1;
                switch (commasFound) {
                case(1): { LableName = currentString; break; }
                case(2): { height = stoi(currentString); break; }
                case(3): { width = stoi(currentString); break; }
                }
                currentString = "";
            }
            else {
                currentString += s[currentIndex];
            }
            currentIndex += 1;
        }

        totalInput = height * width;
        data = new double[totalInput];

        // Get the data string
        for (int i = 0; i < totalInput; i += 1) {
            data[i] = (s[currentIndex + i] == '1' ? 1 : 0);
        }

        if (!Classes.exists(LableName)) { Classes.Add(LableName); }

        for (int i = 0; i < Classes.GetSize(); i++) {
            if (Classes[i] == LableName) {
                LableId = i;
                break;
            }
        }

    }

    void operator=(const Pattern& copy) {
        if (data != nullptr) { delete[] data; }
        LableId = copy.LableId;
        LableName = copy.LableName;
        height = copy.height; 
        width = copy.width; 
        totalInput = copy.totalInput;
        data = new double[totalInput];
        for (int i = 0; i < totalInput; i++) { data[i] = copy.data[i]; }
    }

    Pattern(const Pattern& copy) { (*this) = copy; }

    void PrintPattern(ostream& out) {
        string s = "#################################\n";
        s = s
            + "## Class Name   : " + LableName + "\n"
            + "## Class ID     : " + to_string(LableId) + "\n"
            + "## Dimensions   : " + to_string(height) + "(h)x" + to_string(width) + "(w)\n"
            + "## Total Inputs : " + to_string(totalInput) + "\n"
            + "#################################\n";

        string horizontalBorder = "";
        for (int i = 0; i < width + 2; i++) { horizontalBorder += "-"; }
        s += horizontalBorder;
        for (int h = 0; h < height; h += 1) {
            string t = "|";
            for (int w = 0; w < width; w += 1) {
                t += (data[h * height + w] == 1 ? "*" : " ");
            }
            s += (t + "|\n");
        }
        s += horizontalBorder;
        out << s;
    }

    double* GetOutputAsDoubleArray(unsigned int outputSize) {
        double* outputArray = new double[outputSize];
        for (unsigned int i = 0; i < LableId; i++) {
            outputArray[i] = 0;
        }
        outputArray[LableId] = 1;
        for (unsigned int i = LableId + 1; i < outputSize; i++) {
            outputArray[i] = 0;
        }
        return outputArray;
    }

    bool GuessedCorrectly(prediction answer) const {
        unsigned int maxIndex = 0;
        double maxValue = answer[0];
        for (unsigned int i = 1; i < answer.getSize(); i++) {
            if (maxValue < answer[i]) {
                maxValue = answer[i];
                maxIndex = i;
            }
        }
        return maxIndex == LableId;
    }

    static string classIdToClassName(unsigned int i) {
        if (i < Classes.GetSize()) {
            return Classes[i];
        }
        return "";
    }

};

Array<string> Pattern::Classes;

Pattern* GetDataArray(const string& dataFile, int& dataCount);
void SplitTestTrain(unsigned int training[], unsigned int testing[], const unsigned int& numData, const unsigned int& checker, unsigned int& trainingCount, unsigned int& testingCount);

#define ENV "WINDOWS"

int main(int argc, char** argv) {

    ifstream inputFile;

    if (argc == 2) {
        inputFile.open(argv[1]);
    }

    if (!inputFile.is_open()) {
        cout << "Unable to open input params file..." << endl;
        return 1;
    }

    string dataFile;
    getline(inputFile, dataFile);

    string temp;
    getline(inputFile, temp);
    unsigned int numThreads = stoi(temp);

    getline(inputFile, temp);
    unsigned int epochs = stoi(temp);

    getline(inputFile, temp);
    unsigned int numHiddenLayers = stoi(temp) + 1;

    getline(inputFile, temp);
    double learningRate = stod(temp);

    getline(inputFile, temp);
    bool useThreads = (temp == "true") ? true : false;

    getline(inputFile, temp);
    string actFunc = temp;

    getline(inputFile, temp);
    bool softMax = (temp == "true") ? true : false;

    unsigned int* hiddenLayers = new unsigned int[numHiddenLayers];
    for (unsigned int i = 0; i < numHiddenLayers - 1; i++) {
        getline(inputFile, temp);
        hiddenLayers[i] = stoi(temp);
    }

    inputFile.close();

    int inputSize = 160*160;
    int outputSize = 20;

    Pattern* dataArray = nullptr;
    cout << "Start: Preparing Input Data" << endl;
    int dataCount = 12500;
    double testingSplit = 0.20;

    unsigned int testingCount = 0;
    unsigned int trainingCount = 0;
    unsigned int checker = testingSplit * 100;
    unsigned int* trainingData = new unsigned int[dataCount];
    unsigned int* testingData = new unsigned int[dataCount];

    cout << trainingCount << " , " << testingCount << endl;
    cout << dataFile << endl;
    SplitTestTrain(trainingData, testingData, dataCount, checker, trainingCount, testingCount);
    dataArray = GetDataArray(dataFile, dataCount);
    if (dataArray == nullptr) {
        cout << "No input data to train on." << endl;
        return 1;
    }
    else {
        cout << "Read " << dataCount << " patterns." << endl;
    }

    for (unsigned int i = 0; i < outputSize; i++) {
        cout << "Class " << i << " : " << Pattern::classIdToClassName(i) << endl;
    }

    bool testEfficiency = true;

    //------------------------------------
    if (dataCount > 0 && dataArray != nullptr) {
        unsigned int inputLayerSize = inputSize;
        unsigned int outputLayerSize = outputSize;
        hiddenLayers[numHiddenLayers - 1] = outputLayerSize;

        if (!testEfficiency) {
            cout << "Start: Creating Neural Network of size: " << numHiddenLayers << ". Structure: { ";
            for (unsigned int i = 0; i < numHiddenLayers; i++) {
                cout << hiddenLayers[i] << " ";
            }
            cout << "} with input size: " << inputLayerSize << ". Will train with " << numThreads << " threads." << endl;
            fcnn fcnn(inputLayerSize, numHiddenLayers, hiddenLayers, numThreads, useThreads, actFunc, softMax);
            cout << "Finish: Creating Neural Network" << endl;


            cout << "Training instances: " << trainingCount << ", Testing instances: " << testingCount << endl;

            // Make the input data
            unsigned int numData = trainingCount;
            double** inputData = new double* [numData];
            double** labels = new double* [numData];

            for (unsigned int i = 0; i < numData; i++) {
                inputData[i] = dataArray[trainingData[i]].data;
                labels[i] = dataArray[trainingData[i]].GetOutputAsDoubleArray(outputSize);

            }
            cout << "Finish: Preparing Input Data" << endl;

            cout << "Start: Training" << endl;
            fcnn.train(inputData, labels, trainingCount, epochs, learningRate, true);
            cout << "Finish: Training" << endl;

            /*************************************************************************/
            /*************************************************************************/
            /*************************************************************************/
            if (trainingCount > 0) {
                cout << "Validating Training:" << endl;
                fcnn.validate(inputData, labels, trainingCount, &cout);
            }
            /*************************************************************************/
            /*************************************************************************/
            /*************************************************************************/
            if (testingCount > 0) {
                cout << "Validating Testing:" << endl;
                double** inputDataTest = new double* [testingCount];
                double** labelsTest = new double* [testingCount];
                for (unsigned int i = 0; i < testingCount; i++) {
                    inputDataTest[i] = dataArray[testingData[i]].data;
                    labelsTest[i] = dataArray[testingData[i]].GetOutputAsDoubleArray(outputSize);
                }
                fcnn.validate(inputDataTest, labelsTest, testingCount, &cout);

                for (unsigned int i = 0; i < testingCount; i++) {
                    delete[] labelsTest[i];
                }
                delete[] inputDataTest;
                delete[] labelsTest;
            }
            /*************************************************************************/
            /*************************************************************************/
            /*************************************************************************/

            for (unsigned int i = 0; i < numData; i++) {
                //delete[] inputData[i];
                delete[] labels[i];
            }
            delete[] inputData;
            delete[] labels;
        }
        else {
            fcnn::determineMostEfficientModel(inputLayerSize, numHiddenLayers, hiddenLayers, outputSize, actFunc, softMax, true);
        }
    }

    // Here
    cout << "Start: Deallocation" << endl;
    delete[] hiddenLayers;
    delete[] trainingData;
    delete[] testingData;
    delete[] dataArray;

    cout << "Finish: Deallocation" << endl;
    return 0;
}


Pattern* GetDataArray(const string& dataFile, int &dataCount) {
    Pattern* data = nullptr;
    ifstream inFile;
    inFile.open(dataFile);
    int lineCounter = 0;

    bool dataCountKnown = (dataCount > 0);
    unsigned int lineSplitter = dataCount / 20;

    if (inFile.is_open()) {
        string s;
        Array<string> dataStrings;
        if (!dataCountKnown) {
            cout << "Reading Data In..." << endl;
            while (getline(inFile, s)) {
                dataStrings.Add(s);
                lineCounter += 1;
                if (lineCounter % 1000 == 0) {
                    cout << lineCounter << endl;
                }
            }
            dataCount = dataStrings.GetSize();
        }
        data = new Pattern[dataCount];
        cout << "Creating Pattern Objects..." << endl;
        lineCounter = 0;
        for (int i = 0; i < dataCount; i++) {
            if (dataCountKnown) { getline(inFile, s); }
            else { s = dataStrings[i]; }
            data[i] = Pattern(s);
            lineCounter += 1;
            if (lineCounter % lineSplitter == 0) {
                cout << lineCounter << endl;
            }
        }
    }
    else {
        cout << "ERROR: Could not open data file" << endl;
    }
    return data;
}


void SplitTestTrain(unsigned int* training, unsigned int* testing, const unsigned int& numData, const unsigned int& checker, unsigned int& trainingCount, unsigned int& testingCount) {
    trainingCount = 0;
    testingCount = 0;
    //cout << training << " , " << testing << endl;
    for (unsigned int i = 0; i < numData; i++) {
        //cout << i;
        if ((i % 100) > checker) {
            //cout << " (training)" << endl;
            training[trainingCount++] = i;
        }
        else {
            //cout << " (testing)" << endl;
            testing[testingCount++] = i;
        }
    }

    cout << trainingCount << " , " << testingCount << endl;

}
