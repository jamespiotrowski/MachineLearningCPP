
#include "MNIST.h"
using namespace std;

#define TRAINING_IMAGES 1
#define TRAINING_LABELS 2
#define TESTING_IMAGES  3
#define TESTING_LABELS  4
#define NUM_THREADS     5
#define NUM_EPOCHS      6
#define NUM_LAYERS      7
#define LEARNING_RATE   8
#define LAYER_START     9

#define ENV             "WINDOWS"
//#define ENV             "LINUX"

int main(int argc, char** argv) {
    ifstream inputFile;
    if (ENV == "WINDOWS") {
        inputFile.open("C:\\Users\\james\\Code\\CPP\\MNIST\\win_InputParameters.txt");
    }
    else {
        inputFile.open("/mnt/c/Users/james/Code/CPP/MNIST/lin_InputParameters.txt");
    }

    if (!inputFile.is_open()) {
        cout << "Unable to open input params file..." << endl;
        return 1;
    }

    string images, labels, imagesTest, labelsTest;
    getline(inputFile, images);
    getline(inputFile, labels);
    getline(inputFile, imagesTest);
    getline(inputFile, labelsTest);

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

    unsigned int * hiddenLayers = new unsigned int[numHiddenLayers];
    for (unsigned int i = 0; i < numHiddenLayers - 1; i++) {
        getline(inputFile, temp);
        hiddenLayers[i] = stoi(temp);
    }
    
    cout << images << endl << labels << endl;

    MNISTDataSet mnistDataSet(images, labels);
    //------------------------------------
    if (mnistDataSet.getSize() > 0) {
        unsigned int inputLayerSize = mnistDataSet[0]->getImageHeight() * mnistDataSet[0]->getImageWidth();
        unsigned int outputLayerSize = mnistDataSet[0]->getLabelSize();
        hiddenLayers[numHiddenLayers - 1] = outputLayerSize;

        cout << "Start: Creating Neural Network of size: " << numHiddenLayers << ". Structure: { ";
        for (unsigned int i = 0; i < numHiddenLayers; i++) {
            cout << hiddenLayers[i] << " ";
        }
        cout << "} with input size: " << inputLayerSize << ". Will train with " << ((useThreads) ? numThreads : 1) << " threads." << endl;
        fcnn mnistNN(inputLayerSize, numHiddenLayers, hiddenLayers, numThreads, useThreads, actFunc);
        cout << "Finish: Creating Neural Network" << endl;

        mnistDataSet[10]->printImage();
        //mnistDataSet[10]->printValues();
        mnistDataSet[21]->printImage();
        mnistDataSet[32]->printImage();
        mnistDataSet[43]->printImage();
        mnistDataSet[54]->printImage();

        /*******************/
        // Make the input data
        unsigned int numData = mnistDataSet.getSize();
        double** inputData = new double* [numData];
        double** labels = new double* [numData];

        cout << "Start: Preparing Input Data" << endl;
        for (unsigned int i = 0; i < numData; i++) {
            inputData[i] = mnistDataSet[i]->convertToInputArray();
            labels[i] = mnistDataSet[i]->getLabelAsDoubleArray();
        }
        cout << "Finish: Preparing Input Data" << endl;

        /* Train? */
        cout << "Start: Training" << endl;
        mnistNN.train(inputData, labels, numData, epochs, learningRate);
        cout << "Finish: Training" << endl;

        /* Validate training */
        cout << "Validating Training:" << endl;
        unsigned int guessedCorrectly = 0;
        for (unsigned int i = 0; i < mnistDataSet.getSize(); i++) {
            double* testNumber = mnistDataSet[i]->convertToInputArray();
            double* label = mnistDataSet[i]->getLabelAsDoubleArray();
            if (mnistDataSet[i]->guessedCorrectly(mnistNN.predict(testNumber))) {
                guessedCorrectly += 1;
            }
            delete[] testNumber;
            delete[] label;
        }
        cout << guessedCorrectly << " out of " << mnistDataSet.getSize() << " correct!" << endl;
        cout << ((double)guessedCorrectly / (double)(mnistDataSet.getSize())) * 100 << "%" << endl;

        /* Validate */
        cout << "Validating Testing:" << endl;
        MNISTDataSet testMnistDataSet(imagesTest, labelsTest);
        if (testMnistDataSet.getSize() > 0) {
            unsigned int guessedCorrectly = 0;
            for (unsigned int i = 0; i < testMnistDataSet.getSize(); i++) {
                double* testNumber = testMnistDataSet[i]->convertToInputArray();
                double* label = testMnistDataSet[i]->getLabelAsDoubleArray();
                if (testMnistDataSet[i]->guessedCorrectly(mnistNN.predict(testNumber))) {
                    guessedCorrectly += 1;
                }
                delete[] testNumber;
                delete[] label;
            }
            cout << guessedCorrectly << " out of " << testMnistDataSet.getSize() << " correct!" << endl;
            cout << ((double)guessedCorrectly / (double)(testMnistDataSet.getSize())) * 100 << "%" << endl;
        }

        /* Clean Up */
        cout << "Start: Deallocation" << endl;
        delete[] hiddenLayers;
        for (unsigned int i = 0; i < numData; i++) {
            delete[] inputData[i];
            delete[] labels[i];
        }
        delete[] inputData;
        delete[] labels;
        cout << "Finish: Deallocation" << endl;
    }
    return 0;
}
