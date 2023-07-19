#pragma once
#include <random>
#include <iostream>
#include <pthread.h>
#include <chrono>
#include <ctime>   

using namespace std;
/*************************
 * Fully Connected Neural Network Class for C++
 * Code Written by James Piotrowski
 * University of Nevada Las Vegas
 * December 2021
*************************/
class prediction {
private:
    unsigned int size;
    double* arr = nullptr;
public:
    unsigned int getSize() const { return size; }
    prediction() {}
    prediction(unsigned int s) {
        size = s;
        arr = new double[size];
    }
    ~prediction() {
        if (arr != nullptr) {
            delete[] arr;
        }
    }
    prediction(const prediction& copy) {
        arr = nullptr;
        size = 0;
        (*this) = copy;
    }
    void operator=(const prediction& copy) {
        if (arr != nullptr) {
            delete[] arr;
            arr = nullptr;
        }
        size = copy.size;
        if (size == 0) {
            return;
        }
        arr = new double[size];
        for (unsigned int i = 0; i < size; i++) {
            arr[i] = copy.arr[i];
        }
    }
    double& operator[](unsigned int ind) {
        if (ind >= size) {
            exit(1);
        }
        return arr[ind];
    }
};



class fcnn {
private:
    // Size of the input
    unsigned int inputSize = 0;
    // Number of layers
    unsigned int numHiddenLayers = 0;
    // Size of the layers
    unsigned int* layerSize = nullptr;
    // Weights for each layer 
    double*** w = nullptr;
    // Node Bias
    double** b = nullptr;
    // Node Equation (s = i[0]*w[0] + i[1]*w[i] ... i[n]*w[n] + b) 
    double** s = nullptr;
    // Activation Node
    double** y = nullptr;
    // Learning Rate
    double e = 0.1;
    // Place to store error calculations
    double* error = nullptr;
    double** errStore = nullptr;
    // Place to store NEW weights
    double*** wN = nullptr;
    double** bN = nullptr;
    // Parallelization
    unsigned int numThreads = 1;
    pthread_t* threads = NULL;
    static const unsigned int maxThreads = 32000;
    bool useThreads = false;

    enum activationFunction {
        SIGMOID
        , TANH
        , RELU
        , LEAKY_RELU
        , SM_SIGMOID
        , SM_TANH
        , SM_RELU
        , SM_LEAKY_RELU
    };

    activationFunction af = SIGMOID;

    enum functionToCall {
        TRAIN_THREAD
    };

    void deleteAll();

    double Init() const {
        return ((rand() % 2) == 0 ? (((double)(rand() % 1000)) / 100000.0) : (-1.0 * (((double)(rand() % 1000)) / 100000.0)));
    }

    double computeActivationFunction(double d);
    double computeActivationFunctionDerivative(double d);

    double SoftMax(int i) {
        double sum = 0.0;
        for (unsigned int j = 0; j < layerSize[numHiddenLayers - 1]; j++) {
            sum += exp(s[numHiddenLayers - 1][j]);
        }
        return exp(s[numHiddenLayers - 1][i]) / sum;
    }

    int delta(int i, int j) { return (i == j ? 1 : 0); }

    unsigned int getPreviousSize(unsigned int i) const {
        if (i == 0) {
            return inputSize;
        }
        return layerSize[i - 1];
    }

    unsigned int getNextSize(unsigned int i) const {
        if (i == numHiddenLayers - 1) {
            exit(1);
        }
        return layerSize[i + 1];
    }

    struct helperNode {
        functionToCall ftc;
        fcnn* objectReference;
        void* arguments;
    };

    struct computeForwardNodeStruct {
        int currentLayer;
        unsigned int minNode;
        unsigned int maxNode;
        double* input = nullptr;
    };

    struct trainingResources {
        prediction *p = nullptr;            
        double sumError = 0;
        double divisor = 0;
    };

    struct threadArguments {
        unsigned int threadId;
        pthread_barrier_t* barrier;
        double** dataInput;
        double** dataOutput; 
        unsigned int numData; 
        unsigned int epochs;
        trainingResources* tr;
    };

    void setThreadArguments(threadArguments &ta, unsigned int threadId, pthread_barrier_t* barrier, double** dataInput, double** dataOutput, unsigned int numData, unsigned int epochs, trainingResources *tr) {
        ta.threadId = threadId;
        ta.barrier = barrier;
        ta.dataInput = dataInput;
        ta.dataOutput = dataOutput;
        ta.numData = numData;
        ta.epochs = epochs;
        ta.tr = tr;
    }

    static void computeThreadRange(unsigned int& minNode, unsigned int& maxNode, const unsigned int& threadNum, const unsigned int& numThreads, const unsigned int& layerSize) {
        
        string msg = "{" + to_string(threadNum) + "," + to_string(numThreads) + "," + to_string(layerSize) + "}\n";
        cout << msg;
        minNode = (unsigned int)((double)threadNum * (((double)layerSize) / ((double)numThreads)));
        maxNode = (unsigned int)((double)(threadNum + 1) * (((double)layerSize) / ((double)numThreads)));

        /* old... and stupid :)
        unsigned int slice = layerSize / numThreads;
        if (!((layerSize % numThreads) != 0)) { slice += 1; }
        minNode = threadNum * slice;
        maxNode = minNode + slice;
        if (maxNode > layerSize) { maxNode = layerSize; }
        */
    }

    computeForwardNodeStruct* threadArgs = nullptr;
    helperNode* threadHelpers = nullptr;

    void trainMaster(double** dataInput, double** dataOuput, unsigned int numData, unsigned int epochs);
    void trainIndividual(void* args);

public:
    fcnn() {}
    fcnn(unsigned int is, unsigned int nhl, unsigned int hls[], unsigned int nt, bool ut, string actFunc);
    ~fcnn();
    fcnn(const fcnn& copy);
    void operator=(const fcnn& copy);
    prediction predict(double* input);
    void predict(double* input, unsigned int** threadRange, prediction& p, pthread_barrier_t* barrier);
    void computeForwardNode(unsigned int layer, unsigned int node, double* inp, unsigned int inpSize);
    void computeBackwardNode(int layer, unsigned int node, double* inp, unsigned int inpSize);


    void train(double** dataInput, double** dataOuput, unsigned int numData, unsigned int epochs, double lr);


    static void* callMemberFunctionForThread(void* args);
    void* computeForwardNodeT(void* args);
    void* computeBackwardNodeT(void* args);

};


fcnn::fcnn(unsigned int is, unsigned int nhl, unsigned int hls[], unsigned int nt, bool ut, string actFunc) {
    numThreads = nt;
    useThreads = ut;
    threadArgs = new computeForwardNodeStruct[nt];
    threadHelpers = new helperNode[nt];
    if (numThreads > maxThreads) {
        numThreads = maxThreads;
    }
    threads = new pthread_t[nt];
    inputSize = is;
    numHiddenLayers = nhl;
    layerSize = new unsigned int[numHiddenLayers];
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        layerSize[i] = hls[i];
    }
    // Init each layer
    s = new double* [numHiddenLayers];
    b = new double* [numHiddenLayers];
    y = new double* [numHiddenLayers];
    errStore = new double* [numHiddenLayers];
    bN = new double* [numHiddenLayers];
    w = new double** [numHiddenLayers];
    wN = new double** [numHiddenLayers];
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        s[i] = new double[layerSize[i]];
        b[i] = new double[layerSize[i]];
        y[i] = new double[layerSize[i]];
        for (unsigned int j = 0; j < layerSize[i]; j++) { y[i][j] = 0.0; s[i][j] = 0.0; b[i][j] = 0.0; }
        errStore[i] = new double[layerSize[i]];
        bN[i] = new double[layerSize[i]];
        unsigned int lowerSize = inputSize;
        if (i > 0) {
            lowerSize = layerSize[i - 1];
        }
        w[i] = new double* [lowerSize];
        wN[i] = new double* [lowerSize];
        for (unsigned int j = 0; j < lowerSize; j++) {
            w[i][j] = new double[layerSize[i]];
            wN[i][j] = new double[layerSize[i]];
            for (unsigned int k = 0; k < layerSize[i]; k++) {
                w[i][j][k] = Init();
                wN[i][j][k] = w[i][j][k];
            }
        }
        for (unsigned int j = 0; j < layerSize[i]; j++) {
            b[i][j] = Init();
        }
    }
    error = new double[layerSize[numHiddenLayers - 1]];

    if (actFunc == "SIGMOID") { af = SIGMOID; }
    else if (actFunc == "TANH") { af = TANH; }
    else if (actFunc == "RELU") { af = RELU; }
    else if (actFunc == "LEAKY_RELU") { af = LEAKY_RELU; }
    else if (actFunc == "SM_SIGMOID") { af = SM_SIGMOID; }
    else if (actFunc == "SM_TANH") { af = SM_TANH; }
    else if (actFunc == "SM_RELU") { af = SM_RELU; }
    else if (actFunc == "SM_LEAKY_RELU") { af = SM_LEAKY_RELU; }
  
}

fcnn::~fcnn() {
    deleteAll();
}

void fcnn::deleteAll() {
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        unsigned int lowerSize = inputSize;
        if (i > 0) {
            lowerSize = layerSize[i - 1];
        }
        for (unsigned int j = 0; j < lowerSize; j++) {
            delete[] w[i][j];
            delete[] wN[i][j];
        }
        delete[] s[i];
        delete[] b[i];
        delete[] y[i];
        delete[] errStore[i];
        delete[] bN[i];
        delete[] w[i];
        delete[] wN[i];
    }
    delete[] layerSize;
    // Init each layer
    delete[] s;
    delete[] b;
    delete[] y;
    delete[] errStore;
    delete[] bN;
    delete[] w;
    delete[] wN;
    delete[] error;

    delete[] threads;
    delete[] threadArgs;
    delete[] threadHelpers;
    threadHelpers = nullptr;
    threadArgs = nullptr;
    threads = nullptr;

    inputSize = 0;
    numHiddenLayers = 0;
    layerSize = nullptr;
    w = nullptr;
    b = nullptr;
    s = nullptr;
    y = nullptr;
    errStore = nullptr;
    wN = nullptr;
    bN = nullptr;
    error = nullptr;
}

fcnn::fcnn(const fcnn& copy) {
    (*this) = copy;
}

void fcnn::operator=(const fcnn& copy) {
    deleteAll();
    af = copy.af;
    useThreads = copy.useThreads;
    inputSize = copy.inputSize;
    numHiddenLayers = copy.numHiddenLayers;
    numThreads = copy.numThreads;
    threads = new pthread_t[numThreads];
    threadArgs = new computeForwardNodeStruct[numThreads];
    threadHelpers = new helperNode[numThreads];
    layerSize = new unsigned int[numHiddenLayers];
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        layerSize[i] = copy.layerSize[i];
    }
    e = copy.e;
    // Init each layer
    s = new double* [numHiddenLayers];
    b = new double* [numHiddenLayers];
    y = new double* [numHiddenLayers];
    errStore = new double* [numHiddenLayers];
    bN = new double* [numHiddenLayers];
    w = new double** [numHiddenLayers];
    wN = new double** [numHiddenLayers];
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        s[i] = new double[layerSize[i]];
        b[i] = new double[layerSize[i]];
        y[i] = new double[layerSize[i]];
        for (unsigned int j = 0; j < layerSize[i]; j++) { y[i][j] = 0.0; s[i][j] = 0.0; b[i][j] = 0.0; }
        errStore[i] = new double[layerSize[i]];
        bN[i] = new double[layerSize[i]];
        unsigned int lowerSize = inputSize;
        if (i > 0) {
            lowerSize = layerSize[i - 1];
        }
        w[i] = new double* [lowerSize];
        wN[i] = new double* [lowerSize];
        for (unsigned int j = 0; j < lowerSize; j++) {
            w[i][j] = new double[layerSize[i]];
            wN[i][j] = new double[layerSize[i]];
            for (unsigned int k = 0; k < layerSize[i]; k++) {
                w[i][j][k] = copy.w[i][j][k];
            }
        }
        for (unsigned int j = 0; j < layerSize[i]; j++) {
            b[i][j] = copy.b[i][j];
        }
    }
    error = new double[layerSize[numHiddenLayers - 1]];
}


/******************************************************************************
 * PREDICTION FUNCTION
-------------------------------------------------------------------------------
 * This member function takes input for the neural network
   and makes a prediction.
*******************************************************************************/
prediction fcnn::predict(double* input) {
    prediction p(layerSize[numHiddenLayers - 1]);
    bool lastLoop = false;
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        if (i == (numHiddenLayers - 1)) {
            lastLoop = true;
        }
        double* inp = (i > 0) ? y[i - 1] : input;
        unsigned int inpSize = (i > 0) ? layerSize[i - 1] : inputSize;
        for (unsigned int j = 0; j < layerSize[i]; j++) {
            computeForwardNode(i, j, inp, inpSize);
            if (lastLoop) {
                p[j] = y[i][j];
            }
        }
    }
    return p;
}

/******************************************************************************
 * TRAIN FUNCTION
-------------------------------------------------------------------------------
 *
*******************************************************************************/
void fcnn::train(double** dataInput, double** dataOutput, unsigned int numData, unsigned int epochs = 50, double lr = 0.1) {
    e = lr; 
    if (useThreads) {
        trainMaster(dataInput, dataOutput, numData, epochs);
    }
    else {
        std::chrono::duration<double> elapsed_seconds;
        unsigned int totalTime = 0;
        for (unsigned int epc = 0; epc < epochs; epc++) {
            auto start = std::chrono::system_clock::now();
            cout << "Start: Epoch " << epc << endl;
            double sumError = 0;
            double divisor = 0;
            for (unsigned int d = 0; d < numData; d++) {
                if (d % (numData / 100) == 0) {
                    cout << ((double)d / (double)numData) * 100.0 << "%" << endl;
                }
                // Make a prediction
                prediction p = predict(dataInput[d]);

                for (unsigned int i = 0; i < layerSize[numHiddenLayers - 1]; i++) {
                    error[i] = p[i] - dataOutput[d][i];
                    sumError += (error[i] * error[i]);
                    divisor += 1;
                    errStore[numHiddenLayers - 1][i] = (e * error[i] * computeActivationFunctionDerivative(y[numHiddenLayers - 1][i]));
                    bN[numHiddenLayers - 1][i] = b[numHiddenLayers - 1][i] - errStore[numHiddenLayers - 1][i];
                }

                if (numHiddenLayers > 1) {
                    for (int h = numHiddenLayers - 2; h >= 0; h--) {
                        for (unsigned int j = 0; j < layerSize[h]; j++) {
                            computeBackwardNode(h, j, y[h], layerSize[h + 1]);
                        }
                    }
                }

                for (unsigned int j = 0; j < inputSize; j++) {
                    computeBackwardNode(-1, j, dataInput[d], layerSize[0]);
                }
               
                // Reassignment
                for (unsigned int i = 0; i < numHiddenLayers; i++) {
                    unsigned int lowerSize = inputSize;
                    if (i > 0) {
                        lowerSize = layerSize[i - 1];
                    }
                    for (unsigned int k = 0; k < layerSize[i]; k++) {
                        for (unsigned int j = 0; j < lowerSize; j++) {
                            w[i][j][k] = wN[i][j][k];
                        }
                        b[i][k] = bN[i][k];
                    }
                }
            }

            auto end = std::chrono::system_clock::now();
            elapsed_seconds = end - start;
            totalTime += (unsigned int)elapsed_seconds.count();
            std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
            cout << "Finish: Epoch " << epc << endl;
            std::cout << epc << " (training): " << sqrt(sumError) / (double)numData << " RMSE per sample" << std::endl;
            std::cout << epc << " (training): " << sqrt(sumError) / divisor << " RMSE per output" << std::endl;
        }

        cout << "Average Epoch Time: " << (double)totalTime / (double)epochs << endl;
    }
}

void fcnn::computeForwardNode(unsigned int layer, unsigned int node, double* input, unsigned int inputSize) {
    s[layer][node] = 0;
    for (unsigned int j = 0; j < inputSize; j++) {
        s[layer][node] += w[layer][j][node] * input[j];
    }
    s[layer][node] += b[layer][node];
    y[layer][node] = computeActivationFunction(s[layer][node]);
}

/******************************************************************************
 * computeBackwardNode FUNCTION
-------------------------------------------------------------------------------
 *
*******************************************************************************/
void fcnn::computeBackwardNode(int layer, unsigned int node, double* input, unsigned int inputSize) {
    if (layer >= 0) {
        double sumZ = 0;
        for (unsigned int i = 0; i < inputSize; i++) {
            wN[layer + 1][node][i] = w[layer + 1][node][i] - ((errStore[layer + 1][i]) * input[node]);
            sumZ += errStore[layer + 1][i] * w[layer + 1][node][i];
        }
        errStore[layer][node] = sumZ * computeActivationFunctionDerivative(input[node]); //(input[node] * (1.0 - input[node])); // Sigmoid
        bN[layer][node] = b[layer][node] - errStore[layer][node];
    }
    else {
        for (unsigned int i = 0; i < inputSize; i++) {
            wN[0][node][i] = w[0][node][i] - ((errStore[0][i]) * input[node]);
        }
    }
}


/******************************************************************************
 * callMemberFunctionForThread FUNCTION
-------------------------------------------------------------------------------
 * Calling member functions using threads is difficult,
    so we use a static function and an object reference instead.

    This function will be used to call multiple member functions,
    thus a switch statement is used to determine which function to call.
*******************************************************************************/
void* fcnn::callMemberFunctionForThread(void* args) {
    helperNode* hn = (helperNode*)args;
    switch (hn->ftc) {
    case(TRAIN_THREAD): {
        hn->objectReference->trainIndividual(hn->arguments);
        break;
    }
    default: {
        break;
    }
    }
    pthread_exit(NULL);
    return nullptr;
}

double fcnn::computeActivationFunction(double d) {
    switch(af) {
    case(SIGMOID) : {
        return (1.0 / (1.0 + exp(-d)));
        break;
    }
    case(RELU): { 
        return (d > 0.0) ? d : 0.0; // Is this right?
        break;
    }
    case(TANH): {
        return tanh(d); 
        break;
    }
    default: {
        cout << "ERROR: ACTIVATION FUNCTION NOT SELECTED." << endl;
        return 0;
    }
    }
}

double fcnn::computeActivationFunctionDerivative(double d) {
    switch (af) {
    case(SIGMOID): {
        return (d * (1.0 - d));
        break;
    }
    case(RELU): {
        return (d > 0.0) ? 1.0 : 0.0; // Is this right?
        break;
    }
    case(TANH): {
        return 1.0 - (d * d); 
        break;
    }
    default: {
        cout << "ERROR: ACTIVATION FUNCTION NOT SELECTED." << endl;
        return 0;
    }
    }
}

/******************************************************************************
############################################################################### 
############################################################################### 
############################################################################### 
############################################################################### 
*******************************************************************************/


void fcnn::trainMaster(double** dataInput, double** dataOuput, unsigned int numData, unsigned int epochs) {
    /* If threads wider than output, we should update. */
    if (numThreads > layerSize[numHiddenLayers - 1]) {
        numThreads = layerSize[numHiddenLayers - 1];
    }
    /* Allocate threads */
    pthread_t* threads = new pthread_t[numThreads];
    /* Variable for output */
    unsigned int threadResult = 0;
    /* Barrier used to keep threads from moving on without one another */
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, numThreads);
    /* Thread Arguments */
    threadArguments* ta = new threadArguments[numThreads];
    helperNode* th = new helperNode[numThreads];
    trainingResources tr;
    tr.p = new prediction(layerSize[numHiddenLayers - 1]);
    cout << "Setting up threads for training...." << endl;
    /* Load up threads and send them off */
    for (unsigned int i = 0; i < numThreads; i++) {
        setThreadArguments(ta[i], i, &barrier, dataInput, dataOuput, numData, epochs, &tr);
        th[i].ftc = TRAIN_THREAD;
        th[i].objectReference = this;
        th[i].arguments = (void*)(ta + i);
        threadResult = pthread_create((threads + i), NULL, callMemberFunctionForThread, (void*)(th + i));
    }
    cout << "Threads are training...." << endl;

    for (unsigned int i = 0; i < numThreads; i++) {
        threadResult = pthread_join(*(threads + i), NULL);
    }


    delete[] ta;
    delete[] threads;
}

void fcnn::predict(double* input, unsigned int **threadRange, prediction &p, pthread_barrier_t* barrier) {
    bool lastLoop = false;
    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        if (i == (numHiddenLayers - 1)) {
            lastLoop = true;
        }
        double* inp = (i > 0) ? y[i - 1] : input;
        unsigned int inpSize = (i > 0) ? layerSize[i - 1] : inputSize;
        for (unsigned int j = threadRange[i][0]; j < threadRange[i][1]; j++) {
            computeForwardNode(i, j, inp, inpSize);
            if (lastLoop) {
                p[j] = y[i][j];
            }
        }
        if (barrier != nullptr) {
            pthread_barrier_wait(barrier);
        }
    }
}

void fcnn::trainIndividual(void* args) {
    threadArguments* ta = (threadArguments*)args;
    unsigned int threadId = ta->threadId;
    double** dataInput = ta->dataInput;
    double** dataOutput = ta->dataOutput;
    unsigned int numData = ta->numData;
    unsigned int epochs = ta->epochs;
    pthread_barrier_t *barrier = ta->barrier;
    trainingResources* tr = ta->tr;
    unsigned int** threadRange = new unsigned int*[numHiddenLayers + 1];
    threadRange[numHiddenLayers] = new unsigned int[2]; // Start (1) & End (2)

    computeThreadRange(threadRange[numHiddenLayers][0], threadRange[numHiddenLayers][1], threadId, numThreads, inputSize);

    string msg = "[Thread " + to_string(threadId) + "]: Reporting In! I am responsible for the following nodes: [ {" 
        + to_string(threadRange[numHiddenLayers][0]) + "," + to_string(threadRange[numHiddenLayers][1]) + "}";

    for (unsigned int i = 0; i < numHiddenLayers; i++) {
        threadRange[i] = new unsigned int[2]; // Start (1) & End (2)
        computeThreadRange(threadRange[i][0], threadRange[i][1], threadId, numThreads, layerSize[i]);
        msg += " {" + to_string(threadRange[i][0]) + "," + to_string(threadRange[i][1]) + "}";
    }
    msg += "]\n";
    cout << msg;

    for (unsigned int epc = 0; epc < epochs; epc++) {
        for (unsigned int d = 0; d < numData; d++) {
            //cout << d << endl;
            if (d % (numData / 10) == 0) {
                msg = "[Thread " + to_string(threadId) + "]: " + to_string(((double)d / (double)numData) * 100.0) + "%\n";
                cout << msg;
            }

            /* Forward Computation (Prediction) */
            predict(dataInput[d], threadRange, *(tr->p), barrier);

            /* Set up back prop using the prediction */
            for (unsigned int i = threadRange[numHiddenLayers - 1][0]; i < threadRange[numHiddenLayers - 1][1]; i++) {
                error[i] = tr->p->operator[](i) - dataOutput[d][i];
                errStore[numHiddenLayers - 1][i] = (e * error[i] * computeActivationFunctionDerivative(y[numHiddenLayers - 1][i]));
                bN[numHiddenLayers - 1][i] = b[numHiddenLayers - 1][i] - errStore[numHiddenLayers - 1][i];
            }
            pthread_barrier_wait(barrier);

            /* Back prop for all hidden layers */
            if (numHiddenLayers > 1) {
                for (int h = numHiddenLayers - 2; h >= 0; h--) {
                    for (unsigned int j = threadRange[h][0]; j < threadRange[h][1]; j++) {
                        computeBackwardNode(h, j, y[h], layerSize[h + 1]);
                    }
                    pthread_barrier_wait(barrier);
                }
            }

            unsigned int size = layerSize[0];
            /* For input layer */
            for (unsigned int j = threadRange[numHiddenLayers][0]; j < threadRange[numHiddenLayers][1]; j++) {
                for (unsigned int i = 0; i < size; i++) {
                    wN[0][j][i] = w[0][j][i] - ((errStore[0][i]) * dataInput[d][j]);
                }
            }
            pthread_barrier_wait(barrier);

 

            // Reassignment
            for (unsigned int i = 0; i < numHiddenLayers; i++) {
                /* Determine the size of the previous layer */
                unsigned int lowerSize = inputSize;
                if (i > 0) {
                    lowerSize = layerSize[i - 1];
                }
                /* Reassign weights */
                for (unsigned int k = threadRange[i][0]; k < threadRange[i][1]; k++) {
                    for (unsigned int j = 0; j < lowerSize; j++) { 
                        w[i][j][k] = wN[i][j][k];
                    }
                    b[i][k] = bN[i][k];
                }
            } 
            pthread_barrier_wait(barrier);
        }
        msg = "[Thread " + to_string(threadId) + "]: Finished epoch " + to_string(epc) + "\n";
        cout << msg;
    }

    msg = "[Thread " + to_string(threadId) + "]: On the way out!\n";
    cout << msg;
}

