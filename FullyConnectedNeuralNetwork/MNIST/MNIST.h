#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include "FCNN.h"
using namespace std;

/*****************************************
 * MNIST Letter
*****************************************/
class MNISTImage {
private:
    unsigned int imageWidth = 0, imageHeight = 0, labelSize = 0;
    unsigned char** image = nullptr;
    char* label = nullptr;
    //double** image_n = nullptr;
    //const double maxPixelValue = 255.0;
public:
    MNISTImage(unsigned int imageHeight, unsigned int imageWidth, unsigned char** image, unsigned int labelSize, char* label);
    ~MNISTImage();
    MNISTImage(const MNISTImage& copy);
    void operator=(const MNISTImage& copy);
    unsigned char* operator[](unsigned int height);
    char* getLabel();
    double* getLabelAsDoubleArray() const;
    unsigned int getLabelAsInt() const;
    char getLabelAsChar() const;
    void printImage() const;
    void printValues() const;
    double* convertToInputArray() const;
    unsigned int getImageWidth() const { return imageWidth; }
    unsigned int getImageHeight() const { return imageHeight; }
    unsigned int getLabelSize() const { return labelSize; }
    bool guessedCorrectly(prediction p);
};

MNISTImage::MNISTImage(unsigned int imageHeight, unsigned int imageWidth, unsigned char** image
    , unsigned int labelSize, char* label)
    : imageHeight(imageHeight), imageWidth(imageWidth), labelSize(labelSize) {

    this->image = new unsigned char* [imageHeight];
    //this->image_n = new double* [imageHeight];
    for (unsigned int i = 0; i < imageHeight; i += 1) {
        this->image[i] = new unsigned char[imageWidth];
        //this->image_n[i] = new double[imageWidth];
        for (unsigned int j = 0; j < imageWidth; j += 1) {
            this->image[i][j] = image[i][j];
            //this->image_n[i][j] = ((double)image[i][j]) / maxPixelValue;
        }
    }

    this->label = new char[labelSize];
    for (unsigned int i = 0; i < labelSize; i += 1) {
        this->label[i] = label[i];
    }
}

MNISTImage::~MNISTImage() {
    if (image != nullptr) {
        for (unsigned int i = 0; i < imageHeight; i += 1) {
            if (image[i] != nullptr) {
                delete[] image[i];
                //delete[] image_n[i];
            }
        }
        //delete[] image_n;
        delete[] image;
    }

    if (label != nullptr) {
        delete[] label;
    }
}

unsigned char* MNISTImage::operator[](unsigned int i) {
    if (i < imageHeight) {
        return image[i];
    }
    cout << "ERROR <char *MNISTImage::operator[](unsigned int i)> : First index out of bounds." << endl;
    return nullptr;
}

char* MNISTImage::getLabel() {
    return label;
}

unsigned int MNISTImage::getLabelAsInt() const {
    for (unsigned int i = 0; i < labelSize; i += 1) {
        if (label[i]) {
            return i;
        }
    }
    return -1;
}

char MNISTImage::getLabelAsChar() const {
    for (unsigned int i = 0; i < labelSize; i += 1) {
        if (label[i]) {
            return i + 48;
        }
    }
    return '\0';
}

void MNISTImage::printValues() const {
    cout << "/*-----------------------*\\" << endl
        << "| Digit: " << getLabelAsChar() << endl
        << "|*-----------------------*|" << endl
        << endl;

    for (unsigned int r = 0; r < imageHeight; r += 1) {
        for (unsigned int c = 0; c < imageWidth; c += 1) {
            cout << (int)image[r][c] << " , ";
        }
        cout << endl;
    }
    cout << endl
        << "\\*-----------------------*/" << endl << endl;
}


void MNISTImage::printImage() const {
    cout << "/*-----------------------*\\" << endl
        << "| Digit: " << getLabelAsChar() << endl
        << "|*-----------------------*|" << endl
        << endl;

    for (unsigned int r = 0; r < imageHeight; r += 1) {
        for (unsigned int c = 0; c < imageWidth; c += 1) {
            cout << ((image[r][c] > 0) ? '*' : ' ');
        }
        cout << endl;
    }
    cout << endl
        << "\\*-----------------------*/" << endl << endl;
}

MNISTImage::MNISTImage(const MNISTImage& copy) {
    image = nullptr;
    label = nullptr;
    (*this) = copy;
}

void MNISTImage::operator=(const MNISTImage& copy) {
    if (image != nullptr) {
        for (unsigned int i = 0; i < imageHeight; i += 1) {
            if (image[i] != nullptr) {
                delete[] image[i];
            }
        }
        delete[] image;
    }

    if (label != nullptr) {
        delete[] label;
    }

    imageHeight = copy.imageHeight;
    imageWidth = copy.imageWidth;
    labelSize = copy.labelSize;

    this->image = new unsigned char* [imageHeight];
    for (unsigned int i = 0; i < imageHeight; i += 1) {
        this->image[i] = new unsigned char[imageWidth];
        for (unsigned int j = 0; j < imageWidth; j += 1) {
            this->image[i][j] = copy.image[i][j];
        }
    }

    this->label = new char[labelSize];
    for (unsigned int i = 0; i < labelSize; i += 1) {
        this->label[i] = copy.label[i];
    }
}

double* MNISTImage::convertToInputArray() const {
    double normalizedValue = 255.0;
    double* returnArray = new double[imageHeight * imageWidth];
    for (unsigned int i = 0; i < imageHeight; i++) {
        for (unsigned int j = 0; j < imageWidth; j++) {
            returnArray[(imageHeight * i) + j] = ((double)(image[i][j])) / normalizedValue;
        }
    }
    return returnArray;
}

double* MNISTImage::getLabelAsDoubleArray() const {
    double* returnArray = new double[labelSize];
    for (unsigned int i = 0; i < labelSize; i++) {
        returnArray[i] = (double)label[i];
    }
    return returnArray;
}

bool MNISTImage::guessedCorrectly(prediction p) {
    if (p.getSize() > 0) {
        double val = p[0];
        unsigned int bestVal = 0;
        for (unsigned int i = 1; i < p.getSize(); i++) {
            if (p[i] > val) {
                val = p[i];
                bestVal = i;
            }
        }
        return (bestVal == getLabelAsInt());
    }
    return false;
}

/*****************************************
 * MNIST DataSet
*****************************************/
class MNISTDataSet {
private:
    struct ImageHeader {
        uint32_t magicNumber;
        uint32_t maxImages;
        uint32_t Width;
        uint32_t Height;
    };

    struct ImageLabel {
        uint32_t magicNumber;
        uint32_t maxLabels;
    };

    ImageHeader imageHeader;
    ImageLabel imageLabel;
    static int flipBytes(int n);

    MNISTImage** MNISTimages = nullptr;

    static const unsigned int labelSize = 10;
public:
    MNISTDataSet(string images, string labels);
    ~MNISTDataSet();
    MNISTImage* operator[](unsigned int i);
    MNISTImage at(unsigned int i);
    uint32_t getSize() const { return imageHeader.maxImages; }
};

int MNISTDataSet::flipBytes(int n) {
    int b0, b1, b2, b3;
    b0 = (n & 0x000000ff) << 24u;
    b1 = (n & 0x0000ff00) << 8u;
    b2 = (n & 0x00ff0000) >> 8u;
    b3 = (n & 0xff000000) >> 24u;
    return (b0 | b1 | b2 | b3);
}

// "train-images-idx3-ubyte"
// "train-labels-idx1-ubyte"
MNISTDataSet::MNISTDataSet(string images, string labels) {
    // Open data files
    ifstream Imagefile(images, std::ios::binary);
    ifstream Labelfile(labels);
    // If files are open
    if (Imagefile.is_open() && Labelfile.is_open()) {
        // Read in header for images
        imageHeader.magicNumber = 0;
        imageHeader.maxImages = 0;
        imageHeader.Height = 0;
        imageHeader.Width = 0;
        Imagefile.read((char*)&imageHeader.magicNumber, sizeof(imageHeader.magicNumber));
        imageHeader.magicNumber = flipBytes(imageHeader.magicNumber);
        Imagefile.read((char*)&imageHeader.maxImages, sizeof(imageHeader.maxImages));
        imageHeader.maxImages = flipBytes(imageHeader.maxImages);
        Imagefile.read((char*)&imageHeader.Height, sizeof(imageHeader.Height));
        imageHeader.Height = flipBytes(imageHeader.Height);
        Imagefile.read((char*)&imageHeader.Width, sizeof(imageHeader.Width));
        imageHeader.Width = flipBytes(imageHeader.Width);
        // Read in header for labels
        Labelfile.read((char*)&imageLabel.magicNumber, sizeof(imageHeader.magicNumber));
        imageLabel.magicNumber = flipBytes(imageLabel.magicNumber);
        Labelfile.read((char*)&imageLabel.maxLabels, sizeof(imageLabel.maxLabels));
        imageLabel.maxLabels = flipBytes(imageLabel.maxLabels);
        // Show Details 
        cout << "Max Images: " << imageHeader.maxImages << endl;
        cout << "Height: " << imageHeader.Height << endl;
        cout << "Width: " << imageHeader.Width << endl;
        cout << "Magic Number: " << imageHeader.magicNumber << endl;

        // Allocate images 
        MNISTimages = new MNISTImage * [imageHeader.maxImages];
        // Temp data for Images
        char* label = new char[labelSize];
        for (unsigned int i = 0; i < labelSize; i += 1) {
            label[i] = 0;
        }
        unsigned char** image = new unsigned char* [imageHeader.Height];
        for (unsigned int i = 0; i < imageHeader.Height; i++) {
            image[i] = new unsigned char[imageHeader.Width];
        }
        char t_label;
        // For all images
        for (unsigned int p = 0; p < (unsigned int)imageHeader.maxImages; p++) {
            if ((p % ((unsigned int)imageHeader.maxImages / 100)) == 0) {
                cout << p << "....." << endl;
            }
            Labelfile.read(&t_label, 1);
            label[t_label] = 1; //Set truth vector
            //***Read the picture*******
            for (unsigned int r = 0; r < imageHeader.Height; r += 1) {
                for (unsigned int c = 0; c < imageHeader.Width; c += 1) {
                    Imagefile.read((char*)(image[r] + c), sizeof(char));
                }
            }
            MNISTimages[p] = new MNISTImage(imageHeader.Height, imageHeader.Width, image, labelSize, label);
            label[t_label] = 0;
        }
        Imagefile.close();
        Labelfile.close();
        delete[] label;
        for (unsigned int i = 0; i < imageHeader.Height; i++) {
            delete[] image[i];
        }
        delete[] image;
        cout << "Read " << imageHeader.maxImages << " total images!" << endl;
    } 
    else {
        cout << "Failed to open MNIST files!" << endl;
    }
}

MNISTDataSet::~MNISTDataSet() {
    for (unsigned int i = 0; i < imageHeader.maxImages; i += 1) {
        if (MNISTimages[i] != nullptr) {
            delete MNISTimages[i];
        }
    }
    delete[] MNISTimages;
}

MNISTImage* MNISTDataSet::operator[](unsigned int i) {
    if (i < imageHeader.maxImages) {
        return MNISTimages[i];
    }
    cout << "ERROR <MNISTImage *MNISTDataSet::operator=(unsigned int i)> : Index out of bounds." << endl;
    return nullptr;
}

MNISTImage MNISTDataSet::at(unsigned int i) {
    if (i < imageHeader.maxImages) {
        return (*MNISTimages[i]);
    }
    cout << "ERROR <MNISTImage MNISTDataSet::at(unsigned int i)> : Index out of bounds." << endl;
    return MNISTImage(0, 0, nullptr, 0, nullptr);
}




/*************************************************





*************************************************/
