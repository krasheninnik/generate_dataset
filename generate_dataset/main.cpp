#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <ctime>
#include "task.h"

const std::string inputGridFile = "grid.txt";
const std::string magnetizationFolderName = "magnetizationInputs/";
const std::string inputMagnetizationFile = "_magnetization.txt";
const std::string outputFolder = "resultOutputs/";
const std::string outputFile = "_result.txt";
const std::string generateParams = "params.txt";
const std::string logsFolderName = "logs/";
const std::string logFileName = "_log.txt";

// grid 5 x 10
// values: 0, 0.5, 1

const std::vector<float> referenceValues {0.5, 1};
uint64_t generatedFilesAmount = 1;
std::ofstream logFile;

std::vector<int> getGridSize() {
	std::ifstream fin(inputGridFile);
	int temp = 0; 
	double temp1 = 0;
	int xsize = 0, ysize = 0, zsize = 0;
	fin >> temp >> temp >> temp >> temp;
	fin >> temp1 >> temp1 >> temp1;
	fin >> temp >> temp >> xsize;
	fin >> temp >> temp >> ysize;
	fin >> temp >> temp >> zsize;
	return std::vector<int>{ xsize, ysize, zsize };
}

std::string getMagnetizationFileName(int number) { return magnetizationFolderName + std::to_string(number) + inputMagnetizationFile; }
std::string getOutputFileName(int number) { return outputFolder + std::to_string(number) + outputFile; }


void generateTxtFile(int number, std::vector<float> &values) {
	std::ofstream fout(getMagnetizationFileName(number));
	std::for_each(values.begin(), values.end(), [&](double x) {fout << "0 0 " + std::to_string(x) << std::endl; });
	fout.close();
}

void caluclateDirectTask(int number)
{
	Task().launch(inputGridFile, getMagnetizationFileName(number), getOutputFileName(number));
}

void debugPrintValues(std::vector<float> values, int xsize, int fileNumber) {
	logFile << "fileNumber: " << fileNumber << " with values:";
	for (int i = 0; i < values.size(); i++) {
		if ((i) % xsize == 0) logFile << std::endl;
		logFile << values[i] << " ";
	}
	logFile << std::endl;
}

void setValueInSquare(int squareSize, int zpos, int xsize, int xpos, float currentValueInSquare, std::vector<float>& values) {
	for (int zi = 0; zi < squareSize; zi++) {
		int zOffset = (zpos + zi) * xsize;
		for (int xi = 0; xi < squareSize; xi++) {
			values[zOffset + xpos + xi] = currentValueInSquare;
		}
	}
}

void makeMagicThings(int xpos, int zpos, int squareSize, int &fileNumber, std::vector<float>& values, int xsize) {
	logFile << "squareSize: " << squareSize << " xpos,zpos: " << xpos << " " << zpos << " fileNumber: " << fileNumber << std::endl;
	for (float currentValueInSquare : referenceValues) {
		// set values
		setValueInSquare(squareSize, zpos, xsize, xpos, currentValueInSquare, values); // calculate values

		generateTxtFile(fileNumber, values);
		caluclateDirectTask(fileNumber);
		debugPrintValues(values, xsize, fileNumber);

		// clear values
		setValueInSquare(squareSize, zpos, xsize, xpos, 0, values); // clear values
		fileNumber++;
	}
}

int calculateAmoutOfGeneratedFiles(int xsize, int zsize, int maxSquareSize) {
	int result = 1; // with all 0-s

	for (int i = 0; i < maxSquareSize; i++) {
		result += (xsize - i) * (zsize - i) * referenceValues.size();
	}

	return result;
}

const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d_%H_%M_%S", &tstruct);

	return buf;
}

int main() {
	logFile.open(logsFolderName + currentDateTime() + logFileName);

	auto gridSizes = getGridSize();
	int xsize = gridSizes[0];
	int zsize = gridSizes[2];
	auto gridSize = xsize * zsize;
	auto values = std::vector<float>(gridSize);

	int maxSquareSize = std::min(gridSizes[0], gridSizes[2]); // min from X or Z <3
	int fileNumber = 0;

	int amountOfGeneratedFiles = calculateAmoutOfGeneratedFiles(xsize, zsize, maxSquareSize);
	std::cout << "Will be generated and calculated " << amountOfGeneratedFiles << " files" << std::endl;
	logFile << "Will be generated and calculated " << amountOfGeneratedFiles << " files" << std::endl;

	// calculate with zero
	generateTxtFile(fileNumber, values);
	caluclateDirectTask(fileNumber);
	debugPrintValues(values, xsize, fileNumber);
	fileNumber++;

	for (int squareSize = 1; squareSize <= maxSquareSize; squareSize++)
	{
		int zposition = 0;
		while (zsize - zposition - squareSize >= 0) {
			int xposition = 0;
			while (xsize - xposition - squareSize >= 0) {
				makeMagicThings(xposition, zposition, squareSize, fileNumber, values, xsize);
				xposition++;
				std::cout << "\r progress: " << fileNumber / (float)amountOfGeneratedFiles * 100 << "%";
			}
			zposition++;	
		}
		logFile << "---------------------" << std::endl;
	}

	return 0;
}
