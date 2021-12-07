#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <omp.h>
#include <iomanip>

#include "task.h"

const std::string inputGridFile = "grid.txt";
const std::string magnetizationFolderName = "magnetizationInputs/";
const std::string inputMagnetizationFile = "_magnetization.txt";
const std::string outputFolder = "resultOutputs/";
const std::string outputFile = "_result.txt";
const std::string generateParams = "params.txt";
const std::string logsFolderName = "logs/";
const std::string logFileName = "_log.txt";
const std::string calculateInfoFile = "calculateInfoFile.txt";

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

const int maxUsedValuesAmount = 3;
const int allCountFilter = 2; // позволяет выбирать один из заданного значения
int xsize = 0; // yep, wanna global variable. why not?

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
	int nonZeroValues = 0;
	for (double value : values) if (value != 0) nonZeroValues++;

	logFile << "fileNumber: " << fileNumber << " with values (amount " << nonZeroValues << "):";
	for (int i = 0; i < values.size(); i++) {
		if ((i) % xsize == 0) logFile << std::endl;
		logFile << values[i] << " ";
	}
	logFile << std::endl;
	logFile << "---------------------------------------------" << std::endl;
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


uint64_t calculateAmoutOfGeneratedFilesRecur(int gridSize) {
	uint64_t result = 0;
	int minsize = std::min(gridSize, maxUsedValuesAmount);

	for (int j = minsize; j > 0; j--) {
		uint64_t combinationsFromGridResult = 1;
		// сколькими способами можно вытащить j клеток из gridSize 
		// (неупорядочкенная выборка из gridSize по j)
		int diff = gridSize - j + 1;
		for (int i = diff; i <= gridSize; i++) combinationsFromGridResult *= i;
		for (int i = 2; i <= j; i++) combinationsFromGridResult /= i;

		// перестановки с повторением из J значений множества ReferenceValues:
		uint64_t referenceValuesCombinations = 1;
		for (int i = 0; i < j; i++) referenceValuesCombinations *= referenceValues.size();

		result += combinationsFromGridResult * referenceValuesCombinations;
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


void magicShit(int& fileNumber, int& allcount, int currentValuesLvl, int usedValuesCount, std::vector<float> values) {
	if (usedValuesCount > maxUsedValuesAmount) return;

	// if all Values have been formed - generate input file
	if (currentValuesLvl == values.size()) {
		if (usedValuesCount <= maxUsedValuesAmount) {
			if (allcount++ % allCountFilter == 0) {
				auto t1 = high_resolution_clock::now();

				generateTxtFile(fileNumber, values);
				debugPrintValues(values, xsize, fileNumber);

				auto t2 = high_resolution_clock::now();
				/* Getting number of milliseconds as an integer. */
				auto oneTookMs = duration_cast<milliseconds>(t2 - t1);

				if (fileNumber % 50 == 0) {
					long leftToDo = generatedFilesAmount - fileNumber;
					int remainigSeconds = leftToDo * oneTookMs.count() / 1000;
					float progressPercentage = (float)fileNumber / generatedFilesAmount * 100;

					std::cout << "\r generated " << fileNumber << " file. One Took: " << oneTookMs.count();
					std::cout << "ms. Progerss: " << std::setprecision(4) << progressPercentage << "%, remaining seconds: " << remainigSeconds;
				}

				fileNumber++;
			}
		}
		return;
	}

	// continute forming Values:
	int nextLvl = currentValuesLvl + 1;
	int nextUsedValuesCount = usedValuesCount + 1;

	// forming with 0 in currentValuesLvl place
	magicShit(fileNumber, allcount, nextLvl, usedValuesCount, values);

	for (int i = 0; i < referenceValues.size(); i++) {
		values[currentValuesLvl] = referenceValues[i];
		magicShit(fileNumber, allcount, nextLvl, nextUsedValuesCount, values);
	}
}


void calculateDirectTasks() {
	int currentCount = 0;
	omp_lock_t writelock;
	omp_init_lock(&writelock);

	auto calculatingStartTime = high_resolution_clock::now();

#pragma omp parallel for
	for (int i = 0; i < generatedFilesAmount; i++) {
		int numOfThreads = omp_get_num_threads();
		int fileNumber = i;

#pragma omp atomic
		currentCount++;

		// -------------- main work ---------------------
		caluclateDirectTask(fileNumber);
		// -------------- main work ---------------------

		if (currentCount % 50 == 0)
		{
			omp_set_lock(&writelock);
			auto currentCalculationTime = high_resolution_clock::now();
			auto alreadyTookMs = duration_cast<milliseconds>(currentCalculationTime - calculatingStartTime);

			int tookSeconds = alreadyTookMs.count() / 1000;
			float progressPercentage = (float)currentCount / generatedFilesAmount * 100;

			std::cout << "\r calculated " << currentCount << " files. With " << numOfThreads << " threads. ";
			std::cout << "Progeress: " << std::setprecision(4) << progressPercentage << "%, it took: " << alreadyTookMs.count() / 1000  << " sec";
			omp_unset_lock(&writelock);
		}
	}
}

int main() {
	logFile.open(logsFolderName + currentDateTime() + logFileName);

	auto gridSizes = getGridSize();
	xsize = gridSizes[0];
	int zsize = gridSizes[2];
	auto gridSize = xsize * zsize;
	auto values = std::vector<float>(gridSize);

	int maxSquareSize = std::min(gridSizes[0], gridSizes[2]); // min from X or Z <3
	int fileNumber = 0;
	int allcount = 0;

	generatedFilesAmount = calculateAmoutOfGeneratedFilesRecur(gridSize);
	generatedFilesAmount /= allCountFilter;
	std::cout << "Will be generated and calculated " << generatedFilesAmount << " files" << std::endl;
	logFile << "Will be generated and calculated " << generatedFilesAmount << " files" << std::endl;

	std::cout << "start generating files" << std::endl;
	magicShit(fileNumber, allcount, 0, 0, values);

	std::cout << std::endl << "start calculating direct tasks" << std::endl;
	calculateDirectTasks();

	return 0;
}
