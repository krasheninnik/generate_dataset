#include "task.h"
#define TASK_DIM 3
#define ELEM_DIM 8


void Task::launch(std::string inputAreaFile, std::string inputParametersFile, std::string outputFile) {
	std::ifstream fin;
	fin.open(inputAreaFile);
	double hxMeasure, hyMeasure, x0Grid, x1Grid, y0Grid, y1Grid, z0Grid, z1Grid;
	int nxMeasure, nyMeasure, xStepsGrid, yStepsGrid, zStepsGrid;
	Point p0Measure;

	fin >> hxMeasure >> nxMeasure >> hyMeasure >> nyMeasure;
	fin >> p0Measure.x >> p0Measure.y >> p0Measure.z;

	fin >> x0Grid >> x1Grid >> xStepsGrid;
	fin >> y0Grid >> y1Grid >> yStepsGrid;
	fin >> z0Grid >> z1Grid >> zStepsGrid;

	init(hxMeasure, nxMeasure, hyMeasure, nyMeasure, p0Measure,
		x0Grid, x1Grid, xStepsGrid,
		y0Grid, y1Grid, yStepsGrid,
		z0Grid, z1Grid, zStepsGrid);

	fin.close();

	fin.open(inputParametersFile);
	for (int i = 0; i < elems.size(); i++) {
		fin >> elems[i].p.x >> elems[i].p.y >> elems[i].p.z;
	}
	fin.close();

	getB();

	std::ofstream fout;
	fout.open(outputFile);
	for (int i = 0; i < measures.size(); i++) {
		fout << measures[i].point.x << "\t" << measures[i].point.y << "\t" << measures[i].point.z << "\t";
		fout << measures[i].B.x << "\t" << measures[i].B.y << "\t" << measures[i].B.z << "\t";
		fout << std::endl;
	}
	fout.close();
}

void Task::init(double hxMeasure, int nxMeasure, double hyMeasure, int nyMeasure, Point p0Measure, 
	double x0Grid, double x1Grid, int xStepsGrid,
	double y0Grid, double y1Grid, int yStepsGrid,
	double z0Grid, double z1Grid, int zStepsGrid) {

	fillAxisGrid(xAxisGrid, x0Grid, x1Grid, xStepsGrid);
	fillAxisGrid(yAxisGrid, y0Grid, y1Grid, yStepsGrid);
	fillAxisGrid(zAxisGrid, z0Grid, z1Grid, zStepsGrid);

	int pointsInX = xAxisGrid.size();
	int pointsInY = yAxisGrid.size();
	int pointsInZ = zAxisGrid.size();
	int pointsInXY = pointsInX * pointsInY;

	nodes.resize(pointsInX * pointsInY * pointsInZ);
	int ni = 0;
	for (int zi = 0; zi < pointsInZ; zi++)
	{
		for (int yi = 0; yi < pointsInY; yi++)
		{
			for (int xi = 0; xi < pointsInX; xi++)
			{
				nodes[ni++] = Point(xAxisGrid[xi], yAxisGrid[yi], zAxisGrid[zi]);
			}
		}
	}

	elems.resize((pointsInX - 1) * (pointsInY - 1) * (pointsInZ - 1));

	int ei = 0;
	for (int zi = 0; zi < pointsInZ - 1; zi++)
	{
		for (int yi = 0; yi < pointsInY - 1; yi++)
		{
			for (int xi = 0; xi < pointsInX - 1; xi++)
			{
				int n0 = zi * pointsInXY + yi * pointsInX + xi;
				int n1 = zi * pointsInXY + yi * pointsInX + xi + 1;
				int n2 = zi * pointsInXY + (yi + 1) * pointsInX + xi;
				int n3 = zi * pointsInXY + (yi + 1) * pointsInX + xi + 1;
				int n4 = (zi + 1) * pointsInXY + yi * pointsInX + xi;
				int n5 = (zi + 1) * pointsInXY + yi * pointsInX + xi + 1;
				int n6 = (zi + 1) * pointsInXY + (yi + 1) * pointsInX + xi;
				int n7 = (zi + 1) * pointsInXY + (yi + 1) * pointsInX + xi + 1;
				int value = (xi + 1) * (yi + 1) * (zi + 1);

				Point p = Point(0, 0, 0);
				elems[ei++] = FiniteElem(std::vector<int>{n0, n1, n2, n3, n4, n5, n6, n7}, p);
			}
		}
	}


	measures.clear();
	measures.reserve(yAxisMeasures.size() * xAxisMeasures.size());
	fillAxisMeasures(xAxisMeasures, p0Measure.x, nxMeasure, hxMeasure);
	fillAxisMeasures(yAxisMeasures, p0Measure.y, nyMeasure, hyMeasure);
	zMeasure = p0Measure.z;
	for (int j = 0; j < yAxisMeasures.size(); j++) {
		for (int i = 0; i < xAxisMeasures.size(); i++) {
			measures.push_back(Measure(Point(xAxisMeasures[i], yAxisMeasures[j], zMeasure), Point (0,0,0)));
		}
	}

	gaussPoints = { -0.774596669, 0.0 , 0.774596669 };
	gaussWeights = { 0.555555556, 0.888888889, 0.555555556 };


	double dim = TASK_DIM * elems.size();
	p.resize(dim);

	// save grid info:
	gridInfo.elemsInX = pointsInX - 1;
	gridInfo.elemsInY = pointsInY - 1;
	gridInfo.elemsInZ = pointsInZ - 1;
	gridInfo.elemsSize = gridInfo.elemsInX * gridInfo.elemsInY * gridInfo.elemsInZ;
	gridInfo.pointsSize = nodes.size();
	gridInfo.yResultsLayersSize = yAxisGrid.size() - 1;
	gridInfo.yMeasureLayersSize = yAxisMeasures.size();
	gridInfo.xMeasureLayersSize = xAxisMeasures.size();

	gridInfo.dx = xAxisGrid[1] - xAxisGrid[0];
	gridInfo.xStart = xAxisGrid.front();
	gridInfo.xEnd = xAxisGrid.back();
	gridInfo.dz = zAxisGrid[1] - zAxisGrid[0];
	gridInfo.zStart = zAxisGrid.front();
	gridInfo.zEndl = zAxisGrid.back();
}

void Task::fillAxisGrid(std::vector<double>& axis, double a, double b, int steps) {
	double step = (b - a) / steps;
	double point = a;

	axis.clear();

	if (axis.empty()) axis.push_back(point);

	for (int i = 0; i < steps; i++) {
		point += step;
		axis.push_back(point);
	}
}

void Task::fillAxisMeasures(std::vector<double>& axis, double a, int steps, double h) {
	double point = a;
	axis.clear();
	for (int i = 0; i <= steps; i++) {
		axis.push_back(point);
		point += h;
	}
}

void Task::getB() {
	for (int i = 0; i < measures.size(); i++) {
		measures[i].B = calculateB(i);
	}
}

Point Task::calculateB(int i) {
	Point res = Point(0, 0, 0);
	for (FiniteElem& elem : elems) {
		Point b = B(i, elem);
		res.x += b.x;
		res.y += b.y;
		res.z += b.z;
	}
	return res;
}

Point Task::B(int i, FiniteElem elem) {
	double sumx = 0;
	double sumy = 0;
	double sumz = 0;
	for (int gx = 0; gx < gaussPoints.size(); gx++) {
		for (int gy = 0; gy < gaussPoints.size(); gy++) {
			for (int gz = 0; gz < gaussPoints.size(); gz++) {
				Point gaussPoint = Point(
					variableChange(gaussPoints[gx], nodes[elem.nodes[0]].x, nodes[elem.nodes[elem.nodes.size() - 1]].x),
					variableChange(gaussPoints[gy], nodes[elem.nodes[0]].y, nodes[elem.nodes[elem.nodes.size() - 1]].y),
					variableChange(gaussPoints[gz], nodes[elem.nodes[0]].z, nodes[elem.nodes[elem.nodes.size() - 1]].z));
				double rj = r(measures[i].point, gaussPoint);
				Point pj = calculateCoords(measures[i].point, gaussPoint);

				sumx += gaussWeights[gx] * gaussWeights[gy] * gaussWeights[gz] / (4.0 * M_PI * pow(rj, 3)) *
					(elem.p.x * (3.0 * pow(pj.x, 2) / pow(rj, 2) - 1.0) +
						elem.p.y * (3.0 * pj.x * pj.y / pow(rj, 2)) +
						elem.p.z * (3.0 * pj.x * pj.z / pow(rj, 2)));

				sumy += gaussWeights[gx] * gaussWeights[gy] * gaussWeights[gz] / (4.0 * M_PI * pow(rj, 3)) *
					(elem.p.y * (3.0 * pow(pj.y, 2) / pow(rj, 2) - 1.0) +
						elem.p.x * (3.0 * pj.x * pj.y / pow(rj, 2)) +
						elem.p.z * (3.0 * pj.y * pj.z / pow(rj, 2)));

				sumz += gaussWeights[gx] * gaussWeights[gy] * gaussWeights[gz] / (4.0 * M_PI * pow(rj, 3)) *
					(elem.p.z * (3.0 * pow(pj.z, 2) / pow(rj, 2) - 1.0) +
						elem.p.x * (3.0 * pj.x * pj.z / pow(rj, 2)) +
						elem.p.y * (3.0 * pj.y * pj.z / pow(rj, 2)));
			}
		}
	}
	return Point(mes(elem) * sumx, mes(elem) * sumy, mes(elem) * sumz);
}

double Task::mes(FiniteElem elem) {

	double a = abs(nodes[elem.nodes[elem.nodes.size() - 1]].x - nodes[elem.nodes[0]].x);
	double b = abs(nodes[elem.nodes[elem.nodes.size() - 1]].y - nodes[elem.nodes[0]].y);
	double c = abs(nodes[elem.nodes[elem.nodes.size() - 1]].z - nodes[elem.nodes[0]].z);
	return a * b * c;
}

double Task::r(Point p1, Point p2) {
	return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2) + pow(p2.z - p1.z, 2));
}

Point Task::calculateCoords(Point p1, Point p2) {
	return Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

double Task::variableChange(double var, double a, double b) {
	return (b + a) / 2.0 + (b - a) / 2.0 * var;
}
