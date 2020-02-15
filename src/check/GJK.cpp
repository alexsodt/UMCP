#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "functions.h"

using namespace std;
/*
double** a;
int numA;
double** b;
int numB;

void input_data(const char*& file1, const char*& file2) {
    ifstream ifs;
    ifs.open(file1);
    string line;
    getline(ifs, line);
    istringstream hold(line);
    hold >> numA;

    // comment line
    getline(ifs, line);
    istringstream hold2(line);
    a = new double*[numA];
	for (int i = 0; i < numA; i++) {
		a[i] = new double[3];
	}

    for (int i = 0; i < numA; i++) {
        getline(ifs, line);
        istringstream hold3(line);
        string name;
        hold3 >> name;
        hold3 >> a[i][0];
        hold3 >> a[i][1];
        hold3 >> a[i][2];
    }   

    ifs.close();

    ifstream ifs2;
    ifs2.open(file2);
    getline(ifs2, line);
    istringstream hold4(line);
    hold4 >> numB;

    // comment line
    getline(ifs2, line);
    istringstream hold5(line);
    b = new double*[numB];
	for (int i = 0; i < numB; i++) {
		b[i] = new double[3];
	}

    for (int i = 0; i < numB; i++) {
        getline(ifs2, line);
        istringstream hold6(line);
        string name;
        hold6 >> name;
        hold6 >> b[i][0];
        hold6 >> b[i][1];
        hold6 >> b[i][2];
    }   

    ifs2.close();
}

void write_Minkowski() {
	ofstream file;
	file.open("Minkowski_Shape.xyz");
	file << num << "\n";
	file << "//comment here" << "\n";
	for (int i = 0; i < num; i++) {
		file << "C" << " " <<  md[i].x << " " << md[i].y << " " << md[i].z << "\n";
	}
	file.close();
}
*/

bool gjk_algorithm(double* a, int numA, double* b, int numB) {
	Functions test(a, b, numA, numB);

	for (int i = 0; i < numA; i++) {
//		cout << a[3*i + 0] << " " << a[3*i + 1] << " " << a[3*i + 2] << endl;
	}
//	cout << "next" << endl;
	for (int i = 0; i < numB; i++) {
//		cout << b[3*i + 0] << " " << b[3*i + 1] << " " << b[3*i + 2] << endl;
	}
	bool result = test.run_GJK();
//	cout << "hi" << result << endl;
	return result;
}

/*
int main(int argc, const char* argv[]) {

	const char* input1;
	const char* input2;
	if (argc > 2) {
		input1 = argv[1];
		input2 = argv[2];
	} else {
		return -1;
	}

	input_data(input1, input2);
	cout << test_GJK(a, b, numA, numB) << endl;

}
*/
