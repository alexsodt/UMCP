#ifndef POINT_H
#define POINT_H

struct Point {
	double x, y, z;
};

bool operator==(const Point& a, const Point& b);
Point operator-(const Point& a, const Point& b);
Point operator+(const Point& a, const Point& b);
Point operator*(const Point& a, const double b);
Point operator/(const Point& a, const double b);
Point negative(const Point& a);

#ifdef min_distance_cpp

bool operator==(const Point& a, const Point& b) {
    if ((a.x == b.x) && (a.y == b.y) && (a.z == b.z)) {
		return true;
	}
	return false;
}

Point operator-(const Point& a, const Point& b) {
	Point new_point;
	new_point.x = a.x - b.x;
	new_point.y = a.y - b.y;
	new_point.z = a.z - b.z;
	return new_point;
}

Point operator+(const Point& a, const Point& b) {
	Point new_point;
	new_point.x = a.x + b.x;
	new_point.y = a.y + b.y;
	new_point.z = a.z + b.z;
	return new_point;
}

Point operator*(const Point& a, const double b) {
	Point new_point;
	new_point.x = a.x * b;
	new_point.y = a.y * b;
	new_point.z = a.z * b;
	return new_point;
}

Point operator/(const Point& a, const double b) {
	Point new_point;
	new_point.x = a.x / b;
	new_point.y = a.y / b;
	new_point.z = a.z / b;
	return new_point;
}


Point negative(const Point& a) {
	Point neg_point;
	neg_point.x = -a.x;
	neg_point.y = -a.y;
	neg_point.z = -a.z;
	return neg_point;
}

#endif

#endif
