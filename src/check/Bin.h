#ifndef BIN_H
#define BIN_H

struct Bin {
	int *tri_list;
	int num_tri;
	int space;
	double search[3];
	double vmin, vmax;
};

#endif
