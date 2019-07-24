#ifndef __primitiveh__
#define __primitiveh__

typedef struct ring
{
	int *link;
	int len;
	struct ring *next;
} ring;

struct primitive
{
	ring *rings;
};



#endif
