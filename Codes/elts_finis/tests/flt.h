#ifndef __FLT_H__
#define __FLT_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct Flt{
	double 	mnt;
	int 		exp;
} Flt;

Flt*	crtFlt(double d){
	Flt* c = malloc(sizeof(Flt));
	c->exp = floor(log(d) / log(10));
	c->mnt = d / pow(10, c->exp);
	return c;
};

Flt* mltFlt(Flt* c1, Flt* c2){
	Flt* c = malloc(sizeof(Flt));
	c->exp = c1->exp + c2->exp;
	int tmp = floor(log(c1->mnt * c2->mnt) / log(10));
	c->mnt = c1->mnt * c2->mnt / pow(10, tmp);
	c->exp += tmp;
	return c;
};

Flt* divFlt(Flt* c1, Flt* c2){
	Flt* c = malloc(sizeof(Flt));
	c->exp *= -c2->exp;
	c->mnt = 1 / c2->mnt;
	return mltFlt(c1, c);
};

Flt* addFlt(Flt* c1, Flt* c2){
	if(c2->exp > c1->exp){
		Flt* tmp = c1;
		c1 = c2;
		c2 = tmp;
	}
	Flt* c = malloc(sizeof(Flt));
	c->exp = c1->exp;
	int tmp = floor(log(c1->mnt + c2->mnt * pow(10, c2->exp - c1->exp)) / log(10));
	c->mnt = (c1->mnt + c2->mnt * pow(10, c2->exp - c1->exp)) / tmp;
	c->exp += tmp;
	return c;
}

Flt* subFlt(Flt* c1, Flt* c2){
	Flt* c = malloc(sizeof(Flt));
	c->exp = c2->exp;
	c->mnt = -c2->mnt;
	return addFlt(c1, c);
}

void prtFlt(Flt* c){
	printf("%f e %d", c->mnt, c->exp);
};

#endif
