#ifndef __ADRC_H__
#define __ADRC_H__
#include <stdint.h>
#include "Ano_Math.h"
#include "Math.h"
//#define h 0;
//#define h0 0;
//#define r 0;
//#define bata1 0;
//#define bata2 0;
//#define bata3 0;


typedef struct ADRC
{
	float Z1,Z2,Z3,e,e1,e2,u,u0,v1,v2,v;

	
} _ADRC_;

void ADRC_Init(_ADRC_ *sptr);
float ADRC_han(float vin,float y , _ADRC_ *sptr);

float ADRC_sun(float vin,float y , _ADRC_ *sptr);
#endif
