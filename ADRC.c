#include"ADRC.h"

float h=0;
float h0=0;
float r=0;
float b=0;
float b0=0;

float delta = 0;
float belta01 = 0;
float belta02 = 0;
float belta03 = 0;

/**************NLSEF*******/
float alpha1 = 0.5,//
      alpha2 = 0.25,//
      belta1 = 0,//
      belta2 = 0;//
float uu = 0;
int Sign_ADRC(float Input)
{
    int output=0;
    if(Input>1E-6) output=1;
    else if(Input<-1E-6) output=-1;
    else output=0;
    return output;
}

int Fsg_ADRC(float x,float d)
{
  int output=0;
  output=(Sign_ADRC(x+d)-Sign_ADRC(x-d))/2;
  return output;
}

float sat(float x , float  delta)
{
   if(ABS(x)>delta)
    return Sign_ADRC(x);
	 else 
		 return x/delta;
}


float fhan(float x1,float x2,float r,float h)
{

	float deltaa  =0,
		  deltaa0 =0,
	      y       =0,
	      a0      =0,
	      a       =0,
	      fhan    =0;
	
	deltaa = r*h;
	deltaa0 = deltaa*h;
	y=x1+x2*h;
	a0 = sqrtf(deltaa*deltaa+8*r*fabsf(y));
	if(fabsf(y)<=deltaa0)
		a=x2+y/h;
	else
		a=x2+0.5*(a0-deltaa)*Sign_ADRC(y);
	if(fabsf(a)<=deltaa)
		fhan = -r*a/deltaa;
	else
		fhan = -r*Sign_ADRC(a);
	
  return fhan;
}

float fsun(float x1,float x2,float r,float h) // 
{

	float deltaa  =0,
		  deltaa0 =0,
	      y       =0,
	      //a0      =0,
	      
	      fsun    =0,
	      kpie   = 0,
	      k = 0;
	
	deltaa = r*h;
	deltaa0 = deltaa*h;
	y=x1+x2*h;
	kpie = 0.5 * (1+sqrt(1+8*fabsf(y)/deltaa0));
	k = Sign_ADRC(kpie - (int)kpie)+(int)kpie;
	if(fabsf(y)<=deltaa0)
		fsun = -r * sat (x2 + y/h , deltaa); 
	else
		fsun = -r * sat((1-0.5*k)* Sign_ADRC(y) - (x1+k*h*x2)/((k-1)*deltaa0),1);

	
  return fsun;
}

float fal(float e,float alpha,float delta)
{
  float result = 0,fabsf_e = 0;
  
  fabsf_e = fabsf(e);
  
  if(delta>=fabsf_e)
    result = e/powf(delta,1.0-alpha);
  else //if(delta<fabsf_e)
    result = powf(fabsf_e,alpha)*Sign_ADRC(e);
 
 return result;     
}


void TD_process(_ADRC_ *sptr)
{

	
  sptr->v1 = sptr->v1+sptr->v2*h;
	sptr->v2 = sptr->v2+sptr->u*h;
	
}
void ADRC_Init(_ADRC_ *sptr)
{
	sptr->v = 0;
  sptr->v1 = 0;
	sptr->v2 = 0;
	sptr->u = 0;
	sptr->u0 = 0;
	sptr->e = 0;
	sptr->e1 = 0;
	sptr->e2 = 0;
  sptr->Z1 = 0;
  sptr->Z2 = 0;
  sptr->Z3 = 0;


}

float ADRC_han(float vin,float y , _ADRC_ *sptr)
{
   //the process of TD
  sptr->v1 = sptr->v1+sptr->v2*h0;
	sptr->v2 = sptr->v2+sptr->u*h0;
	sptr->u = fhan(sptr->v1 - vin,sptr->v2, r, h);
	////ESO
	sptr->e = sptr->Z1 - y;
  sptr->Z1 = sptr->Z1 + h*(sptr->Z2-belta01*sptr->e);
  sptr->Z2 = sptr->Z2 + h*(sptr->Z3-belta02*fal(sptr->e,alpha1,delta)+b*sptr->u);
  sptr->Z3 = sptr->Z3 + h*(-belta03*fal(sptr->e,alpha2,delta));
	//NLSEF
	sptr->e1 = sptr->v1 - sptr->Z1;
  sptr->e2 = sptr->v2 - sptr->Z2;
  
  sptr->u0 = belta1*fal(sptr->e1,alpha1,delta) + belta2*fal(sptr->e2,alpha2,delta);//??0<alpha1<1<alpha2
  
  sptr->u = sptr->u0 - sptr->Z3/b;
 
	return sptr->u;
}

float ADRC_sun(float vin,float y , _ADRC_ *sptr)
{
   //the process of TD
  sptr->v1 = sptr->v1+sptr->v2*h0;
	sptr->v2 = sptr->v2+sptr->u*h0;
	sptr->u = fsun(sptr->v1 - vin,sptr->v2, r, h);
	////ESO
	sptr->e = sptr->Z1 - y;
  sptr->Z1 = sptr->Z1 + h*(sptr->Z2-belta01*sptr->e);
  sptr->Z2 = sptr->Z2 + h*(sptr->Z3-belta02*fal(sptr->e,alpha1,delta)+b*sptr->u);
  sptr->Z3 = sptr->Z3 + h*(-belta03*fal(sptr->e,alpha2,delta));
	//NLSEF
	sptr->e1 = sptr->v1 - sptr->Z1;
  sptr->e2 = sptr->v2 - sptr->Z2;
  
  sptr->u0 = belta1*fal(sptr->e1,alpha1,delta) + belta2*fal(sptr->e2,alpha2,delta);//??0<alpha1<1<alpha2
  
  sptr->u = sptr->u0 - sptr->Z3/b;
 
	return sptr->u;
}
//之前出了点问题，本版本新增fsun函数替换fhan函数，也是按公式扒的，期待验证

