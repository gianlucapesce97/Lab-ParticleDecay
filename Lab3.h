#ifndef Lab3_h
#define Lab3_h

#include "TLorentzVector.h"
#include <iostream>
#include "TMath.h"
#include <cmath>
#define _USE_MATH_DEFINES

using namespace std;



double ScalarProduct (const TLorentzVector& tn,const TLorentzVector& ln) {

  double sp=tn.Px()*ln.Px() + tn.Py()*ln.Py() + tn.Pz()*ln.Pz();
  //cout<<"Prodotto scalare: "<<sp<<endl;
  return sp;
}


double ComputeAngle (const TLorentzVector& tn, const  TLorentzVector& ln) {
  double sp=ScalarProduct(tn,ln);
  double mag1=sqrt(pow(tn.Px(),2)+pow(tn.Py(),2)+pow(tn.Pz(),2));
  double mag2=sqrt(pow(ln.Px(),2)+pow(ln.Py(),2)+pow(ln.Pz(),2));
  double angle=acos(sp/(mag1*mag2));
  cout<<"Angolo tra i due vettori: "<<angle<<endl;
  return angle*(180/M_PI);
}


#endif
