#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#define _USE_MATH_DEFINES

#include "Lab3.h"

using namespace std;



int main() {

  const float c=300000000.;
  
  TLorentzVector B,pi,k;
  double x,y,z;
  
  TRandom* r1=new TRandom();
  TRandom* gen=new TRandom();
  vector <double> vecxpion,vecypion,veczpion,vecxkaon,vecykaon,veczkaon;
  
  double m_B=5279.; //MeV
  double m_pi=140.; //MeV
  double m_k=500.; //MeV

  double invmpi,invmk,truemass;
  double resol=0.03;
  
  TH1F h1 ("h1","TrueMass",100,4500,6000);
  TH1F h2 ("h2","Scattering Angle pi-k",1000,-100,300);
TH1F h3 ("h3","NewTrueMass",100,4500,6000);
  TCanvas canv("canv","canvas for plotting",1280,1024);
  
  B.SetPxPyPzE(300.,0.,0.,sqrt(m_B*m_B + 300.*300.));
  
  double p_pi=sqrt(pow(m_B,4)+pow(m_k,4)+pow(m_pi,4)-2*m_B*m_B*m_pi*m_pi-2*m_B*m_B*m_k*m_k-2*m_pi*m_pi*m_k*m_k)/(2*m_B);
  double p_k = p_pi;
  double p_k_0,p_pi_0,p_k_meas,p_pi_meas=0.;
  
  
  
  for(int iexp=0;iexp<10000;iexp++) {
    r1->SetSeed(0);
    r1->Sphere(x,y,z,p_pi);
    vecxpion.push_back(x);
    vecypion.push_back(y);
    veczpion.push_back(z);
    
    //pi.SetVect(TVector3(x,y,z));
    pi.SetPxPyPzE(x,y,z,sqrt(p_pi*p_pi + m_pi*m_pi));
    
    r1->SetSeed(0);
    r1->Sphere(x,y,z,p_k);
    
    vecxkaon.push_back(x);
    vecykaon.push_back(y);
    veczkaon.push_back(z);
    
    k.SetPxPyPzE(x,y,z,sqrt(p_k*p_k + m_k*m_k));
    
    p_k_0=sqrt(pi.Px()*pi.Px() + pi.Py()*pi.Py() + pi.Pz()*pi.Pz());
    p_pi_0=sqrt(k.Px()*k.Px() + k.Py()*k.Py() + k.Pz()*k.Pz());
    
    
    
    //cout<<"------->Boosting info"<<endl;
    pi.BoostVector();
    k.BoostVector();
    
    cout<<"\n--------->Boosting per pi e k nel sistema di riferimento del Lab"<<endl;
    pi.Boost(pi.BoostVector());
    k.Boost(k.BoostVector());
    /*pi.Print();
      k.Print();*/

    //invmpi=pi.Gamma()*m_pi; 
    truemass=sqrt(p_pi*p_pi + m_pi*m_pi) + sqrt(p_k*p_k + m_k*m_k);
    cout<<"\nValore truemass: "<<setprecision(10)<<truemass<<endl;
    
    h1.Fill(truemass);
    h2.Fill(ComputeAngle(pi,k));
    
    
    //Parte Seconda----> Non so se si deve fare così

    
    p_pi_meas=gen->Gaus(p_pi_0,p_pi_0*resol);
    p_k_meas=gen->Gaus(p_k_0,p_k_0*resol);
    

    /* k.SetPxPyPzE(x,y,z,sqrt(p_k_meas*p_k_meas + m_k*m_k));
    pi.SetPxPyPzE(x,y,z,sqrt(p_pi_meas*p_pi_meas + m_pi*m_pi));
    k.BoostVector();
    pi.BoostVector();
    pi.Boost(pi.BoostVector());
    k.Boost(k.BoostVector());
    */

    double newtruemass=sqrt(p_k_meas*p_k_meas + m_k*m_k)+sqrt(p_pi_meas*p_pi_meas + m_pi*m_pi);
    h3.Fill(newtruemass);

    cout<<endl<<endl<<endl; 
  }
  
  cout<<"Gamma per il pione: "<<setprecision(4)<<pi.Gamma()<<endl<<endl;
  // cout<<"Beta per il pione: "<<setprecision(9)<<pi.Beta()<<endl<<endl;
  
  cout<<"\nPer il pione:"<<endl;
  for(int i=0;i<vecxpion.size();i++) {
    cout<<"Vettore px: "<<vecxpion[i]<<"\tVettore py: "<<vecypion[i]<<"\tVettore pz: "<<veczpion[i]<<endl;
  } 

  cout<<"\nPer il kaone:"<<endl;
for(int i=0;i<vecxkaon.size();i++) {
    cout<<"Vettore px: "<<vecxkaon[i]<<"\tVettore py: "<<vecykaon[i]<<"\tVettore pz: "<<veczkaon[i]<<endl;
  } 

 cout<<endl<<"Modulo impulso pione: "<<p_pi<<"\tModulo impulso kaone: "<<p_k<<endl;

 cout<<"Per il pione--> Theta= "<<pi.Theta()<<"\tPhi= "<<pi.Phi()<<endl;
  
  delete r1;
  delete gen;
  
  h1.GetXaxis()->SetTitle("m [MeV]");
  h1.GetYaxis()->SetTitle("A.U.");
  h1.Draw();
  canv.SaveAs("TrueMass.pdf");

  h2.GetXaxis()->SetTitle("Angle[degree]");
  h2.GetYaxis()->SetTitle("A.U.");
  h2.Fit("gaus");
  h2.Draw();
  canv.SaveAs("Scattering Angle pi-k.pdf");
  //h1.Write();
  
  h3.GetXaxis()->SetTitle("m [MeV]");
  h3.GetYaxis()->SetTitle("A.U.");
  h3.Fit("gaus");
  h3.Draw();
  canv.SaveAs("NewTrueMass.pdf");
  
  return 0;
}
