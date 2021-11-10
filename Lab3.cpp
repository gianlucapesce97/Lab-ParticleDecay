#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iomanip>
#include <iostream>
#include <vector>
//#include <cmath>
#define _USE_MATH_DEFINES

#include "Lab3.h"

using namespace std;



int main() {


  //Dichiarazioni Variabili

  const float c=300000000.;  // km/s
  int N=0;
  cout<<"Dimmi il numero di processi di scattering che vuoi analizzare: "<<endl;
  cin>>N;
  
  TLorentzVector B,pi,k;
  double x,y,z=0.;
  
  TRandom* r1=new TRandom();
  TRandom* gen1=new TRandom();
  TRandom* gen2=new TRandom();
  TRandom* gen3=new TRandom();
  TRandom* gen4=new TRandom();
  vector <double> vecxpion,vecypion,veczpion,vecxkaon,vecykaon,veczkaon;
  
  double m_B=5279.; //MeV
  double m_pi=140.; //MeV
  double m_k=500.; //MeV

  double invmpi,invmk,truemass=0.;
  double resol1=0.03;
  double resol2=0.01;
  double resol3=0.05;
  double resol4=0.10;
  
  TH1F h1 ("h1","TrueMass",100,4500,6000);
  TH1F h2 ("h2","Scattering Angle pi-k",1000,2,5);
  TH1F h3 ("h3","NewTrueMass-Resolution 3%",100,2000,7000);
  TH1F h4 ("h4","NewTrueMass-Resolution 1%",100,2000,7000);
  TH1F h5 ("h5","NewTrueMass-Resolution 5%",100,2000,7000);
  TH1F h6 ("h6","NewTrueMass-Resolution 10%",100,2000,7000);
  
  TCanvas canv("canv","canvas for plotting",1280,1024);
  //TCanvas canv2("canv2","SuperimposingIstograms",1280,1024);
  
  B.SetPxPyPzE(300.,0.,0.,sqrt(m_B*m_B + 300.*300.)); 
  
  double p_pi=sqrt(pow(m_B,4)+pow(m_k,4)+pow(m_pi,4)-2*m_B*m_B*m_pi*m_pi-2*m_B*m_B*m_k*m_k-2*m_pi*m_pi*m_k*m_k)/(2*m_B);
  double p_k = p_pi;
  double p_k_0,p_pi_0,p_k_meas1,p_pi_meas1,p_k_meas2,p_pi_meas2,p_k_meas3,p_pi_meas3,p_k_meas4,p_pi_meas4=0.;
  
  
  
  for(int iexp=0;iexp<N;iexp++) {
    r1->SetSeed(0);
    r1->Sphere(x,y,z,p_pi);
    vecxpion.push_back(x);
    vecypion.push_back(y);
    veczpion.push_back(z);
    
    
    
    vecxkaon.push_back(-x);
    vecykaon.push_back(-y);
    veczkaon.push_back(-z);
    
    pi.SetPxPyPzE(+x,+y,+z,sqrt(p_pi*p_pi + m_pi*m_pi));
    k.SetPxPyPzE(-x,-y,-z,sqrt(p_k*p_k + m_k*m_k));
    
    cout<<"Informazioni su pi prima del boost\n"<<endl<<endl;
    pi.Print();
    
    cout<<"------->Boosting info"<<endl;
    pi.BoostVector();
    k.BoostVector();
    
    //--------->Boosting per pi e k nel sistema di riferimento del Lab

    pi.Boost(B.BoostVector());
    k.Boost(B.BoostVector());
    pi.Print();
    k.Print();
    
    
    
   
    cout<<"Informazioni su pi dopo il boost\n"<<endl<<endl;
    pi.Print();
    k.Print();
    
    //Questa è la truemass calcolata nel SDR del CDM. Ovviamente viene precisamente 5279 MeV
    //truemass=sqrt(p_pi*p_pi + m_pi*m_pi) + sqrt(p_k*p_k + m_k*m_k);
    truemass=sqrt(pow(pi.E()+k.E(),2)-sqrt(pow(pi.Px()+k.Px(),2) + pow(pi.Py()+k.Py(),2) + pow(pi.Pz()+k.Pz(),2)));
    h1.Fill(truemass);


    h2.Fill(ComputeAngle(pi,k)); 

    
    
    //cout<<"\nValore truemass: "<<setprecision(10)<<truemass<<endl;
    
    
    
    
    
    cout<<endl<<endl<<endl; 
  }
  
  
  
  //Varie Stampe di informazioni
  
  cout<<"Gamma per il pione: "<<setprecision(4)<<pi.Gamma()<<"\tGamma per il kaone: "<<setprecision(4)<<k.Gamma()<<endl<<endl;
  cout<<"Beta per il pione: "<<setprecision(9)<<pi.Beta()<<"\tBeta per il kaone: "<<setprecision(9)<<k.Beta()<<endl<<endl;
  
  cout<<"\nPer il pione:"<<endl;
  for(int i=0;i<vecxpion.size();i++) {
    cout<<"Vettore px: "<<vecxpion[i]<<"\tVettore py: "<<vecypion[i]<<"\tVettore pz: "<<veczpion[i]<<endl;
  } 
  
  cout<<"\nPer il kaone:"<<endl;
  for(int i=0;i<vecxkaon.size();i++) {
    cout<<"Vettore px: "<<vecxkaon[i]<<"\tVettore py: "<<vecykaon[i]<<"\tVettore pz: "<<veczkaon[i]<<endl;
  } 
  
  cout<<endl<<"Modulo impulso pione: "<<p_pi<<"\tModulo impulso kaone: "<<p_k<<endl;
  
  //cout<<"Per il pione--> Theta= "<<pi.Theta()<<"\tPhi= "<<pi.Phi()<<endl;
  
  


  //-------->Parte Seconda
  //Detector di rilevazione impulso con distribuzione gaussiana con risoluzione del 3%/1%/5%/10%
  
  for (int iexp=0;iexp<N;iexp++) {
    
    r1->SetSeed(0);
    r1->Sphere(x,y,z,p_pi);
    
    pi.SetPxPyPzE(x,y,z,sqrt(p_pi*p_pi + m_pi*m_pi));
    
    k.SetPxPyPzE(-x,-y,-z,sqrt(p_k*p_k + m_k*m_k));
    
    
    
    pi.Boost(B.BoostVector());
    k.Boost(B.BoostVector());
    
    p_k_0=sqrt(pi.Px()*pi.Px() + pi.Py()*pi.Py() + pi.Pz()*pi.Pz());
    p_pi_0=sqrt(k.Px()*k.Px() + k.Py()*k.Py() + k.Pz()*k.Pz());
    
    
    
    p_pi_meas1=gen1->Gaus(p_pi_0,p_pi_0*resol1);
    p_k_meas1=gen1->Gaus(p_k_0,p_k_0*resol1);

    p_pi_meas2=gen2->Gaus(p_pi_0,p_pi_0*resol2);
    p_k_meas2=gen2->Gaus(p_pi_0,p_pi_0*resol2);
    
    p_pi_meas3=gen3->Gaus(p_pi_0,p_pi_0*resol3);
    p_k_meas3=gen3->Gaus(p_pi_0,p_pi_0*resol3);
    
    p_pi_meas4=gen4->Gaus(p_pi_0,p_pi_0*resol4);
    p_k_meas4=gen4->Gaus(p_pi_0,p_pi_0*resol4);
    
    double newtruemass1=sqrt(p_k_meas1*p_k_meas1 + m_k*m_k)+sqrt(p_pi_meas1*p_pi_meas1 + m_pi*m_pi);
    double newtruemass2=sqrt(p_k_meas2*p_k_meas2 + m_k*m_k)+sqrt(p_pi_meas2*p_pi_meas2 + m_pi*m_pi);
    double newtruemass3=sqrt(p_k_meas3*p_k_meas3 + m_k*m_k)+sqrt(p_pi_meas3*p_pi_meas3 + m_pi*m_pi);
    double newtruemass4=sqrt(p_k_meas4*p_k_meas4 + m_k*m_k)+sqrt(p_pi_meas4*p_pi_meas4 + m_pi*m_pi);

    h3.Fill(newtruemass1);
    h4.Fill(newtruemass2);
    h5.Fill(newtruemass3);
    h6.Fill(newtruemass4);
    
  }


  
  






    //-------> Grafici
  
  h1.GetXaxis()->SetTitle("m [MeV]");
  h1.GetYaxis()->SetTitle("A.U.");
  h1.Draw();
  canv.SaveAs("TrueMass.pdf");
  
  h2.GetXaxis()->SetTitle("Angle[rad]");
  h2.GetYaxis()->SetTitle("A.U.");
  h2.Draw();
  canv.SaveAs("Scattering Angle pi-k.pdf");
  //h1.Write();
  
  h3.GetXaxis()->SetTitle("m [MeV]");
  h3.GetYaxis()->SetTitle("A.U.");
  h3.Fit("gaus");
  h3.Draw("");
  canv.SaveAs("NewTrueMass-Resolution 3%.pdf");
  
  h4.GetXaxis()->SetTitle("m [MeV]");
  h4.GetYaxis()->SetTitle("A.U.");
  h4.Fit("gaus");
  h4.Draw("pe");
  canv.SaveAs("NewTrueMass-Resolution 1%.pdf");
  
  h5.GetXaxis()->SetTitle("m [MeV]");
  h5.GetYaxis()->SetTitle("A.U.");
  h5.Fit("gaus");
  h5.Draw("pe");
  canv.SaveAs("NewTrueMass-Resolution 5%.pdf");
  
  h6.GetXaxis()->SetTitle("m [MeV]");
  h6.GetYaxis()->SetTitle("A.U.");
  h6.Fit("gaus");
  h6.Draw("pe");
  canv.SaveAs("NewTrueMass-Resolution 10%.pdf");
  
  






  
  //----->Eliminazioni
  delete r1;
  delete gen1;
  delete gen2;
  delete gen3;
  delete gen4;
  
  return 0;
}
