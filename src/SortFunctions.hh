#ifndef SORT_FUNCTIONS_HH
#define SORT_FUNCTIONS_HH

#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"

TH1I* get_TH1I(TDirectory *file, std::string name, std::string title, int nbins, double xlo, double xhi) {
  file->cd();
  TH1I *hist = file->Get<TH1I>(name.c_str());
  if (hist != NULL) { return hist; }
  hist = new TH1I(name.c_str(), title.c_str(), nbins, xlo, xhi);
  return hist;
}

TH2I* get_TH2I(TDirectory *file, std::string name, std::string title, int nbinsx, double xlo, double xhi, int nbinsy, double ylo, double yhi) {
  file->cd();
  TH2I *hist = file->Get<TH2I>(name.c_str());
  if (hist != NULL) { return hist; }
  hist = new TH2I(name.c_str(), title.c_str(), nbinsx, xlo, xhi, nbinsy, ylo, yhi);
  return hist;
}

TH2D* get_TH2D(TDirectory *file, std::string name, std::string title, int nbinsx, double xlo, double xhi, int nbinsy, double ylo, double yhi) {
  file->cd();
  TH2D *hist = file->Get<TH2D>(name.c_str());
  if (hist != NULL) { return hist; }
  hist = new TH2D(name.c_str(), title.c_str(), nbinsx, xlo, xhi, nbinsy, ylo, yhi);
  return hist;
}

void read_cuts(std::string s, TCutG** cuts, int* cut_ids, int &nCuts) {
  FILE *file = fopen(s.c_str(), "ra");
  if (file == NULL) { std::cout << "Cuts file " << s << " does not exist" << std::endl; return; }
  std::stringstream ss;
  char cline[2048];

  int ind = 0;
  while(std::fgets(cline, sizeof cline, file)!=NULL) {
    std::string line(cline);
    if (line.size() == 0) { continue; }
    if (line[0] == '#') { continue; }
    if (line[0] == ';') { continue; }

    ss.clear();
    ss.str(line);

    int id, npoints;
    ss >> id;
    ss >> npoints;
    cuts[id] = new TCutG();
    cut_ids[ind] = id;
    int ct = 0;
    std::cout << "Loading cut ID " << id << std::endl;
    while (ct < npoints) {
      if (std::fgets(cline, sizeof cline, file) == NULL ) { return; }

      std::string line(cline);
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      ss.clear();
      ss.str(line);

      double x, y;
      ss >> x;
      ss >> y;

      cuts[id]->AddPoint(x,y);
      ++ct;
    }
    ++ind;
    ++nCuts;
  }
}

void write_trace(const PIXIE::Measurement &meas, TFile *file, std::string name) {
  file->cd();
  TH1D *trace = new TH1D(name.c_str(), name.c_str(), meas.traceLength, 0, meas.traceLength);
  for (int i=0; i<meas.traceLength; ++i) {
    trace->SetBinContent(i+1,meas.trace[i]);
  }
  trace->Write(); 
}

void write_traces(const PIXIE::Measurement &meas1,
    const PIXIE::Measurement &meas2, TFile *file, std::string name) {
  file->cd();
  TH1D *trace1 = new TH1D((name+"_1").c_str(), (name+"_1").c_str(), 
      meas1.traceLength, 0, meas1.traceLength);
  TH1D *trace2 = new TH1D((name+"_2").c_str(), (name+"_2").c_str(), 
      meas2.traceLength, 0, meas2.traceLength);
  if (meas1.traceLength != meas2.traceLength) {
    std::cout << "Severe warning: write_traces called with unequal trace lengths" << std::endl;
    return;
  }
  for (int i=0; i<meas1.traceLength; ++i) {
    trace1->SetBinContent(i+1,meas1.trace[i]);
    trace2->SetBinContent(i+1,meas2.trace[i]);
  }
  trace1->Write();
  trace2->Write();
}

float interpolate(float *yvals, int n, float x) {
  int ix = (int)x;
  if (ix < 0) { return 0; }
  if (ix >= n-1) { return 0; }
  float y1 = yvals[ix];
  float y2 = yvals[ix+1];
  return (y2-y1)*x + y1 - (y2-y1)*(float)ix;
}

void read_pca(double* basis_vectors[600], std::string bvecname) {
  std::ifstream bvecfile(bvecname);

  std::string line;
  while (std::getline(bvecfile, line)) {
    if (line.size() == 0) { continue; }
    if (line[0] == '#') { continue; }
    if (line[0] == ';') { continue; }
    std::stringstream ss(line);
    int ID;
    int nPoints;
    ss >> ID >> nPoints;

    std::cout << "Loading basis vector for ID " << ID <<", " << nPoints << " points" << std::endl;
    basis_vectors[ID] = new double[nPoints];
    for (int i=0; i<nPoints; ++i) {
      std::getline(bvecfile,line);

      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      double pc0,pc1,pc2,pc3,pc4;
      std::stringstream ss2(line);
      ss2 >> pc0 >> pc1 >> pc2 >> pc3 >> pc4;
      basis_vectors[ID][i] = pc0;
    }
  }
}

TH2D *convert(TH2I *hist, std::string name) {  
  //std::cout << "Converting " << hist->GetName() << " to TH2D " << name << std::endl;
  TH2D *outHist = new TH2D(name.c_str(), hist->GetTitle(),
      hist->GetNbinsX(), hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()),
      hist->GetNbinsY(), hist->GetYaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinUpEdge(hist->GetNbinsY()));
  //std::cout << "X axis " << hist->GetNbinsX() << "   " << hist->GetXaxis()->GetBinLowEdge(1) << "   " << hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()) << std::endl;
  //std::cout << "Y axis " << hist->GetNbinsY() << "   " << hist->GetYaxis()->GetBinLowEdge(1) << "   " << hist->GetYaxis()->GetBinUpEdge(hist->GetNbinsY()) << std::endl;
  for (int i=0; i<=hist->GetNbinsX(); ++i) {
    for (int j=0; j<=hist->GetNbinsY(); ++j) {
      outHist->SetBinContent(i, j, hist->GetBinContent(i, j));
      outHist->SetBinError(i, j, hist->GetBinError(i, j));
    }
  }
  return outHist;
}


#endif
