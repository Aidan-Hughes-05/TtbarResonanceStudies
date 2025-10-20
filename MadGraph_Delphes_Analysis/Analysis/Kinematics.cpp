#include "Kinematics.h"
#include "Constants.h"
#include "Logger.h"
#include <TMath.h>
#include <TLorentzVector.h>
#include <vector>
#include <algorithm>
#include <string>

namespace Kinematics{

double del_phi(double phi_1, double phi_2){
    double pi=TMath::Pi();
    double phi_1_norm, phi_2_norm;
    if (phi_1<0.0){
        phi_1_norm=phi_1+2*pi;
    }else {
        phi_1_norm=phi_1;
    }

    if (phi_2<0.0){
        phi_2_norm=phi_2+2*pi;
    }else {
        phi_2_norm=phi_2;
    }
    double delta=std::abs(phi_1_norm-phi_2_norm);
    if (delta>pi){
        delta=2*pi-delta;
        delta=std::abs(delta);
    }

    return delta;
}

double Mass(const std::vector<TLorentzVector>& particles)
{
    double sum{0};
    for (const auto& particle1 : particles)
    {
        for (const auto& particle2 : particles)
        {
            if (particle1 == particle2) { continue;} // Skip the same particle
            sum += particle1.Dot(particle2);
        }
    }
    return sqrt(sum);
}

bool isChargeCorrect(const std::string& regionName, float q1, float q2){
    if (regionName.size() < 2){
        Constants::g_LOG(LogLevel::ERROR, "Region name is not valid.");
        exit(1);
    }

    bool isOS = regionName.find("OS") != std::string::npos;
    bool isSS = regionName.find("SS") != std::string::npos;

    if (isOS) return q1!=q2;
    if (isSS) return q1==q2;
    else {
        Constants::g_LOG(LogLevel::ERROR, "Region name is not valid. The region should end with OS or SS.");
        exit(1);
    }
}

double TransverseMass(const TLorentzVector* particle, const TLorentzVector* MET)
{
    double mt = sqrt(2*particle->Pt()*MET->Pt() * (1-cos(del_phi(particle->Phi(),MET->Phi()))));
    return mt;
}

int is_inside_jets(TLorentzVector* test_jet,TLorentzVector* j1, TLorentzVector* j2){
  double delta_y_j1j2 = abs(j1->Rapidity()-j2->Rapidity());
  double delta_y_j1test = abs(j1->Rapidity()-test_jet->Rapidity());
  double delta_y_j2test = abs(j2->Rapidity()-test_jet->Rapidity());
  if(delta_y_j1test>delta_y_j1j2 || delta_y_j2test>delta_y_j1j2) return 0;
  return 1;
}

double deltaR(TLorentzVector* test_particle, const std::vector<UInt_t>& bool_vector_container,const std::vector<TLorentzVector*>& jet_container){

  std::vector<double> delta_Rs{};

  for (size_t index{0};index<jet_container.size();index++){
    if (bool_vector_container[index]!=0){
      delta_Rs.push_back(jet_container[index]->DeltaR(*test_particle));
    }
    else {break;}
  }

  double min_dR=*std::min_element(delta_Rs.begin(),delta_Rs.end());
  return min_dR;
}

int getNumberOfGapJets(const std::vector<TLorentzVector*>& jetsVector, const std::vector<UInt_t>& jetsBoolVector){
    // First check that the two tagging jets are present in the event.
    if (jetsBoolVector[0] == 0 || jetsBoolVector[1] == 0){
        Constants::g_LOG(LogLevel::ERROR, "Function called with less than two particles to compare!");
        Constants::g_LOG(LogLevel::ERROR, "At least one of the two tagging jets are not present in the event.");
        exit(1);
    }

    // Count the number of jets in the gap.
    int nGapJets{0};
    for (size_t index{2};index<jetsVector.size();index++){
        if (jetsBoolVector[index] == 0) continue;
        nGapJets += is_inside_jets(jetsVector[index], jetsVector[0], jetsVector[1]);
    }
    return nGapJets;
}

double getPtBalance(const std::vector<TLorentzVector*>& particles){
    double scalarSum{0};
    TLorentzVector vectorSum{};
    for (const auto& particle : particles){
        scalarSum += particle->Pt();
        vectorSum += *particle;
    }

    // Calculate the pT balance.
    double pt_bal= vectorSum.Pt()/scalarSum;
    return pt_bal;
}

double getSignedCentrality(const TLorentzVector* taggingJet1, const TLorentzVector* taggingJet2, const TLorentzVector* Tau, const TLorentzVector* Lepton){
    double lepton_xi = ((*Tau) + (*Lepton)).Rapidity();
    double dijet_xi = taggingJet1->Rapidity() + taggingJet2->Rapidity();
    double delta_y = (taggingJet1->Rapidity() - taggingJet2->Rapidity());
    double centrality = (lepton_xi - 0.5*dijet_xi)/delta_y;
    return centrality;
}
}