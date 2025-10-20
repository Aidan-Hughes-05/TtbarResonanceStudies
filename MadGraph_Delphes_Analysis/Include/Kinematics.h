#pragma once

#include <string>
#include <TLorentzVector.h>
#include <vector>

namespace Kinematics {

/**
 * @brief Function to calculate delta phi between two angles.
 * @param phi_1 angle 1.
 * @param phi_2 angle 2.
 * 
 */
double del_phi(double phi_1, double phi_2);

/**
 * @brief This function calculates the invariant mass of a list of particles
 * @param particles A list of TLorentzVector objects pointers
 * 
 */
double Mass(const std::vector<TLorentzVector>& particles);
/**
 * @brief This function returns true if the charge of the particles is consistent with the region name.
 * @param q1
 * @param q2
 * @param regionName
 */
bool isChargeCorrect(const std::string& regionName, float q1, float q2);

/**
 * @brief This function calculates the transverse mass of a particle.
 * @param particle The particle to calculate the transverse mass for.
 * @param MET The missing transverse energy.
 * @return The transverse mass of the particle.
 */
double TransverseMass(const TLorentzVector* particle, const TLorentzVector* MET);

/**
 * @brief Function to see if a jet is inside the rapidity gap spanned by two jets.
 * @param test_jet Jet to test.
 * @param j1 leading jet.
 * @param j2: subleading jet.
 */
int is_inside_jets(TLorentzVector* test_jet,TLorentzVector* j1, TLorentzVector* j2);

/**
 * @brief Function to calculate the minimum delta R between a test particle and a container of particles.
 * @param test_particle: particle to test.
 * @param bool_vector_container: container of booleans to select particles.
 * @param @param jet_container: container of particles to test against.
 */
double min_deltaR(TLorentzVector* test_particle, const std::vector<UInt_t>& bool_vector_container,const std::vector<TLorentzVector*>& jet_container);

/**
 * @brief Get the number of jets in the gap defined by two tagging jets.
 * @param jetsVector The two first jets should be the tagging jets.
 * @param jetsBoolVector The boolean vector indicating if the jet is present in the event or not.
 * @return The number of jets in the gap.
 */
int getNumberOfGapJets(const std::vector<TLorentzVector*>& jetsVector, const std::vector<UInt_t>& jetsBoolVector);

/**
 * @brief Construct the pT balance variable.
 * @param particles The particle container.
 * @return The pT balance of the event.
 */
double getPtBalance(const std::vector<TLorentzVector*>& particles);

/**
 * @brief Construct the signed-centrality variable.
 * @param taggingJet1 The leading tagging jet.
 * @param taggingJet2 The sub-leading tagging jet.
 * @param Tau The tau jet.
 * @param Lepton The lepton.
 * @return The centrality of the event.  
 */
double getSignedCentrality(const TLorentzVector* taggingJet1, const TLorentzVector* taggingJet2, const TLorentzVector* Tau, const TLorentzVector* Lepton);

}