#pragma once
#include <string>
#include <vector>
#include <filesystem>

#include <TChain.h>


void PartonHistograms(const std::vector<std::filesystem::path>& filesToProcess, const std::string& outputFileName, float lumiWeight, bool doWeights = false);
void PartonNtuples(const std::vector<std::filesystem::path>& filesToProcess, const std::string& outputFileName, float lumiWeight, bool doWeights = false);
void DelphesHistograms(const std::string& fileName, const std::string& outputFileName, bool doWeights = false);

float CalculateLumiWeight(const std::vector<std::filesystem::path>& filesToProcess);
void getChainFromUserInput(const std::vector<std::filesystem::path>& filesToProcess, TChain& chain);