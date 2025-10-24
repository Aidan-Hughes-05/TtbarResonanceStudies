#pragma once
#include <string>

void PartonHistograms(const std::string& fileName, const std::string& outputFileName, float lumiWeight, bool doWeights = false);
void PartonNtuples(const std::string& fileName, const std::string& outputFileName, float lumiWeight, bool doWeights = false);

void DelphesHistograms(const std::string& fileName, const std::string& outputFileName, bool doWeights = false);

float CalculateLumiWeight(const std::string& fileName);