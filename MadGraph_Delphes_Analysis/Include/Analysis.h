#pragma once
#include <string>

void PartonAnalysisHistograms(const std::string& fileName, const std::string& outputFileName, float lumiWeight, bool doWeights = false);

void DelphesAnalysis(const std::string& fileName, const std::string& outputFileName, bool doWeights = false);

float CalculateLumiWeight(const std::string& fileName);