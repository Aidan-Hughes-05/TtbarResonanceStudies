#pragma once
#include <string>
#include <vector>
#include <filesystem>

#include <TChain.h>

typedef std::map<std::string, std::tuple<float, long, float>> MetadataObject;

void PartonHistograms(const std::vector<std::filesystem::path>& filesToProcess, const std::string& outputFileName, const MetadataObject& metadataObject, bool doWeights = false);
void PartonNtuples(const std::vector<std::filesystem::path>& filesToProcess, const std::string& outputFileName, const MetadataObject& metadataObject, bool doWeights = false);
void DelphesHistograms(const std::string& fileName, const std::string& outputFileName, bool doWeights = false);

MetadataObject FormMetadataObject(std::string pathToMetadataFile);
float CalculateLumiWeight(const std::filesystem::path& filePath, const MetadataObject& metadataObject);