/********************************************************************************/
/* Copyright 2009-2011 -- The Regents of the University of California           */
/* This code is provided for research purposes to scientists at non-profit		*/
/*  organizations.  All other use is strictly prohibited.  For further			*/
/*  details please contact University of California, Santa Cruz or				*/
/*	Five3 Genomics, LLC (http://five3genomics.com).								*/
/********************************************************************************/

#include <fstream>
#include <string>
#include <Rcpp.h>
#include "common.h"
#include "evidencesource.h"
using namespace std; //not sure about this one, can be deleted
using namespace Rcpp;
#define THROW(msg) throw std::runtime_error(msg)

EvidenceFactorGen::EvidenceFactorGen(const PropertySet& p) : _params()
{
  _params.reserve(9);
  if (p.hasKey("factorParams")) {
    vector<string> paramsStr;
    Tokenize(p.getStringAs<string>("factorParams"),paramsStr,";");
    if (paramsStr.size() != 9) {
      THROW("must have 9 elements in factorParams");
    }
    for(size_t i = 0; i < paramsStr.size(); i++)
    {
      _params.push_back(atof(paramsStr[i].c_str()));
    }
    return;
  }

  if (!p.hasKey("epsilon") || !p.hasKey("epsilon0")) {
    THROW("must specify either factorParams, or epsilon and epsilon0");
  }

  double epsilon = p.getStringAs<double>("epsilon");
  double epsilon0 = p.getStringAs<double>("epsilon0");
  Real major = 1 - epsilon;
  Real minor = epsilon / 2;
  Real major0 = 1 - epsilon0;
  Real minor0 = epsilon0 / 2;

  bool flip = false;
  if (p.hasKey("reverse")) {
    if (p.getStringAs<string>("reverse") == "true") {
      flip = true;
    } else if (p.getStringAs<string>("reverse") == "false") {
      flip = false;
    } else {
      THROW("The 'reverse' option for evidence must be 'true' or 'false'");
    }
  }

  if (flip) {
    _params.push_back(minor);
    _params.push_back(minor0);
    _params.push_back(major);

    _params.push_back(minor);
    _params.push_back(major0);
    _params.push_back(minor);

    _params.push_back(major);
    _params.push_back(minor0);
    _params.push_back(minor);
  } else {
    _params.push_back(major);
    _params.push_back(minor0);
    _params.push_back(minor);

    _params.push_back(minor);
    _params.push_back(major0);
    _params.push_back(minor);

    _params.push_back(minor);
    _params.push_back(minor0);
    _params.push_back(major);
  }
}

void EvidenceFactorGen::generateValues(const vector< string >& edge_types,
		      vector< Real >& outVals) const
{
  assert(edge_types.size() == 1);
  for (size_t i = 0; i < _params.size(); i++) {
    outVals.push_back(_params[i]);
  }
}

EvidenceSource::EvidenceSource(PropertySet &p, string base) :
  cutoffs(),
  options(p),
  attachPoint(),
  _evidenceFile()
{
  if (p.hasKey("disc"))
    setCutoffs(p.getAs<string>("disc"));
  else
    setCutoffs("-1.3,1.3");

  if (p.hasKey("node"))
    attachPoint = p.getAs<string>("node");
  else
    THROW("EvidenceSource conf. is missing the required property \"node\"");

  if (p.hasKey("suffix")) {
    _suffix = p.getAs<string>("suffix");
    _evidenceFile = base + _suffix;
  } else {
    THROW("EvidenceSource conf. is missing the required property \"suffix\"");
  }
}

void EvidenceSource::setCutoffs(string discLimits)
{
  vector<string> cutoffsStr;
  Tokenize(discLimits,cutoffsStr,";");
  for(size_t i = 0; i < cutoffsStr.size(); i++)
    {
      cutoffs.push_back(atof(cutoffsStr[i].c_str()));
    }
}

int EvidenceSource::discCutoffs  (float x)
{
  if(cutoffs.size() == 0)
    return 0;
  if(x < cutoffs[0])
    return 0;
  size_t i = 1;
  while (i < cutoffs.size())
    {
      if (x < cutoffs[i])
	return i;
      i++;
    }
  return i;
}

double stringToDouble(const string& s) {
  stringstream ss(s);
  double result;
  if (!(ss >> result).eof()) {
    THROW("String " + s + " can not be converted to double.");
  }
  return result;
}

void EvidenceSource::loadFromFile(DataFrame evidenceDf_mRNA, Rcpp::CharacterVector genes, PathwayTab& p,
				  map<string, size_t>& sampleMap,
				  vector<Evidence::Observation>& sampleData)
{
  vector<string> header;
  for(int i = 0; i < genes.size(); i++){
    header.push_back(as<std::string>(genes[i]));
  }
  vector<Var> vars;
  vars.reserve(header.size());

  for (size_t h = 0; h < header.size(); h++) {
	  if(p.getEntityType(header[h]) == "protein") // skip adding evidence if it's not in the pathway
	  {
	    vars.push_back(p.addObservationNode(header[h], attachPoint, _suffix));
	  }
  }

  FactorGenerator* fgen = new EvidenceFactorGen(options);
  p.addFactorGenerator("protein", _suffix, fgen);
  CharacterVector samples = evidenceDf_mRNA["V1"];
  int df_size = evidenceDf_mRNA.size();
  
  for(int i = 0; i < samples.size(); i++){
    vector<string> vals;
    
    string sample = std::to_string(i+1);
    _sampleNames.push_back(sample);

    for(int b = 0; b < genes.size(); b++){
      string temp = "V";
      string col_id = temp + std::to_string(b+2);
      CharacterVector a = evidenceDf_mRNA[col_id];
      vals.push_back(as<std::string>(a[i]));
    }
    for(size_t x = 0; x < vals.size(); x++) {
      if(p.getEntityType(header[x]) != "protein") // skip adding evidence if it's not in the pathway
        continue;
      if(strcmp(vals[x].c_str(),"NA")==0)
        continue;

      double evidence = stringToDouble(vals[x]);
      if (sampleMap.count(sample) == 0) {
        sampleMap[sample] = sampleData.size();
        sampleData.push_back(Evidence::Observation());
      }
      size_t sample_idx = sampleMap[sample];
      sampleData[sample_idx][vars[x]] = discCutoffs(evidence);
    }


  }
  return;
}


void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}
