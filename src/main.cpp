#include <Rcpp.h>
#include <gmpxx.h>
//#include <mpirxx.h>

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "properties.h"
#include "smallset.h"
#include "util.h"
#include "alldai.h"

#include <sys/time.h>
//#include <sys/resource.h>

#include "common.h"
#include "configuration.h"
#include "evidencesource.h"

using namespace std;
using namespace dai;
//is it because of objects???
using namespace Rcpp;
#define GENE_COL 0
#define MOLECULE_COL 1
#define EVIDENCE_COL 2

#define VAR_DIM 3

typedef std::map< std::string, SmallSet< std::string > > EMStep;
typedef std::vector< EMStep > EMSteps;


inline double log10odds(double post,double prior)
{
  return std::log( ( post / (1.0 - post) )
                     / ( prior / (1.0 - prior) ) )
  / log(10);
}

// returns a map of all the nodes in a single subnet, for testing against
map< long, bool > nodesInSingleNet(vector< Factor > factors)
{
  // start with a single (random) node, and follow it out, marking down all the nodes we see
  map< long, bool > nodesSeen;
  vector< long > nodesToCheck;
  vector< Factor >::iterator factorIter = factors.begin();
  const VarSet tmpVars = factorIter->vars();
  vector< Var >::const_iterator tmpVarsIter;
  for (tmpVarsIter = tmpVars.begin(); tmpVarsIter != tmpVars.end(); ++tmpVarsIter) {
    const long varLabel = tmpVarsIter->label();
    nodesToCheck.push_back(varLabel);
    nodesSeen[varLabel] = true;
  }
  while(nodesToCheck.size() > 0)
  {
    long currNode = nodesToCheck.back();
    nodesToCheck.pop_back();
    for ( factorIter = factors.begin(); factorIter != factors.end(); ++factorIter) {
      const VarSet tmpVars = factorIter->vars();
      //vector< Var >::const_iterator tmpVarsIter = tmpVars.begin();
      for (tmpVarsIter = tmpVars.begin(); tmpVarsIter != tmpVars.end(); ++tmpVarsIter) {
        const long varLabel = tmpVarsIter->label();
        if(varLabel == currNode)
        {
          vector< Var >::const_iterator tmpVarsIter2 = tmpVars.begin();
          for ( ; tmpVarsIter2 != tmpVars.end(); ++tmpVarsIter2) {
            const long varLabel2 = tmpVarsIter2->label();
            if(!nodesSeen[varLabel2])
            {
              nodesToCheck.push_back(varLabel2);
              nodesSeen[varLabel2] = true;
            }
          }
        }
      }
    }

  }
  return nodesSeen;
}

DataFrame outputFastaPerturbations(DataFrame& df, string sampleName, InfAlg* prior, InfAlg* sample,
                              FactorGraph& fg, map<long,string>& activeNodes)
{
  CharacterVector V1 = df["V1"];
  NumericVector V2 = df["V2"];
  // cout << "> " << sampleName;
  // cout << " loglikelihood=" << (sample->logZ() - prior->logZ())
  //      << endl;
  for (size_t i = 0; i < fg.nrVars(); ++i)
  {
    const Var& v = fg.var(i);
    if(activeNodes.count(v.label()) == 0)
      continue;
    V1.push_back(activeNodes[v.label()]);
	//out << activeNodes[v.label()];
    Factor priorBelief = prior->belief(v);
    Factor belief = sample->belief(v);
    vector<double> priors;
    vector<double> posteriors;
    bool beliefEqualOne = false;

    for (size_t j = 0; j < belief.nrStates(); ++j)
    {
      if(belief[j] == 1 || priorBelief[j] == 1)
      {
        beliefEqualOne = true;
        break;
      }
      priors.push_back(priorBelief[j]);
      posteriors.push_back(belief[j]);
    }

    if(beliefEqualOne)
      V2.push_back(NA_REAL);
    else
    {
      double down = log10odds(posteriors[0],priors[0]);
      double nc = log10odds(posteriors[1],priors[1]);
      double up = log10odds(posteriors[2],priors[2]);

      if (nc > down && nc > up)
        V2.push_back(0);
      else if (down > up)
        //out << (-1.0*down);
		V2.push_back(-1.0*down);
      else
        //out << up;
		V2.push_back(up);
    }
    //out << endl;
  }

  //V1.push_back("GATA");
  //V2.push_back(5.23);
  DataFrame dfr = DataFrame::create( Named("V1") = V1,
                                    Named("V2") = V2);
  return dfr;
}
// [[Rcpp::export]]
DataFrame JTinf(DataFrame nodesDf, DataFrame edgesDf, DataFrame evidenceDf_mRNA, CharacterVector genes)
{



  // /////////////////////////////////////////////////
  // Load pathway
  PropertySet pathwayOptions;
  PathwayTab pathway = PathwayTab::create(nodesDf, edgesDf, pathwayOptions);// ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM 

  // /////////////////////////////////////////////////
  // Read in evidence
  vector<EvidenceSource> evid;
  map<string,size_t> sampleMap;
  vector<Evidence::Observation> sampleData;

  string evidOps = "[suffix=_mRNA.tab,node=mRNA,disc=-1.3;1.3,epsilon=0.01,epsilon0=0.2]";
  PropertySet evidenceOptions;
  evidenceOptions.fromString(evidOps);
  
  EvidenceSource e(evidenceOptions, "PLACEHOLDER");
  e.loadFromFile(evidenceDf_mRNA, genes, pathway, sampleMap, sampleData);// ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM ALARM 
  evid.push_back(e);
  
  // /////////////////////////////////////////////////
  // Construct the factor graph
  vector< Factor > factors;
  vector< MaximizationStep > msteps;
  vector< vector < SharedParameters::FactorOrientations > > var_orders;

  //emSteps PropertySet
  EMSteps steps;
  string emSteps_plain = "[_mRNA.tab=-obs>]";
  PropertySet emSteps;
  emSteps.fromString(emSteps_plain);
  std::set<PropertyKey> keys = emSteps.keys();
	std::set<PropertyKey>::iterator i = keys.begin();
	EMStep es;
	for ( ; i != keys.end(); ++i) {
	  std::vector< std::string > edges;
	  edges=tokenizeString(emSteps.getAs<std::string>(*i),true,";");
	  SmallSet< std::string > s(edges.begin(), edges.end(), edges.size());
	  es[*i] = s;
	}
	steps.push_back(es);
  var_orders = pathway.constructFactors(steps, factors, msteps);
  map< long, string > outNodes = pathway.getOutputNodeMap(); //THIS ONE USES PATHWAY AND IS NEEDED, SOLVE THE input thing for pathway construction without read/write

  // add in additional factors to link disconnected pieces of the pathway
  FactorGraph *testGraphConnected = new FactorGraph(factors);

  while(!testGraphConnected->isConnected())
  {
    map< long, bool > nodesInSubNet = nodesInSingleNet(factors);
    long backNode = nodesInSubNet.rbegin()->first;
    VarSet I_vars;
    I_vars |= Var(backNode, PathwayTab::VARIABLE_DIMENSION);
    bool addedLink = false;
    // now go find ones that aren't connected to this
    vector< Factor >::iterator factorIter;
    for ( factorIter = factors.begin(); !addedLink && factorIter != factors.end(); ++factorIter) {
      const VarSet tmpVars = factorIter->vars();
      vector< Var >::const_iterator tmpVarsIter;
      for (tmpVarsIter = tmpVars.begin(); !addedLink && tmpVarsIter != tmpVars.end(); ++tmpVarsIter) {
	const long varLabel = tmpVarsIter->label();
	if(!nodesInSubNet[varLabel])
	  {
	    I_vars |= Var(varLabel, PathwayTab::VARIABLE_DIMENSION);

	    factors.push_back( Factor( I_vars, 1.0 ) );
	    addedLink = true;
	  }
      }
    }
    break;
    delete testGraphConnected;
    testGraphConnected = new FactorGraph(factors);
  }
  delete testGraphConnected;
  FactorGraph priorFG(factors);

  //PropertySet inferenceOptions = conf.getInferenceProperties(pathwayFilename);
  string inference = "[method=JTREE,updates=HUGIN,verbose=0]";
  PropertySet inferenceOptions;
  inferenceOptions.fromString(inference);
  std::string method = inferenceOptions.getAs<std::string>("method");
  InfAlg* prior = newInfAlg("JTREE", priorFG, inferenceOptions);
  prior->init();

  //Run EM
  string emset = "[max_iters=50,log_z_tol=0.01]";
  PropertySet emOptions;
  emOptions.fromString(emset);
  const PropertySet& em_conf = emOptions;
  Evidence evidence(sampleData);
  EMAlg em(evidence, *prior, msteps, em_conf);
  while(!em.hasSatisfiedTermConditions()) {
    em.iterate();
  }
  em.run();
  CharacterVector v = {"PLACEHOLDER"};
  NumericVector w = {0};
  // Creating DataFrame df
  DataFrame df = DataFrame::create( Named("V1") = v,         // simple assign
                                    Named("V2") = w); // using clone()
  prior->run();
  // Run inference on each of the samples
  map<string, size_t>::iterator sample_iter = sampleMap.begin();
  for ( ; sample_iter != sampleMap.end(); ++sample_iter) {
    InfAlg* clamped = prior->clone();
    Evidence::Observation *e = &sampleData[sample_iter->second];
    for (Evidence::Observation::const_iterator i = e->begin(); i != e->end(); ++i) {
      clamped->clamp( clamped->fg().findVar(i->first), i->second);
    }
    clamped->init();
    clamped->run();
    DataFrame dff = outputFastaPerturbations(df, sample_iter->first, prior, clamped, priorFG,
			     outNodes);
    df = dff;
    delete clamped;
  }

  delete prior;

  return df;
}
