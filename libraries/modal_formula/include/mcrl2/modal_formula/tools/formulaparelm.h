// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/tools/formulaparelm.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_TOOLS_FORMULAPARELM_H
#define MCRL2_MODAL_FORMULAS_TOOLS_FORMULAPARELM_H

#include <algorithm>
#include <fstream>
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/network.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/synchronization_vector.h"
#include "mcrl2/modal_formula/algorithms.h"
#include "mcrl2/modal_formula/parelm.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/utilities/logger.h"

namespace mcrl2 {

namespace state_formulas {

void formulaparelm(const std::string& input_filename,
              const std::string& output_filename,
              const std::string& lps_filename,
              const std::string& network_filename
             )
{
  // load formula file
  state_formulas::state_formula formula;
  lps::specification lps; // LPS containing the data specification used for parsing the formula.
  lps::synchronization_vector v;
  lps::network network;
  std::set<process::action_label> action_labels;
  if (!network_filename.empty())
  {
    // read network of LPSs from file
    mCRL2log(log::verbose) << "reading network from file '" <<  network_filename << "'..." << std::endl;
    network.load(network_filename);
    size_t i = 0;
    if (i < network.lps_filenames().size())
    {
      mCRL2log(log::verbose) << "reading LPS at index " << i << " from file '"
          <<  network.lps_filenames()[i] << "'..." << std::endl;
      load_lps(lps, network.lps_filenames()[i]);
    }
    else
    {
      throw std::runtime_error("Index should be in the range [0 .. n-1], where n is the size of the network.");
    }
    network.synchronization_vector().normalize(lps.data());
    v = network.synchronization_vector();
    action_labels = v.action_labels();
  }
  else
  {
    if (lps_filename.empty())
    {
      throw std::runtime_error("options -n and -p not specified. An LPS is required for parsing the formula.");
    }
    load_lps(lps, lps_filename);
  }
  std::set<process::action_label> lps_labels(lps.action_labels().begin(), lps.action_labels().end());
  std::set<process::action_label> vector_labels;
  std::set_difference(action_labels.begin(), action_labels.end(), lps_labels.begin(), lps_labels.end(),
                        std::inserter(vector_labels, vector_labels.begin()));
  // add action labels
  for(auto it = vector_labels.begin(); it != vector_labels.end(); ++it)
  {
    //mCRL2log(log::debug) << "adding action label: " << pp(*it) << ": " << pp((*it).sorts()) << "." << std::endl;
    lps.action_labels().push_front(*it);
  }

  // load formula file
  if (input_filename.empty())
  {
    mCRL2log(log::verbose) << "reading formula from stdin..." << std::endl;
    formula = state_formulas::algorithms::parse_state_formula(std::cin, lps);
  }
  else
  {
    mCRL2log(log::verbose) << "reading formula from file '" <<  input_filename << "'..." << std::endl;
    std::ifstream instream(input_filename.c_str(), std::ifstream::in);
    if (!instream)
    {
      throw mcrl2::runtime_error("cannot open input file: " + input_filename);
    }
    formula = state_formulas::algorithms::parse_state_formula(instream, lps);
    instream.close();
  }

  // applying parameter elimination
  mCRL2log(log::verbose) << "removing unused parameters..." << std::endl;
  state_formulas::state_formula_parelm_algorithm algorithm;
  algorithm.run(formula);

  // save the result formula
  if (output_filename.empty())
  {
    mCRL2log(log::verbose) << "writing result formula to stdout..." << std::endl;
  }
  else
  {
    mCRL2log(log::verbose) << "writing result formula to file '" <<  output_filename << "'..." << std::endl;
  }
  if (output_filename.empty())
  {
    std::cout << pp(formula) << std::endl;
  }
  else
  {
    std::ofstream outstream(output_filename.c_str(), std::ofstream::out);
    if (!outstream)
    {
      throw mcrl2::runtime_error("cannot open output file: " + output_filename);
    }
    outstream << pp(formula) << std::endl;
    outstream.close();
  }
  mCRL2log(log::verbose) << "done." << std::endl;
}

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULAS_TOOLS_FORMULAQUOTIENT_H
