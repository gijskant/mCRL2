
// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file quotienting_test.cpp
/// \brief Add your file description here.

#include <iostream>
#include <iterator>
#include <boost/test/included/unit_test_framework.hpp>
#include "mcrl2/utilities/test_utilities.h"
#include "mcrl2/lps/parse.h"
#include "mcrl2/modal_formula/parelm.h"
#include "mcrl2/modal_formula/parse.h"

using namespace mcrl2;

BOOST_AUTO_TEST_CASE(test_parelm)
{
  const std::string LPS =
    "act a;                                  \n"
    "    b : Nat;                            \n"
    "                                        \n"
    "proc P = P;                             \n"
    "                                        \n"
    "init P;                                 \n"
    ;

  const std::string FORMULA_TEXT =
    "nu X(c: Pos = 2, d: Nat = 1) . <a>X(c, d) && val(d > 0) && mu Y(e: Nat = d+1) . [b(e)](X(c, e) || Y(e))";

  const std::string EXPECTED_RESULT_TEXT =
    "nu X(d: Nat = 1) . <a>X(d) && val(d > 0) && mu Y(e: Nat = d+1) . [b(e)](X(e) || Y(e))";

  lps::specification spec = lps::parse_linear_process_specification(LPS);

  state_formulas::state_formula formula = state_formulas::parse_state_formula(FORMULA_TEXT, spec);
  state_formulas::state_formula expected_result = state_formulas::parse_state_formula(EXPECTED_RESULT_TEXT, spec);

  std::clog << "formula: " << pp(formula) << std::endl;

  std::clog << "eliminating unused parameters..." << std::endl;
  state_formulas::state_formula_parelm_algorithm algorithm;
  algorithm.run(formula);

  std::clog << "result: " << pp(formula) << std::endl;

  std::string result_text = pp(formula);
  BOOST_CHECK(formula == expected_result);
}

boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[])
{
  return 0;
}
