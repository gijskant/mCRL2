
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
#include "mcrl2/modal_formula/rewrite_quantifiers.h"
#include "mcrl2/modal_formula/parse.h"

using namespace mcrl2;

BOOST_AUTO_TEST_CASE(test_rewrite_quantifiers)
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
    "nu X . forall e: Nat . val(e == 3) => val(e > 1)";

  lps::specification spec = lps::parse_linear_process_specification(LPS);

  state_formulas::state_formula formula = state_formulas::parse_state_formula(FORMULA_TEXT, spec);

  std::clog << "formula: " << pp(formula) << std::endl;

  std::clog << "applying one point rule..." << std::endl;
  state_formulas::state_formula result = state_formulas::rewrite_quantifiers(formula, spec.data());

  std::clog << "result: " << pp(result) << std::endl;
}

boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[])
{
  log::mcrl2_logger::set_reporting_level(log::debug);
  return 0;
}
