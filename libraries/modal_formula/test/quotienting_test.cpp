
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
#include "mcrl2/modal_formula/quotienting.h"
#include "mcrl2/modal_formula/parse.h"

using namespace mcrl2;

BOOST_AUTO_TEST_CASE(test_quotienting)
{
  const std::string LPS1 =
    "act a, b, c, d, e : Nat;                \n"
    "                                        \n"
    "proc P1(n:Nat) = true -> a(n) . P1(n);  \n"
    "                                        \n"
    "init P1(1);                             \n"
    ;

  const std::string LPS2 =
    "act a, b, c, d, e : Nat;                \n"
    "                                        \n"
    "proc P2(m:Nat) = true -> b(m) . P2(m+1) \n"
    "               + true -> e(m) . P2(1);  \n"
    "                                        \n"
    "init P2(1);                             \n"
    ;

  lps::specification spec1 = lps::parse_linear_process_specification(LPS1);
  lps::specification spec2 = lps::parse_linear_process_specification(LPS2);

  std::clog << "spec1: " << pp(spec1) << std::endl;
  std::clog << "spec2: " << pp(spec2) << std::endl;

  lps::synchronization_vector v(2);
  process::action_label_list labels = process::parse_action_declaration("c:Nat;");
  process::action_label label = labels.front();
  v.add_vector(std::vector<std::string> {"a", "b"}, label);
  v.add_vector(std::vector<std::string> {"d", "e"}, label);

  std::clog << "synchronization vector: " << pp(v) << std::endl;

  state_formulas::state_formula formula = state_formulas::parse_state_formula("nu X. exists j:Nat . <c(j)>X", spec1);

  std::clog << "formula: " << pp(formula) << std::endl;

  std::clog << "quotienting formula..." << std::endl;
  //state_formulas::state_formula q_formula = state_formulas::algorithms::quotient(formula, spec1, v, 0);
  state_formulas::quotient_builder f(spec1, v, 0);
  state_formulas::state_formula q_formula = f.apply(formula);
  process::label_generator label_generator = f.label_generator();
  for(auto it = label_generator.generated_labels().begin(); it != label_generator.generated_labels().end(); ++it)
  {
    spec1.action_labels().push_front(*it);
  }

  std::clog << "result: " << pp(q_formula) << std::endl;

  std::string q_formula_text = pp(q_formula);
  state_formulas::state_formula formula2 = state_formulas::parse_state_formula(q_formula_text, spec1);
  BOOST_CHECK(formula2 == q_formula);

  BOOST_CHECK(true);
}

boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[])
{
  return 0;
}
