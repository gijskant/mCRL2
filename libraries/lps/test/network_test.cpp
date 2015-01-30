
// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file network_test.cpp
/// \brief Add your file description here.

#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include <boost/test/included/unit_test_framework.hpp>
#include "mcrl2/utilities/test_utilities.h"
#include "mcrl2/lps/network.h"
#include "mcrl2/lps/parse.h"

using namespace mcrl2;

BOOST_AUTO_TEST_CASE(test_synchronization_vector)
{
  lps::synchronization_vector v(2);
  process::action_label_list labels = process::parse_action_declaration("c:Nat;");
  process::action_label label = labels.front();
  v.add_vector(std::vector<std::string> {"a", "b"}, label);
  v.add_vector(std::vector<std::string> {"d", "e"}, label);

  std::clog << "synchronization vector: " << pp(v) << std::endl;

  std::stringstream s_out;
  v.write(s_out);
  std::string s = s_out.str();
  std::clog << "result of write: " << s << std::endl;

  std::stringstream s_in(s);
  lps::synchronization_vector w;
  w.read(s_in);

  std::clog << "result of read: " << pp(w) << std::endl;

  BOOST_CHECK(pp(v) == pp(w));
}

boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[])
{
  return 0;
}
