// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/quotienting.h
/// \brief The class quotienting.

#ifndef MCRL2_LPS_QUOTIENTING_H_
#define MCRL2_LPS_QUOTIENTING_H_

#include "boost/static_assert.hpp"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/network.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/synchronization_vector.h"
#include "mcrl2/process/label_generator.h"

namespace mcrl2 {

namespace lps {

//FIXME: bug
/// \brief Quotients the synchronization vector for component spec at index i.
/// \param v the synchronization vector
/// \param i the index of the component that is quotiented out
/// \param spec the component that is quotiented out
lps::synchronization_vector quotient( const lps::synchronization_vector& v,
                                      size_t i,
                                      const lps::specification& spec,
                                      process::label_generator& label_generator)
{
  if (v.length() == 0)
  {
    throw std::runtime_error("Length of synchronization vector is 0. Quotienting not possible.");
  }
  lps::synchronization_vector result(v.length() - 1);
  for(auto v_it = v.vector().begin(); v_it != v.vector().end(); ++v_it)
  {
    std::vector<std::string> labels = (*v_it).first;
    std::string a_i = labels[i];
    process::action_label a((*v_it).second);
    bool other_process_active = false;
    // construct new vector
    std::vector<std::string> l;
    size_t j = 0;
    for(auto l_it = labels.begin(); l_it != labels.end(); ++l_it)
    {
      if (i != j)
      {
        if ((*l_it) != lps::inactive_label)
        {
          other_process_active = true;
        }
        l.push_back(*l_it);
      }
      j++;
    }
    if (a_i == lps::inactive_label) {
      result.add_vector(l, a);
    }
    else
    {
      if (other_process_active)
      {
        // case P_i synchronises with other processes
        process::action_label rho_label = label_generator.fresh_action_label(labels, a);
        result.add_vector(l, rho_label);
      }
      else
      {
        // case P_i is the only active process
      }
    }
  }
  return result;
}

lps::network quotient(const lps::network& n,
                      size_t i,
                      process::label_generator& label_generator)
{
  std::vector<std::string> lps_filenames;
  size_t j = 0;
  lps::specification spec;
  for(auto it = n.lps_filenames().begin(); it != n.lps_filenames().end(); ++it)
  {
    if (i == j)
    {
      lps::load_lps(spec, *it);
    }
    else
    {
      lps_filenames.push_back(*it);
    }
    j++;
  }
  lps::synchronization_vector synchronization_vector = quotient(n.synchronization_vector(), i, spec, label_generator);
  process::action_label_list action_labels(label_generator.generated_labels().begin(), label_generator.generated_labels().end());
  lps::network result(lps_filenames, synchronization_vector, action_labels);
  return result;
}

} // namespace lps

} // namespace mcrl2


#endif /* MCRL2_LPS_QUOTIENTING_H_ */
