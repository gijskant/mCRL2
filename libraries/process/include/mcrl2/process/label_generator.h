// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/label_generator.h
/// \brief add your file description here.

#ifndef LIBRARIES_PROCESS_INCLUDE_MCRL2_PROCESS_LABEL_GENERATOR_H_
#define LIBRARIES_PROCESS_INCLUDE_MCRL2_PROCESS_LABEL_GENERATOR_H_

namespace mcrl2 {

namespace process {

struct label_generator
{
  data::set_identifier_generator id_generator;
  std::map<std::pair<std::vector<std::string>, process::action_label>, process::action_label > generated_labels_map;
  std::vector<process::action_label> generated_labels_vector;

  label_generator()
  {
  }

  ~label_generator()
  {
  }

  const std::vector<process::action_label>& generated_labels() const
  {
    return generated_labels_vector;
  }

  std::vector<process::action_label>& generated_labels()
  {
    return generated_labels_vector;
  }

  void add_identifier(core::identifier_string label_name)
  {
    id_generator.add_identifier(label_name);
  }

  process::action_label fresh_action_label(
      std::vector<std::string> labels,
      process::action_label l)
  {
    assert(!labels.empty());
    auto key = std::make_pair(labels, l);
    auto it = generated_labels_map.find(key);
    if (it != generated_labels_map.end())
    {
      return (*it).second;
    }
    else
    {
      core::identifier_string label_name;
      std::string l_name(l.name());
      // special case for tau actions, rename them to t:
      if (l_name.size() >= 3 && l_name.substr(0,3) == "tau")
      {
        core::identifier_string t("t1");
        label_name = id_generator(t);
      }
      else
      {
        label_name = id_generator(l.name());
      }
      process::action_label label(label_name, l.sorts());
      generated_labels_map[key] = label;
      generated_labels_vector.push_back(label);
      return label;
    }
  }

};

} // namespace process

} // namespace mcrl2

#endif /* LIBRARIES_PROCESS_INCLUDE_MCRL2_PROCESS_LABEL_GENERATOR_H_ */
