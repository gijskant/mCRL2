// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file next_state_generator.h

#ifndef MCRL2_LPS_NETWORK_NEXT_STATE_GENERATOR_H
#define MCRL2_LPS_NETWORK_NEXT_STATE_GENERATOR_H

#include <iterator>
#include <string>
#include <vector>

#include "mcrl2/lps/network.h"
#include "mcrl2/lps/next_state_generator.h"
#include "mcrl2/lps/synchronization_vector.h"

namespace mcrl2
{

namespace lps
{

  typedef lps::state state_t;
  typedef next_state_generator::transition_t transition_t;
  /*typedef next_state_generator::substitution_t substitution_t;*/


  std::string state_to_string(next_state_generator* generator, const state_t& state)
  {
    std::stringstream ss;
    ss << "<";
    bool first = true;
    for(auto it = state.begin(); it != state.end(); ++it)
    {
      if (!first)
      {
        ss << ", ";
      }
      ss << data::pp(*it);
      first = false;
    }
    ss << ">";
    return ss.str();
  }

  void explore(const std::vector<lps::stochastic_specification>& specifications,
      const lps::synchronization_vector& synchronization_vector,
      const std::vector<next_state_generator*>& generators,
      /*substitution_t* substitution,*/
      const state_t& state,
      std::list<transition_t> &transitions,
      next_state_generator::enumerator_queue_t* enumeration_queue
      )
  {
    mCRL2log(log::debug) << "Explore global state: " << state_to_string(generators[0], state) <<
        " (" << state.size() << ")" << std::endl;
    // split vector in parts and compute local successors
    std::vector<state_t> local_source;
    std::vector<std::list<transition_t>> local_transitions;
    auto state_it = state.begin();
    size_t index = 0;
    for(auto s_it = specifications.begin(); s_it != specifications.end(); ++s_it)
    {
      mCRL2log(log::debug) << "Exploring local transitions of component " << index << "..." << std::endl;
      lps::stochastic_specification spec = (*s_it);
      size_t size = spec.process().process_parameters().size();
      std::vector<data::data_expression> local_state_vector;
      if (state_it == state.end())
      {
        throw std::runtime_error("Incorrect state length.");
      }
      for(size_t i=0; i < size; i++)
      {
        local_state_vector.push_back(*(state_it++));
      }
      state_t local_state(local_state_vector.begin(), local_state_vector.size());
      mCRL2log(log::debug) << " - local state: " << state_to_string(generators[index], local_state) << std::endl;
      local_source.push_back(local_state);
      next_state_generator::iterator it(generators[index]->begin(local_state, enumeration_queue /*, substitution*/));
      std::list<transition_t> local_trans;
      while (it)
      {
        transition_t trans = (*it++);
        mCRL2log(log::debug) << " - local transitions: " << state_to_string(generators[index], trans.target_state()) << std::endl;
        local_trans.push_back(trans);
      }
      local_transitions.push_back(local_trans);
      index++;
    }
    // Considering synchronization vector entries
    auto vector = synchronization_vector.vector();
    size_t vector_index = 0;
    for(auto it = vector.begin(); it != vector.end(); ++it)
    {
      auto entry = (*it);
      mCRL2log(log::debug) << "Considering synchronization vector entry: " << pp(entry) << std::endl;
      std::vector<std::string> names = entry.first;
      process::action_label label = entry.second;
      std::set<data::data_expression_list> arguments; // set of possible synchronization values
      for(size_t i=0; i < names.size(); i++)
      {
        std::string name = names[i];
        if (name != lps::inactive_label)
        {
          std::list<transition_t> local_trans = local_transitions[i];
          for(auto t_ = local_trans.begin(); t_ != local_trans.end(); ++t_)
          {
            transition_t t = (*t_);
            lps::multi_action m = t.action();
            if (m.actions().size() == 0) // tau
            {

            }
            else if (m.actions().size() == 1)
            {
              process::action action = m.actions().front();
              if (action.label().name() == core::identifier_string(name))
              {
                arguments.insert(action.arguments());
              }
            }
            else
            {
              throw std::runtime_error("Unexpected multi-action.");
            }
          }
        }

      }
      // Considering synchronization values
      for(auto a_ = arguments.begin(); a_ != arguments.end(); ++a_)
      {
        mCRL2log(log::debug) << "Considering synchronizing on (" << pp(*a_) << ")" << std::endl;
        std::set<std::vector<state_t> > local_targets;
        std::set<std::vector<state_t> > local_targets_next;
        local_targets.insert(std::vector<state_t>());
        bool valid_transition = true;
        data::data_expression_list args = (*a_);
        for(size_t i=0; i < names.size(); i++)
        {
          bool valid = false;
          std::string name = names[i];
          if (name == lps::inactive_label)
          {
            for(auto tgt_it = local_targets.begin(); tgt_it != local_targets.end(); ++tgt_it)
            {
              std::vector<state_t> local_target = *tgt_it;
              local_target.push_back(local_source[i]);
              local_targets_next.insert(local_target);
            }
            valid = true;
          }
          else
          {
            std::list<transition_t> local_trans = local_transitions[i];
            for(auto t_ = local_trans.begin(); t_ != local_trans.end(); ++t_)
            {
              transition_t t = (*t_);
              lps::multi_action m = t.action();
              if (m.actions().size() == 0) // tau
              {

              }
              else if (m.actions().size() == 1)
              {
                process::action action = m.actions().front();
                if (action.label().name() == core::identifier_string(name))
                {
                  if(action.arguments() == *a_)
                  {
                    for(auto tgt_it = local_targets.begin(); tgt_it != local_targets.end(); ++tgt_it)
                    {
                      std::vector<state_t> local_target = *tgt_it;
                      local_target.push_back(t.target_state());
                      local_targets_next.insert(local_target);
                    }
                    valid = true;
                  }
                }
              }
              else
              {
                throw std::runtime_error("Unexpected multi-action.");
              }

            }
          }
          local_targets.clear();
          local_targets = local_targets_next;
          local_targets_next = std::set<std::vector<state_t> >();
          valid_transition &= valid;
        }
        // Synchronization on a_ is possible, construct transition
        if (valid_transition)
        {
          mCRL2log(log::debug) << "Synchronizing on (" << pp(*a_) << ")" << std::endl;
          // construct target state and transition
          for(auto targets_it = local_targets.begin(); targets_it != local_targets.end(); ++ targets_it)
          {
            std::vector<state_t> local_target = *targets_it;
            std::vector<data::data_expression> target_vector;
            for(auto target_it = local_target.begin(); target_it != local_target.end(); ++target_it)
            {
              target_vector.insert(target_vector.end(), target_it->begin(), target_it->end());
            }
            state_t target(target_vector.begin(), target_vector.size());
            transition_t transition;
            transition.set_target_state(target);
            transition.set_action(multi_action(process::action(label, args)));
            transition.set_summand_index(vector_index);
            transitions.push_back(transition);
            mCRL2log(log::debug) << "Target: " << state_to_string(generators[0], target) <<
                " (" << target.size() << ")" << std::endl;
          }
        }
      }
      vector_index++;
    }
  }

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_NETWORK_NEXT_STATE_GENERATOR_H
