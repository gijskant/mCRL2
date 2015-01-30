// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_LTS_DETAIL_NETWORK_EXPLORE_H
#define MCRL2_LTS_DETAIL_NETWORK_EXPLORE_H

#include <string>
#include <limits>
#include <memory>

#include "mcrl2/atermpp/indexed_set.h"
#include "mcrl2/atermpp/aterm_balanced_tree.h"
#include "mcrl2/lps/network.h"
#include "mcrl2/lps/next_state_generator.h"
#include "mcrl2/lts/lts_lts.h"
#include "mcrl2/lts/detail/bithashtable.h"
#include "mcrl2/lts/detail/queue.h"
#include "mcrl2/lts/detail/exploration_strategy.h"

namespace mcrl2
{

namespace lts
{

struct network_explore_options
{
  network_explore_options() :
    strat(mcrl2::data::jitty),
    expl_strat(es_breadth),
    suppress_progress_messages(false),
    outformat(mcrl2::lts::lts_none),
    outinfo(true),
    detect_deadlock(false)
  {  }

  lps::network network;
  std::vector<lps::stochastic_specification> specifications;

  mcrl2::data::rewriter::strategy strat;
  exploration_strategy expl_strat;
  bool suppress_progress_messages;

  mcrl2::lts::lts_type outformat;
  bool outinfo;
  std::string lts;

  bool detect_deadlock;

};


class network2lts_algorithm
{
  private:
  typedef lps::next_state_generator next_state_generator;
  typedef lps::state state_t;

  private:
    network_explore_options m_options;
    std::vector<next_state_generator*> m_generators;
    /*next_state_generator::substitution_t m_substitution;*/

    atermpp::indexed_set<state_t> m_state_numbers;
    bit_hash_table m_bit_hash_table;

    lts_lts_t m_output_lts;
    atermpp::indexed_set<atermpp::aterm_appl> m_action_label_numbers;
    std::ofstream m_aut_file;

    size_t m_num_states;
    size_t m_num_transitions;
    size_t m_initial_state_number;
    size_t m_level;

    volatile bool m_must_abort;

  public:
    network2lts_algorithm() :
      m_num_states(0),
      m_num_transitions(0),
      m_initial_state_number(0),
      m_level(0),
      m_must_abort(false)
    {
    }

    virtual ~network2lts_algorithm()
    {
      for(auto it=m_generators.begin(); it != m_generators.end(); ++it)
      {
        delete (*it);
      }
    }

    virtual bool initialise_lts_generation(network_explore_options* options);
    virtual bool generate_lts();
    virtual bool finalise_lts_generation();

    virtual void abort()
    {
      // Stops the exploration algorithm if it is running by making sure
      // not a single state can be generated anymore.
      if (!m_must_abort)
      {
        m_must_abort = true;
        mCRL2log(log::warning) << "state space generation was aborted prematurely" << std::endl;
      }
    }

  private:
    void save_deadlock(const state_t& state);
    bool add_transition(const state_t& state, next_state_generator::transition_t &transition);
    void get_transitions(const state_t &state,
                         std::list<next_state_generator::transition_t> &transitions,
                         next_state_generator::enumerator_queue_t& enumeration_queue
                         );

    void generate_lts_breadth();
    void generate_lts_depth(const state_t &initial_state);
};

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LTS_DETAIL_NETWORK_EXPLORE_H
