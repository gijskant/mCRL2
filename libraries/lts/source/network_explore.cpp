// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <time.h>

//#include "mcrl2/atermpp/aterm_balanced_tree.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/network_next_state_generator.h"
#include "mcrl2/lps/resolve_name_clashes.h"
#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lts/detail/network_explore.h"
#include "mcrl2/lts/lts_io.h"
#include "mcrl2/trace/trace.h"
#include <iomanip>

using namespace mcrl2;
using namespace mcrl2::log;
using namespace mcrl2::lps;
using namespace mcrl2::lts;

bool network2lts_algorithm::initialise_lts_generation(network_explore_options *options)
{
  m_options=*options;
  m_state_numbers = atermpp::indexed_set<state_t>();
  m_num_states = 0;
  m_num_transitions = 0;
  m_level = 1;

  lps::network network(m_options.network);
  std::vector<data::variable> parameters;
  for(auto s_it = network.lps_filenames().begin(); s_it != network.lps_filenames().end(); ++s_it)
  {
    lps::stochastic_specification spec;
    load_lps(spec, *s_it);
    resolve_summand_variable_name_clashes(spec);
    m_options.specifications.push_back(spec);
    parameters.insert(parameters.end(), spec.process().process_parameters().begin(), spec.process().process_parameters().end());

    data::rewriter rewriter;
    rewriter = data::rewriter(spec.data(), m_options.strat);

    m_generators.push_back(new next_state_generator(spec, rewriter));
  }

  if (m_options.outformat == lts_aut)
  {
    mCRL2log(verbose) << "writing state space in AUT format to '" << m_options.lts << "'." << std::endl;
    m_aut_file.open(m_options.lts.c_str());
    if (!m_aut_file.is_open())
    {
      mCRL2log(error) << "cannot open '" << m_options.lts << "' for writing" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if (m_options.outformat == lts_none)
  {
    mCRL2log(verbose) << "not saving state space." << std::endl;
  }
  else
  {
    mCRL2log(verbose) << "writing state space in " << mcrl2::lts::detail::string_for_type(m_options.outformat)
                      << " format to '" << m_options.lts << "'." << std::endl;
    m_output_lts.set_data(m_options.specifications[0].data());
    m_output_lts.set_process_parameters(data::variable_list(parameters.begin(), parameters.end()));
    m_output_lts.set_action_label_declarations(network.action_labels());
  }

  if (m_options.detect_deadlock)
  {
    mCRL2log(verbose) << "Detect deadlocks.\n" ;
  }

  return true;
}

bool network2lts_algorithm::generate_lts()
{
  std::vector<data::data_expression> initial_state_vector;
  for(auto g = m_generators.begin(); g != m_generators.end(); ++g)
  {
    auto initial_states = (*g)->initial_states();
    if (initial_states.empty())
    {
      mCRL2log(error) << "list of initial states is empty. " << std::endl;
    }
    state_t local_state = initial_states.front().state();
    initial_state_vector.insert(initial_state_vector.end(), local_state.begin(), local_state.end());
  }
  state_t initial_state(initial_state_vector.begin(), initial_state_vector.size());
  m_initial_state_number=0;

  m_state_numbers.put(initial_state);

  if (m_options.outformat == lts_aut)
  {
    // HACK: this line will be overwritten once generation is finished.
    m_aut_file << "                                                             " << std::endl;
  }
  else if (m_options.outformat != lts_none)
  {
    m_initial_state_number = m_output_lts.add_state(state_label_lts(initial_state));
    m_output_lts.set_initial_state(m_initial_state_number);
  }
  m_num_states = 1;

  mCRL2log(verbose) << "generating state space with '" << m_options.expl_strat << "' strategy...\n";

  if (m_options.expl_strat == es_breadth || m_options.expl_strat == es_value_prioritized)
  {
    generate_lts_breadth();

    mCRL2log(verbose) << "done with state space generation ("
                      << m_level-1 << " level" << ((m_level==2)?"":"s") << ", "
                      << m_num_states << " state" << ((m_num_states == 1)?"":"s")
                      << " and " << m_num_transitions << " transition" << ((m_num_transitions==1)?"":"s") << ")" << std::endl;
  }
  else if (m_options.expl_strat == es_depth)
  {
    generate_lts_depth(initial_state);

    mCRL2log(verbose) << "done with state space generation ("
                      << m_num_states << " state" << ((m_num_states == 1)?"":"s")
                      << " and " << m_num_transitions << " transition" << ((m_num_transitions==1)?"":"s") << ")" << std::endl;
  }
  else if (m_options.expl_strat == es_random || m_options.expl_strat == es_value_random_prioritized)
  {
    mCRL2log(error) << "exploration strategy not implemented" << std::endl;
  }
  else
  {
    mCRL2log(error) << "unknown exploration strategy" << std::endl;
    return false;
  }

  return true;
}

bool network2lts_algorithm::finalise_lts_generation()
{
  if (m_options.outformat == lts_aut)
  {
    m_aut_file.flush();
    m_aut_file.seekp(0);
    m_aut_file << "des (" << m_initial_state_number << "," << m_num_transitions << "," << m_num_states << ")";
    m_aut_file.close();
  }
  else if (m_options.outformat != lts_none)
  {
    if (!m_options.outinfo)
    {
      m_output_lts.clear_state_labels();
    }

    switch (m_options.outformat)
    {
      case lts_lts:
      {
        m_output_lts.save(m_options.lts);
        break;
      }
      case lts_fsm:
      {
        lts_fsm_t fsm;
        detail::lts_convert(m_output_lts, fsm);
        fsm.save(m_options.lts);
        break;
      }
#ifdef USE_BCG
      case lts_bcg:
      {
        lts_bcg_t bcg;
        detail::lts_convert(m_output_lts, bcg);
        bcg.save(m_options.lts);
        break;
      }
#endif
      case lts_dot:
      {
        lts_dot_t dot;
        detail::lts_convert(m_output_lts, dot);
        dot.save(m_options.lts);
        break;
      }
      default:
        assert(0);
    }
  }

  return true;
}

void network2lts_algorithm::save_deadlock(const network2lts_algorithm::state_t &state)
{
  size_t state_number = m_state_numbers.index(state);
  mCRL2log(info) << "deadlock-detect: deadlock found (state index: " << state_number <<  ").\n";
}

bool network2lts_algorithm::add_transition(const network2lts_algorithm::state_t &source, next_state_generator::transition_t &transition)
{
  state_t destination = transition.target_state();

  size_t source_state_number;
  std::pair<size_t, bool> destination_state_number;

  source_state_number = m_state_numbers[source];
  destination_state_number = m_state_numbers.put(destination);

  if (destination_state_number.second)
  {
    m_num_states++;

    if (m_options.outformat != lts_none && m_options.outformat != lts_aut)
    {
      size_t state_number = m_output_lts.add_state(state_label_lts(transition.target_state()));
      assert(state_number == destination_state_number.first);
      static_cast <void>(state_number);
    }
  }

  if (m_options.outformat == lts_aut)
  {
    m_aut_file << "(" << source_state_number << ",\"" << lps::pp(transition.action()) << "\"," << destination_state_number.first << ")" << std::endl;
  }
  else if (m_options.outformat != lts_none)
  {
    std::pair<size_t, bool> action_label_number = m_action_label_numbers.put(lps::detail::multi_action_to_aterm(transition.action()));
    if (action_label_number.second)
    {
      size_t action_number = m_output_lts.add_action(action_label_lts(transition.action()));
      assert(action_number == action_label_number.first);
      static_cast <void>(action_number); // Avoid a warning when compiling in non debug mode.
    }

    m_output_lts.add_transition(mcrl2::lts::transition(source_state_number, action_label_number.first, destination_state_number.first));
  }

  m_num_transitions++;

  return destination_state_number.second;
}

// std::list<lps2lts_algorithm::next_state_generator::transition_t> lps2lts_algorithm::get_transitions(const lps2lts_algorithm::generator_state_t &state)
void network2lts_algorithm::get_transitions(const network2lts_algorithm::state_t &state,
                                        std::list<network2lts_algorithm::next_state_generator::transition_t> &transitions,
                                        next_state_generator::enumerator_queue_t& enumeration_queue
                                        )
{
  assert(transitions.empty());

  // std::list<next_state_generator::transition_t> transitions;
  try
  {
    lps::explore(m_options.specifications,
        m_options.network.synchronization_vector(),
        m_generators,
        /* &m_substitution,*/
        state,
        transitions,
        &enumeration_queue
    );
  }
  catch (mcrl2::runtime_error& e)
  {
    mCRL2log(error) << "Error while exploring state space: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }

  if (transitions.empty() && m_options.detect_deadlock)
  {
    save_deadlock(state);
  }

  // return transitions;
}

void network2lts_algorithm::generate_lts_breadth()
{
  size_t current_state = 0;
  size_t start_level_seen = 1;
  size_t start_level_transitions = 0;
  time_t last_log_time = time(NULL) - 1, new_log_time;
  next_state_generator::enumerator_queue_t enumeration_queue;

  while (!m_must_abort && (current_state < m_state_numbers.size()))
  {
    state_t state(m_state_numbers.get(current_state));
    std::list<next_state_generator::transition_t> transitions;
    get_transitions(state,transitions, enumeration_queue);

    for (std::list<next_state_generator::transition_t>::iterator i = transitions.begin(); i != transitions.end(); i++)
    {
      add_transition(state, *i);
    }

    current_state++;
    if (current_state == start_level_seen)
    {
      m_level++;
      start_level_seen = m_num_states;
      start_level_transitions = m_num_transitions;
    }

    if (!m_options.suppress_progress_messages && time(&new_log_time) > last_log_time)
    {
      last_log_time = new_log_time;
      size_t lvl_states = m_num_states - start_level_seen;
      size_t lvl_transitions = m_num_transitions - start_level_transitions;
      mCRL2log(status) << std::fixed << std::setprecision(2)
                       << m_num_states << "st, " << m_num_transitions << "tr"
                       << ", explored " << 100.0 * ((float)current_state / m_num_states)
                       << "%. Last level: " << m_level << ", " << lvl_states << "st, " << lvl_transitions
                       << "tr.\n";
    }
  }
}


void network2lts_algorithm::generate_lts_depth(const state_t &initial_state)
{
  std::list<state_t> stack;
  stack.push_back(initial_state);

  size_t current_state = 0;
  next_state_generator::enumerator_queue_t enumeration_queue;

  while (!m_must_abort && (!stack.empty()))
  {
    state_t state = stack.back();
    stack.pop_back();
    std::list<next_state_generator::transition_t> transitions;
    get_transitions(state,transitions, enumeration_queue);

    for (std::list<next_state_generator::transition_t>::iterator i = transitions.begin(); i != transitions.end(); i++)
    {
      if (add_transition(state, *i))
      {
        stack.push_back(i->target_state());
      }
    }

    current_state++;
    if (!m_options.suppress_progress_messages && ((current_state % 1000) == 0))
    {
      mCRL2log(verbose) << "monitor: currently explored "
                        << current_state << " state" << ((current_state==1)?"":"s")
                        << " and " << m_num_transitions << " transition" << ((m_num_transitions==1)?"":"s")
                        << " (stacksize is " << stack.size() << ")" << std::endl;
    }
  }
}

