// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file network2lts.cpp

#include <string>
#include <cassert>
#include <signal.h>

#include "boost/lexical_cast.hpp"

#include "mcrl2/utilities/logger.h"

#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/data/rewriter_tool.h"

#include "mcrl2/lps/multi_action.h"
#include "mcrl2/process/action_parse.h"

#include "mcrl2/lts/lts_io.h"
#include "mcrl2/lts/detail/network_explore.h"

#define __STRINGIFY(x) #x
#define STRINGIFY(x) __STRINGIFY(x)

using namespace std;
using namespace mcrl2::data::tools;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::utilities;
using namespace mcrl2::core;
using namespace mcrl2::lts;
using namespace mcrl2::lps;
using namespace mcrl2::log;

typedef mcrl2::data::tools::rewriter_tool< input_output_tool > network2lts_base;
class network2lts_tool : public network2lts_base
{
  protected:
    mcrl2::lts::network2lts_algorithm m_network2lts;
    network_explore_options m_options;
    std::string m_filename;

  public:
    network2lts_tool() :
      network2lts_base("network2lts","Gijs Kant",
                   "generate an LTS from a network of LPSs",
                   "Generate an LTS from the network in INFILE and save the result to OUTFILE. "
                   "If INFILE is not supplied, stdin is used. "
                   "If OUTFILE is not supplied, the LTS is not stored.\n"
                   "\n"
                   "If the 'jittyc' rewriter is used, then the MCRL2_COMPILEREWRITER environment "
                   "variable (default value: 'mcrl2compilerewriter') determines the script that "
                   "compiles the rewriter, and MCRL2_COMPILEDIR (default value: '.') determines "
                   "where temporary files are stored.\n"
                   "\n"
                   "Note that network2lts can deliver multiple transitions with the same label between"
                   "any pair of states. If this is not desired, such transitions can be removed by"
                   "applying a strong bisimulation reducton using for instance the tool ltsconvert.\n"
                   "\n"
                   "The format of OUTFILE is determined by its extension (unless it is specified "
                   "by an option). The supported formats are:\n"
                   "\n"
                   +mcrl2::lts::detail::supported_lts_formats_text()+"\n"
                   "If the jittyc rewriter is used, then the MCRL2_COMPILEREWRITER environment "
                   "variable (default value: mcrl2compilerewriter) determines the script that "
                   "compiles the rewriter, and MCRL2_COMPILEDIR (default value: '.') "
                   "determines where temporary files are stored."
                   "\n"
                   "Note that network2lts can deliver multiple transitions with the same "
                   "label between any pair of states. If this is not desired, such "
                   "transitions can be removed by applying a strong bisimulation reducton "
                   "using for instance the tool ltsconvert."
                  )
    {
    }

    void abort()
    {
      m_network2lts.abort();
    }

    bool run()
    {
      m_options.network.load(m_filename);

      if (!m_network2lts.initialise_lts_generation(&m_options))
      {
        return false;
      }

      try
      {
        m_network2lts.generate_lts();
      }
      catch (mcrl2::runtime_error& e)
      {
        mCRL2log(error) << e.what() << std::endl;
        m_network2lts.finalise_lts_generation();
        return false;
      }

      m_network2lts.finalise_lts_generation();

      return true;
    }

  protected:
    void add_options(interface_description& desc)
    {
      network2lts_base::add_options(desc);

      desc.
      add_option("deadlock",
                 "detect deadlocks (i.e. for every deadlock a message is printed)", 'D').
      add_option("strategy", make_enum_argument<exploration_strategy>("NAME")
                 .add_value_short(es_breadth, "b", true)
                 .add_value_short(es_depth, "d")
                 .add_value_short(es_value_prioritized, "p")
                 .add_value_short(es_value_random_prioritized, "q")
                 .add_value_short(es_random, "r")
                 , "explore the state space using strategy NAME:"
                 , 's').
      add_option("out", make_mandatory_argument("FORMAT"),
                 "save the output in the specified FORMAT", 'o').
      add_option("no-info", "do not add state information to OUTFILE"
                 "Without this option lps2lts adds state vector to the LTS. This "
                 "option causes this information to be discarded and states are only "
                 "indicated by a sequence number. Explicit state information is useful "
                 "for visualisation purposes, for instance, but can cause the OUTFILE "
                 "to grow considerably. Note that this option is implicit when writing "
                 "in the AUT format.").
      add_option("suppress","in verbose mode, do not print progress messages indicating the number of visited states and transitions. "
                 "For large state spaces the number of progress messages can be quite "
                 "horrendous. This feature helps to suppress those. Other verbose messages, "
                 "such as the total number of states explored, just remain visible.");
    }

    void parse_options(const command_line_parser& parser)
    {
      network2lts_base::parse_options(parser);
      m_options.detect_deadlock = parser.options.count("deadlock") != 0;
      m_options.outinfo         = parser.options.count("no-info") == 0;
      m_options.suppress_progress_messages = parser.options.count("suppress") !=0;
      m_options.strat           = parser.option_argument_as< mcrl2::data::rewriter::strategy >("rewriter");

      m_options.expl_strat = parser.option_argument_as<exploration_strategy>("strategy");

      if (parser.options.count("out"))
      {
        m_options.outformat = mcrl2::lts::detail::parse_format(parser.option_argument("out"));

        if (m_options.outformat == lts_none)
        {
          parser.error("format '" + parser.option_argument("out") + "' is not recognised");
        }
      }

      if (parser.options.count("suppress") && !mCRL2logEnabled(verbose))
      {
        parser.error("option --suppress requires --verbose (of -v)");
      }

      if (2 < parser.arguments.size())
      {
        parser.error("too many file arguments");
      }
      if (0 < parser.arguments.size())
      {
        m_filename = parser.arguments[0];
      }
      if (1 < parser.arguments.size())
      {
        m_options.lts = parser.arguments[1];
      }

      if (!m_options.lts.empty() && m_options.outformat == lts_none)
      {
        m_options.outformat = mcrl2::lts::detail::guess_format(m_options.lts);

        if (m_options.outformat == lts_none)
        {
          mCRL2log(warning) << "no output format set or detected; using default (mcrl2)" << std::endl;
          m_options.outformat = lts_lts;
        }
      }
    }

};

network2lts_tool *tool_instance;

static
void premature_termination_handler(int)
{
  // Reset signal handlers.
  signal(SIGABRT,NULL);
  signal(SIGINT,NULL);
  signal(SIGTERM,NULL);
  tool_instance->abort();
}

int main(int argc, char** argv)
{
  int result;
  tool_instance = new network2lts_tool();

  signal(SIGABRT,premature_termination_handler);
  signal(SIGINT,premature_termination_handler);
  signal(SIGTERM,premature_termination_handler); // At ^C print a message.

  try
  {
    result = tool_instance->execute(argc, argv);
  }
  catch (...)
  {
    delete tool_instance;
    throw;
  }
  delete tool_instance;
  return result;
}
