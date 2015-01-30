// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file formulaquotient.cpp
/// \brief Add your file description here.

#include <iostream>
#include <string>
#include "mcrl2/pbes/tools/formulaquotient.h"
#include "mcrl2/utilities/text_utility.h"
#include "mcrl2/utilities/input_output_tool.h"

using namespace mcrl2;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::log;

class formulaquotient_tool : public input_output_tool
{
    typedef input_output_tool super;

  protected:
    std::string input_network_filename;
    std::string output_network_filename;
    std::string lps_filename;
    std::string synchronization_vector_filename;
    size_t i;
    state_formulas::formulaquotient_options tool_options;
    state_formulas::quotienting_options algorithm_options;

    std::string synopsis() const
    {
      return "[OPTION]... [--network=NETWORK [--out-network=OUTNETWORK]| --lps=LPS --vector=VECTOR] --index=INDEX [INFILE [OUTFILE]]\n";
    }

    void add_options(interface_description& desc)
    {
      super::add_options(desc);
      desc.add_option("network", make_file_argument("NETWORK"),
                      "use the network from file NETWORK", 'n');
      desc.add_option("output-network", make_file_argument("OUTNETWORK"),
                      "write the result network to file OUTNETWORK", 'o');
      desc.add_option("lps", make_file_argument("LPS"),
                      "use the linear process from file LPS", 'p');
      desc.add_option("vector", make_file_argument("VECTOR"),
                      "use the synchronization vector from file VECTOR", 's');
      desc.add_option("index", make_mandatory_argument("INDEX"),
                      "the index of the LPS in the network (default: 0)", 'i');
      desc.add_option("iterative",
                      "iteratively applies quotienting (output will be a PBES)", 'I');
      desc.add_option("parelm",
                      "applies unused parameter elimination after quotienting", 'P');
      desc.add_option("simplify",
                      "simplifies formulae after quotienting", 'S');
      desc.add_option("unfold-recursion",
                      "unfold unguarded recursion", 'U');
      desc.add_option("write-formulas",
                      "saves intermediate formulas to file", 'F');
      desc.add_option("approximate",
		      "solves approximation of intermediate formulas, requires -F", 'A');
      desc.add_option("write-networks",
                      "saves intermediate networks to file", 'N');
      desc.add_option("disable-cache-rewriter",
                      "disables caching of rewriter results", 'D');
      desc.add_option("use-vector-map",
                      "filter synchronization vector for relevant actions", 'M');
    }

    void parse_options(const command_line_parser& parser)
    {
      super::parse_options(parser);
      if (parser.options.count("network"))
      {
        input_network_filename = parser.option_argument("network");
        if (parser.options.count("output-network"))
        {
          output_network_filename = parser.option_argument("output-network");
        }
      }
      else
      {
        if (parser.options.count("lps"))
        {
          lps_filename = parser.option_argument("lps");
        }
        if (parser.options.count("vector"))
        {
          synchronization_vector_filename = parser.option_argument("vector");
        }
      }
      tool_options.iterative = parser.options.count("iterative") != 0;
      tool_options.parelm = parser.options.count("parelm") != 0;
      tool_options.simplify = parser.options.count("simplify") != 0;
      tool_options.unfold_recursion = parser.options.count("unfold-recursion") != 0;
      tool_options.write_intermediate_formulas = parser.options.count("write-formulas") != 0;
      if (tool_options.write_intermediate_formulas)
      {
	tool_options.approximate_intermediate_formulas = parser.options.count("approximate") != 0;
      }
      tool_options.write_intermediate_networks = parser.options.count("write-networks") != 0;
      algorithm_options.use_rewrite_cache = parser.options.count("disable-cache-rewriter") == 0;
      algorithm_options.use_vector_map = parser.options.count("use-vector-map") != 0;

      if (parser.options.count("index"))
      {
        i = size_t(std::stoi(parser.option_argument("index")));
      }
      else
      {
        i = 0;
      }
    }

  public:
    formulaquotient_tool() : super(
        "formulaquotient",
        "Gijs Kant",
        "compute a quotient formula from a state formula and an LPS",
        "Compute a quotient formula from the state formula in INFILE "
        "by quotienting out the linear process (LPS) in LPS or NETWORK "
        "with index INDEX using the synchronization vector in VECTOR, "
        "save it to OUTFILE and write the resulting network to "
        "OUTNETWORK. If OUTFILE is not "
        "present, stdout is used. If INFILE is not present, stdin is used."
      ),
      i(0)
    {}

    bool run()
    {
      state_formulas::formulaquotient(input_filename(),
               output_filename(),
               input_network_filename,
               output_network_filename,
               lps_filename,
               synchronization_vector_filename,
               i,
               tool_options,
               algorithm_options
             );
      return true;
    }

};

int main(int argc, char** argv)
{
  return formulaquotient_tool().execute(argc, argv);
}
