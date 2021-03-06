// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file formulaparelm.cpp

#include "mcrl2/modal_formula/tools/formularewr.h"
#include "mcrl2/utilities/input_output_tool.h"

using namespace mcrl2;
using namespace mcrl2::log;
using namespace mcrl2::state_formulas;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;

class formula_rewr_tool: public input_output_tool
{
    typedef input_output_tool super;

protected:
  std::string lps_filename;
  std::string network_filename;
  formularewr_options options;

  std::string synopsis() const
  {
    return "[OPTION]... [--lps=LPS|--network=NETWORK] [INFILE [OUTFILE]]\n";
  }

  void add_options(interface_description& desc)
  {
    super::add_options(desc);
    desc.add_option("lps", make_file_argument("LPS"),
                    "use the linear process from file LPS", 'p');
    desc.add_option("network", make_file_argument("NETWORK"),
                    "use the network of LPSs from file NETWORK", 'n');
    desc.add_option("unfold-recursion",
                    "unfold unguarded recursion", 'U');
  }

  void parse_options(const command_line_parser& parser)
  {
    super::parse_options(parser);
    if (parser.options.count("lps"))
    {
      lps_filename = parser.option_argument("lps");
    }
    if (parser.options.count("network"))
    {
      network_filename = parser.option_argument("network");
    }
    options.unfold_recursion = parser.options.count("unfold-recursion") != 0;
  }


  public:
    formula_rewr_tool()
      : super(
        "formularewr",
        "Gijs Kant",
        "rewrites a modal formula",
        "Reads a file containing a modal formula, and applies the one point rule for quantifiers and "
        "simplifies data expressions. If OUTFILE is not present, standard output is used. "
        "If INFILE is not present, standard input is used. "
        "The data specification in LPS or NETWORK is used for parsing the formula."
      )
    {}

    bool run() /*< The virtual function `run` executes the tool.
                   The user has to override this function to add behavior. >*/
    {
      mCRL2log(verbose) << "formularewr parameters:" << std::endl;
      mCRL2log(verbose) << "  input file:         " << m_input_filename << std::endl;
      mCRL2log(verbose) << "  output file:        " << m_output_filename << std::endl;
      mCRL2log(verbose) << "  lps file:           " << lps_filename << std::endl;
      mCRL2log(verbose) << "  network file:       " << network_filename << std::endl;

      formularewr(input_filename(),
                 output_filename(),
                 lps_filename,
                 network_filename,
                 options
                );

      return true;
    }
};

int main(int argc, char* argv[])
{
  return formula_rewr_tool().execute(argc, argv); /*< The function `execute` first parses the command line
                                       arguments, and then calls the function `run`. >*/
}
//]
