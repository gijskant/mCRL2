// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/rewriter_tool.h
/// \brief Base class for tools that use a data rewriter.

#ifndef MCRL2_UTILITIES_REWRITER_TOOL_H
#define MCRL2_UTILITIES_REWRITER_TOOL_H

#include "mcrl2/data/rewrite_strategy.h"
#include "mcrl2/data/rewriter.h"
#include "mcrl2/utilities/command_line_interface.h"
#include "mcrl2/data/detail/enumerator_variable_limit.h"

namespace mcrl2
{

namespace utilities
{

namespace tools
{

/// \brief Base class for tools that use a rewriter.
template <typename Tool>
class rewriter_tool: public Tool
{
  protected:
    /// The data rewriter strategy
    data::rewrite_strategy m_rewrite_strategy;

    /// \brief Add options to an interface description. Also includes
    /// rewriter options.
    /// \param desc An interface description
    void add_options(interface_description& desc)
    {
      Tool::add_options(desc);

      desc.add_option(
        "rewriter", make_enum_argument<data::rewrite_strategy>("NAME")
            .add_value(data::jitty, true)
#ifdef MCRL2_JITTYC_AVAILABLE
            .add_value(data::jitty_compiling)
#endif
            .add_value(data::jitty_prover),
        "use rewrite strategy NAME:"
        ,'r'
      );

      desc.add_option(
        "qlimit", make_mandatory_argument("NUM"),
        "limit enumeration of quantifiers to NUM variables. (Default NUM=1000, NUM=0 for unlimited).",
        'Q'
      );

    }

    /// \brief Parse non-standard options
    /// \param parser A command line parser
    void parse_options(const command_line_parser& parser)
    {
      Tool::parse_options(parser);
      m_rewrite_strategy = parser.option_argument_as< data::rewrite_strategy >("rewriter");

      if(parser.options.count("qlimit"))
      {
        //Set enumerator limit for quantifier enumeration
        data::detail::set_enumerator_variable_limit(parser.option_argument_as< size_t >("qlimit"));
      }
    }

  public:

    /// \brief Constructor.
    rewriter_tool(const std::string& name,
                  const std::string& author,
                  const std::string& what_is,
                  const std::string& tool_description,
                  std::string known_issues = ""
                 )
      : Tool(name, author, what_is, tool_description, known_issues),
        m_rewrite_strategy(mcrl2::data::jitty)
    {}

    /// \brief Returns the rewrite strategy
    /// \return The rewrite strategy
    data::rewrite_strategy rewrite_strategy() const
    {
      return data::rewrite_strategy(m_rewrite_strategy);
    }

    /// \brief Creates a data rewriter as specified on the command line.
    /// \param data_spec A data specification
    /// \return A data rewriter
    data::rewriter create_rewriter(data::data_specification const& data_spec = data::data_specification())
    {
      return data::rewriter(data_spec, rewrite_strategy());
    }
};

} // namespace tools

} // namespace utilities

} // namespace mcrl2

#endif // MCRL2_UTILITIES_REWRITER_TOOL_H