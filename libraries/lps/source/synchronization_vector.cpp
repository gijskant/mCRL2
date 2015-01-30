// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/synchronization_vector.cpp
/// \brief

#include <string>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/atermpp/aterm_list.h"
#include "mcrl2/data/data_specification.h"
#include "mcrl2/lps/parse.h"
#include "mcrl2/lps/synchronization_vector.h"
#include "mcrl2/process/normalize_sorts.h"

namespace mcrl2
{

namespace lps
{

///////////////////////////////////////////////////////////////////////////////
// synchronization_vector
/// \brief synchronization vector.

synchronization_vector::synchronization_vector() :
  m_length(0)
{
}

/// \brief Copy constructor.
synchronization_vector::synchronization_vector(const synchronization_vector &other) :
  m_synchronization_vector(other.m_synchronization_vector),
  m_action_labels(other.m_action_labels),
  m_length(other.m_length)
{
}

/// \brief Constructor.
synchronization_vector::synchronization_vector(int length) :
  m_length(length)
{
}

/// \brief Constructor.
synchronization_vector::synchronization_vector(synchronization_vector_type synchronization_vector,
    std::set<process::action_label> action_labels,
    int length) :
  m_synchronization_vector(synchronization_vector),
  m_action_labels(action_labels),
  m_length(length)
{
}

enum parse_state {
  list, vector, label
};

std::string print_vector(std::vector<std::string> x)
{
  std::stringstream ss;
  bool first = true;
  for(auto it = x.begin(); it != x.end(); ++it)
  {
    if (first){
      first = false;
    } else
    {
      ss << ", ";
    }
    ss << (*it);
  }
  return ss.str();
}

void synchronization_vector::normalize(const data::data_specification& dataspec)
{
  for (auto it = m_synchronization_vector.begin(); it != m_synchronization_vector.end(); ++it)
  {
    process::action_label label = it->second;
    //mCRL2log(log::verbose) << "normalize: replacing \"" << pp(label.sorts()) << "\"" << std::flush;
    label = process::normalize_sorts(label, dataspec);
    //mCRL2log(log::verbose) << " with \"" << pp(label.sorts()) << "\"" << std::endl;
    it->second = label;
  }
  std::set<process::action_label> label_set;
  for (auto it = m_action_labels.begin(); it != m_action_labels.end(); ++it)
  {
    process::action_label label = *it;
    label = process::normalize_sorts(label, dataspec);
    label_set.insert(label);
  }
  m_action_labels = label_set;
}

/// \brief Reads the synchronization vector from stream.
/// \param in An input stream
void synchronization_vector::read(std::istream& in)
{
  if (in.eof())
  {
    std::clog << "Empty synchronization vector stream." << std::endl;
  }
  else
  {
    in >> m_length;
  }
  mCRL2log(log::debug) << "length: " << m_length << std::endl;
  char c;
  size_t i = 0;
  std::stringstream ss;
  std::vector<std::string> v;
  std::string label_text;
  parse_state state = parse_state::list;
  while (!in.eof())
  {
    in >> c;
    if (c == '{')
    {
      break;
    }
    else if (!std::isspace(c))
    {
      throw std::runtime_error(std::string("Unexpected character: ") + c);
    }
  }
  //mCRL2log(log::debug) << "start reading vector" << std::endl;
  size_t paren_level = 0;
  bool done = false;
  while (!in.eof() && !done)
  {
    in >> c;
    //mCRL2log(log::debug) << "reading char: " << c << std::endl;
    switch(state) {
    case parse_state::list:
      if (c == '}')
      {
        done = true;
      }
      else if (c == '(')
      {
        i = 0;
        ss.str("");
        paren_level = 0;
        //mCRL2log(log::debug) << "start of vector" << std::endl;
        state = parse_state::vector;
      }
      else if (c == ',')
      {
      }
      else if (!std::isspace(c))
      {
        throw std::runtime_error(std::string("Unexpected character: ") + c);
      }
      break;
    case parse_state::vector:
      if (c == ',')
      {
        label_text = ss.str();
        ss.str("");
        if (label_text.empty())
        {
          throw std::runtime_error("Empty label name.");
        }
        else
        {
          //mCRL2log(log::debug) << "adding label: " << label_text << std::endl;
          v.push_back(label_text);
          i++;
          if (i == m_length)
          {
            //mCRL2log(log::debug) << "vector: " << print_vector(v) << std::endl;
            state = parse_state::label;
          }
        }
      }
      else if (c == '(' || c == ')')
      {
        throw std::runtime_error(std::string("Unexpected character: ") + c);
      }
      else if (!std::isspace(c))
      {
        ss << c;
      }
      break;
    case parse_state::label:
      if (c == ')' && paren_level==0)
      {
        label_text = ss.str();
        ss.str("");
        if (label_text.empty())
        {
          throw std::runtime_error("Empty label type.");
        }
        else
        {
          mCRL2log(log::debug) << "parsing label '" << label_text << "'..." << std::endl;
          process::action_label_list labels = process::parse_action_declaration(label_text + ";");
          if (labels.size() != 1)
          {
            throw std::runtime_error("Invalid label declaration, should be singleton: " + label_text);
          }
          process::action_label label = labels.front();
          add_vector(v, label);
          //mCRL2log(log::debug) << "vector read: v = " << print_vector(v) << ", label = " << pp(label) << std::endl;
          v.clear();
          state = parse_state::list;
        }
      }
      else
      {
        if (c == '(')
        {
          paren_level++;
        }
        else if (c == ')')
        {
          paren_level--;
        }
        ss << c;
      }
      break;
    }
  }
}


/// \brief Reads the synchronization vector from file.
/// \param filename A string
/// If filename is nonempty, input is read from the file named filename.
/// If filename is empty, input is read from standard input.
void synchronization_vector::load(const std::string& filename)
{
  std::ifstream instream(filename.c_str(), std::ifstream::in);
  if (!instream)
  {
    throw mcrl2::runtime_error("cannot open input file: " + filename);
  }
  read(instream);
  instream.close();
}

/// \brief Writes the synchronization vector to a stream.
/// \param out An output stream
void synchronization_vector::write(std::ostream& out) const
{
  out << m_length << std::endl;
  /*
  for(auto v_it = m_synchronization_vector.begin(); v_it != m_synchronization_vector.end(); ++v_it)
  {
    std::pair<std::vector<std::string>, action_label > p = (*v_it);
    out << p.second.name() << ":";
    core::detail::apply_printer<data::detail::printer> printer(out);
    printer.print_list(p.second.sorts(), "", "","#");
    out << " ";
    std::copy(p.first.begin(), p.first.end(), std::ostream_iterator<std::string>(out, " "));
    out << std::endl;
  }*/
  out << "{";
  for(auto it = m_synchronization_vector.begin(); it != m_synchronization_vector.end(); ++it)
  {
    out << std::endl << "  ( ";
    synchronization_vector_element element = *it;
    std::vector<std::string> labels = element.first;
    for(auto label_it = labels.begin(); label_it != labels.end(); ++label_it)
    {
        out << *label_it << ", ";
    }
    out << element.second.name();
    out << " )";
  }
  out << std::endl << "}" << std::endl;
}

/// \brief Writes the synchronization vector to file.
/// \param filename A string
/// If filename is nonempty, output is written to the file named filename.
/// If filename is empty, output is written to stdout.
void synchronization_vector::save(const std::string& filename) const
{
  std::clog << "Writing synchronization vector." << std::endl;
  std::ofstream outstream(filename.c_str(), std::ofstream::out);
  if (!outstream)
  {
    throw mcrl2::runtime_error("cannot open output file: " + filename);
  }
  write(outstream);
  outstream.close();
}

/// \brief Destructor
synchronization_vector::~synchronization_vector()
{
}

void synchronization_vector::add_vector(std::vector<std::string> actions, process::action_label a)
{
  if (actions.size()==m_length)
  {
    m_synchronization_vector.push_back(synchronization_vector_element(actions, a));
    m_action_labels.insert(a);
  }
  else
  {
    std::clog << "Vector has length " << actions.size() << ", should be " << m_length << std::endl;
    throw std::runtime_error("vector of actions has invalid length.");
  }
}

const synchronization_vector_type& synchronization_vector::vector() const
{
  return m_synchronization_vector;
}

const std::set<process::action_label>& synchronization_vector::action_labels() const
{
  return m_action_labels;
}

size_t synchronization_vector::length() const
{
  return m_length;
}

} // namespace lps

} // namespace mcrl2
