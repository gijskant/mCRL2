// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/network.cpp
/// \brief The class network.

#include <string>
#include <cassert>
#include <algorithm>
#include <functional>
#include <fstream>
#include "mcrl2/atermpp/aterm_io.h"
#include "mcrl2/lps/network.h"

namespace mcrl2
{

namespace lps
{

///////////////////////////////////////////////////////////////////////////////
// network
/// \brief network.

/// \brief Constructor.
network::network()
{
}

/// \brief Copy constructor.
network::network(const network &other) :
  m_lps_filenames(other.m_lps_filenames),
  m_synchronization_vector(other.m_synchronization_vector)
{
}

/// \brief Constructor.
network::network( std::vector<std::string> lps_filenames,
                  lps::synchronization_vector synchronization_vector,
                  process::action_label_list action_labels) :
  m_lps_filenames(lps_filenames),
  m_synchronization_vector(synchronization_vector)
{
}

/// \brief Reads the network from a stream.
/// \param in A stream
/// Processes are represented by separate files as specified in the file filename.
void network::read(std::istream& in)
{
  size_t length;
  std::string s;
  if (in.eof())
  {
    std::clog << "Empty synchronization vector stream." << std::endl;
  }
  else
  {
    in >> s;
    //mCRL2log(log::debug) << "network::read: " << s << std::endl;
    // read length
    in >> length;
    //mCRL2log(log::debug) << "length: " << length << std::endl;
    in >> s;
    //mCRL2log(log::debug) << "network::read: " << s << std::endl;
  }
  if (!in.eof())
  {
    std::getline(in, s);
    //mCRL2log(log::debug) << "network::read: " << s << std::endl;
  }
  // read lps filenames
  for (size_t i = 0; !in.eof() && i < length; i++)
  {
    std::string filename;
    std::getline(in, filename);
    if (filename.empty())
    {
      throw std::runtime_error("Empty string. Expected filename.");
    }
    mCRL2log(log::debug) << "lps_filename: " << filename << std::endl;
    m_lps_filenames.push_back(filename);
  }
  assert(m_lps_filenames.size() == length);
  if (m_lps_filenames.size() != length)
  {
    throw std::runtime_error("Number of LPSs does not correspond to the network size.");
  }
  if (in.eof())
  {
    throw std::runtime_error("End of file. Expected synchronization vector.");
  }
  else
  {
    in >> s;
    //mCRL2log(log::debug) << "network::read: " << s << std::endl;
    // read synchronization vector
    m_synchronization_vector.read(in);
    if (m_synchronization_vector.length() != length)
    {
      std::clog << "Synchronization vector length: " << m_synchronization_vector.length() << std::endl;
      std::clog << "Network length: " << length << std::endl;
      throw std::runtime_error("Incompatible size of network and synchronization vector.");
    }
  }
}

/// \brief Reads the network from file.
/// \param filename A string
/// If filename is nonempty, input is read from the file named filename.
/// If filename is empty, input is read from standard input.
/// Processes are represented by separate files as specified in the file filename.
void network::load(const std::string& filename)
{
  std::ifstream instream(filename.c_str(), std::ifstream::in);
  if (!instream)
  {
    throw mcrl2::runtime_error("cannot open input file: " + filename);
  }
  read(instream);
  instream.close();
}

/// \brief Writes the network to a stream.
/// \param out A stream
/// The processes are not written, only the filenames of the original processes.
void network::write(std::ostream& out) const
{
  if (m_synchronization_vector.length() != m_lps_filenames.size())
  {
    std::clog << "Synchronization vector length: " << m_synchronization_vector.length() << std::endl;
    std::clog << "Network length: " << m_lps_filenames.size() << std::endl;
    throw std::runtime_error("Incompatible size of network and synchronization vector.");
  }
  out << "length" << std::endl;
  out << m_lps_filenames.size() << std::endl;
  out << "lps_filenames" << std::endl;
  for(auto it = m_lps_filenames.begin(); it != m_lps_filenames.end(); ++it)
  {
    out << (*it) << std::endl;
  }
  out << "synchronization_vector" << std::endl;
  m_synchronization_vector.write(out);
}

/// \brief Writes the network to file.
/// \param filename A string
/// If filename is nonempty, output is written to the file named filename.
/// If filename is empty, output is written to stdout.
/// The processes are not written to file, only the filenames of the original processes.
void network::save(const std::string& filename) const
{
  std::ofstream outstream(filename.c_str(), std::ofstream::out);
  if (!outstream)
  {
    throw mcrl2::runtime_error("cannot open output file: " + filename);
  }
  write(outstream);
  outstream.close();
}

/// \brief Destructor
network::~network()
{
}

/// \brief Returns the vector of processes.
/// \return The vector processes.
const std::vector<std::string>& network::lps_filenames() const
{
  return m_lps_filenames;
}

/// \brief Returns the vector of processes.
/// \return The vector processes.
std::vector<std::string>& network::lps_filenames()
{
  return m_lps_filenames;
}

/// \brief Returns the synchronization vector.
/// \return The synchronization vector.
const lps::synchronization_vector& network::synchronization_vector() const
{
  return m_synchronization_vector;
}

/// \brief Returns the synchronization vector.
/// \return The synchronization vector.
lps::synchronization_vector& network::synchronization_vector()
{
  return m_synchronization_vector;
}

/// \brief Returns the actions labels.
/// \return The action label list.
const process::action_label_list& network::action_labels() const
{
  return m_action_labels;
}

/// \brief Returns the actions labels.
/// \return The action label list.
process::action_label_list& network::action_labels()
{
  return m_action_labels;
}

std::string pp(const network& x)
{
  std::stringstream s;
  x.write(s);
  return s.str();
}


} // namespace lps

} // namespace mcrl2
