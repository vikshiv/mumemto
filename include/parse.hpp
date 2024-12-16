/* pfp-parse - prefix free parsing parse
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file parse.hpp
   \brief parse.hpp define and build the prefix-free parse data structure.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note This is a short version of the parse in https://github.com/maxrossi91/pfp-data-structures
*/

#ifndef _PFP_PARSE_HH
#define _PFP_PARSE_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}

// TODO: Extend it to non-integer alphabets
class parse{
public:
  std::vector<uint32_t> p;
  std::vector<uint_t> saP;
  sdsl::int_vector<40> isaP;
  // std::vector<uint_t> isaP;

  sdsl::int_vector<40> ilist; // Inverted list of phrases of P in BWT_P
  sdsl::bit_vector ilist_s; // The ith 1 is in correspondence of the first occurrence of the ith phrase
  sdsl::bit_vector::select_1_type select_ilist_s;

  size_t alphabet_size;

  typedef size_t size_type;

  // Default constructor for load
  parse() {}

  parse(  std::vector<uint32_t>& p_,
          size_t alphabet_size_):
          p(p_),
          alphabet_size(alphabet_size_)
  {
    assert(p.back() == 0);

    compute_freq();

    build();


  }

  parse(  std::string filename,
          size_t alphabet_size_):
          alphabet_size(alphabet_size_)
  {
    // Building dictionary from file
    std::string tmp_filename = filename + std::string(".parse");
    read_file(tmp_filename.c_str(), p);
    p.push_back(0); // this is the terminator for the sacak algorithm

    // // Uploading the frequency file
    // tmp_filename = filename + std::string(".occ");
    // read_file(tmp_filename.c_str(), freq);
    // freq.insert(freq.begin(), 1);

    // compute_freq();

    build();


  }

  void build(){

    saP.resize(p.size());
    // suffix array of the parsing.
    verbose("Computing SA of the parsing");
    _elapsed_time(
      sacak_int(&p[0],&saP[0],p.size(),alphabet_size);
    );



    // inverted list of the parsing.
    // verbose("Computing ilist");
    // _elapsed_time(
    //   compute_ilist()
    // );


    // inverse suffix array of the parsing.
    verbose("Computing ISA of the parsing");
    _elapsed_time(
      {
        isaP.resize(p.size());
        for(int i = 0; i < saP.size(); ++i){
          isaP[saP[i]] = i;
        }
      }
    );


  }


  void compute_ilist()
  {
    ilist.resize(p.size());
    ilist_s = sdsl::bit_vector(p.size() + 1, 0);
    // computing the bucket boundaries
    size_t j = 0;
    ilist_s[j++] = 1; // this is equivalent to set pf.freq[0]=1;
    ilist_s[j] = 1;
    for (size_t i = 1; i < freq.size(); ++i)
    {
      j += freq[i];
      ilist_s[j] = 1;
      freq[i] = 0;
    }
    std::cerr << "last"<<std::endl;

    select_ilist_s = sdsl::bit_vector::select_1_type(&ilist_s);
    std::cerr << "hi"<<std::endl;

    for(size_t i = 0; i < saP.size(); ++i)
    {
      size_t prec_phrase_index = (saP[i] == 0 ? p.size() : saP[i]) - 1;
      uint_t prec_phrase = p[prec_phrase_index];

      size_t ilist_p = select_ilist_s(prec_phrase + 1) + freq[prec_phrase]++;
      ilist[ilist_p] = i;
    }
    std::cerr << "hi"<<std::endl;
    freq.clear();
    freq.shrink_to_fit();
    // sdsl::int_vector<40>().swap(freq);
  }

  void compute_freq()
  {
    uint32_t max_p = p[0];
    for(size_t i = 1; i < p.size(); ++i)
      max_p = std::max(max_p, p[i]);

    freq.resize(max_p+1);
    for(size_t i = 0; i < p.size(); ++i)
      freq[p[i]]++;

  }

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += my_serialize(p, out, child, "parse");
    written_bytes += my_serialize(saP, out, child, "saP");
    // written_bytes += my_serialize(isaP, out, child, "isaP");
    // written_bytes += my_serialize(ilist, out, child, "ilist");
    written_bytes += ilist_s.serialize(out, child, "ilist_s");
    written_bytes += select_ilist_s.serialize(out, child, "select_ilist_s");
    written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");
    

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    my_load(p, in);
    my_load(saP, in);
    // my_load(isaP, in);
    // my_load(ilist, in);
    ilist_s.load(in);
    select_ilist_s.load(in);
    sdsl::read_member(alphabet_size, in);
  }

private:
  // sdsl::int_vector<40> freq;
  std::vector<uint_t> freq;
};

#endif /* end of include guard: _PFP_PARSE_HH */
