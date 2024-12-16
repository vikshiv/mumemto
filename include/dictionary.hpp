/* pfp-dictionary - prefix free parsing dictionary
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
   \file dictionary.hpp
   \brief dictionary.hpp define and build the prefix-free dictionary data structure.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note This is a short version of the dictionary in https://github.com/maxrossi91/pfp-data-structures
*/

#ifndef _PFP_DICTIONARY_HH
#define _PFP_DICTIONARY_HH

#include <queue>

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}


// TODO: Extend it to integer alphabets
class dictionary{
public:
  std::vector<uint8_t> d;
  std::vector<uint_t> saD;
  sdsl::int_vector<5> isaD;
  std::vector<int_t> lcpD;
  sdsl::rmq_succinct_sct<> rmq_lcp_D;
  sdsl::bit_vector b_d; // Starting position of each phrase in D
  sdsl::bit_vector::rank_1_type rank_b_d;
  sdsl::bit_vector::select_1_type select_b_d;

  std::vector<uint8_t> alphabet;

  typedef size_t size_type;

  // default constructor for load.
  dictionary() {}

  dictionary( std::vector<uint8_t>& d_,
              size_t w ):
              d(d_)
  {
    build();

  }

  dictionary(std::string filename,
             size_t w)
  {
    // Building dictionary from file
    std::string tmp_filename = filename + std::string(".dict");

    read_file(tmp_filename.c_str(), d);
    assert(d[0] == Dollar);

    // Prepending w dollars to d
    // 1. Count how many dollars there are
    int i = 0;
    int n_dollars = 0;

    while(i < d.size() && d[i++] == Dollar)
      ++n_dollars;

    std::vector<uint8_t> dollars(w-n_dollars,Dollar);
    d.insert(d.begin(), dollars.begin(),dollars.end());

    build();
  }

  inline size_t length_of_phrase(size_t id){
    assert(id > 0);
    return select_b_d(id+1)-select_b_d(id) - 1; // to remove the EndOfWord
  }

  inline size_t n_phrases(){
    return rank_b_d(d.size()-1);
  }

  size_t longest_common_phrase_prefix(size_t a, size_t b)
  {
    if(a == 0 || b == 0)
      return 0;
    // Compute the lcp between phrases a and b
    auto a_in_sa = isaD[select_b_d(a)]; // position of the phrase a in saD
    auto b_in_sa = isaD[select_b_d(b)]; // position of the phrase b in saD

    auto lcp_left = std::min(a_in_sa, b_in_sa) + 1;
    auto lcp_right = std::max(a_in_sa, b_in_sa);

    size_t lcp_a_b_i = rmq_lcp_D(lcp_left, lcp_right);
    return lcpD[lcp_a_b_i];

  }

  void build(){
    // Constructing the alphabet
    std::vector<bool> visit(256,false);
    for(auto elem: d)
    {
      if(!visit[elem])
      {
        visit[elem] = true;
        alphabet.push_back(elem);
      }
    }

    // Building the bitvector with a 1 in each starting position of each phrase in D
    b_d.resize(d.size());
    for(size_t i = 0; i < b_d.size(); ++i) b_d[i] = false; // bug in resize
    b_d[0] = true; // Mark the first phrase
    for(size_t i = 1; i < d.size(); ++i )
      b_d[i] = (d[i-1]==EndOfWord);
    b_d[d.size()-1] = true; // This is necessary to get the length of the last phrase

    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    select_b_d = sdsl::bit_vector::select_1_type(&b_d);

    saD.resize(d.size());
    lcpD.resize(d.size());
    // daD.resize(d.size());
    // suffix array, LCP array, and Document array of the dictionary.

    //verbose("Computing SA, LCP, and DA of dictionary");
    //_elapsed_time(
      gsacak(&d[0], &saD[0], &lcpD[0], nullptr, d.size());
      // gsacak(&d[0], &saD[0], &lcpD[0], &daD[0], d.size())
    //);


    // inverse suffix array of the dictionary.
    verbose("Computing ISA of dictionary");
    _elapsed_time(
      {
        isaD.resize(d.size());
        for(int i = 0; i < saD.size(); ++i){
          isaD[saD[i]] = i;
        }
      }
    );

    

    verbose("Computing RMQ over LCP of dictionary");
    // Compute the LCP rank of D
    _elapsed_time(
      rmq_lcp_D = sdsl::rmq_succinct_sct<>(&lcpD)
    );
    
  }

  

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += my_serialize(d, out, child, "dictionary");
    written_bytes += my_serialize(saD, out, child, "saD");
    // written_bytes += my_serialize(isaD, out, child, "isaD");
    written_bytes += my_serialize(lcpD, out, child, "lcpD");
    written_bytes += rmq_lcp_D.serialize(out, child, "rmq_lcp_D");
    written_bytes += b_d.serialize(out, child, "b_d");
    written_bytes += rank_b_d.serialize(out, child, "rank_b_d");
    written_bytes += select_b_d.serialize(out, child, "select_b_d");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;

  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    my_load(d, in);
    my_load(saD, in);
    // my_load(isaD, in);
    my_load(lcpD, in);
    rmq_lcp_D.load(in);
    b_d.load(in);
    rank_b_d.load(in, &b_d);
    select_b_d.load(in, &b_d);
  }
};





#endif /* end of include guard: _PFP_DICTIONARY_HH */
