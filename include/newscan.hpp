#ifndef NEWSCAN_HPP
#define NEWSCAN_HPP

/* ******************************************************************************
 * newscan.hpp
 *
 * parsing algorithm for bwt construction of repetitive sequences based
 * on prefix free parsing. See:
 *   Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini
 *   Prefix-Free Parsing for Building Big BWTs
 *   [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/)
 *
 * Modified by: Vikram Shivakumar
 * Converted into a header file for use in Mumemto
 */
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
#include <stdio.h>
using namespace std;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;


// =============== constants ===================
#define Dollar 2     // special char for the parsing algorithm, must be the highest special char 
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter
#define EXTPARSE "parse"
#define EXTPARS0 "parse_old"
#define EXTDICT  "dict"

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank=0;
};

void die(const string& message) {
  cerr << message << endl;
  exit(1);
}
FILE *open_aux_file(const char *base, const char *ext, const char *mode)
{
  std::string filename = std::string(base) + "." + std::string(ext);
  FILE *f = fopen(filename.c_str(),mode);
  if(f==NULL) die(filename);  
  return f;
}

// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
  int *window;
  int asize;
  const uint64_t prime = 1999999973; //old 1999999973
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime

  KR_window(int w): wsize(w) {
    asize = 256;
    asize_pot = 1;
    for(int i=1;i<wsize;i++)
      asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm
    // alloc and clear window
    window = new int[wsize];
    reset();
  }

  // init window, hash, and tot_char
  void reset() {
    for(int i=0;i<wsize;i++) window[i]=0;
    // init hash value and related values
    hash=tot_char=0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    // complex expression to avoid negative numbers
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution
    hash = (asize*hash + c) % prime;      //  add char i
    window[k]=c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash;
  }
  // debug only
  string get_window() {
    string w = "";
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }

  ~KR_window() {
    delete[] window;
  }

};
// -----------------------------------------------------------

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k=0;k<s.size();k++) {
      int c = (unsigned char) s[k];
      hash = (256*hash + c) % prime;    //  add char k
    }
    return hash;
}


// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a < *b;
}

bool pword_statsCompare(const word_stats *a, const word_stats *b)
{
  return a->str < b->str;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(string &inputFileName, map<uint64_t,word_stats> &wfreq, vector<word_stats *> &sortedDict)
{
  assert(sortedDict.size() == wfreq.size());
  FILE *fdict;
  // open dictionary and occ files
  fdict = open_aux_file(inputFileName.c_str(),EXTDICT,"wb");

  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    const char *word = x->str.data();       // current dictionary word
    size_t len = x->str.size();             // length of word
    size_t s = fwrite(word, 1, len, fdict);
    if(s!=len) die("Error writing to DICT file");
    if(fputc(EndOfWord,fdict)==EOF) die("Error writing EndOfWord to DICT file");
    x->rank = wrank;
    wrank++;
  }
  if(fputc(EndOfDict,fdict)==EOF) die("Error writing EndOfDict to DICT file");
  if(fclose(fdict)!=0) die("Error closing DICT file");
}

void remapParse(std::string &inputFileName, map<uint64_t,word_stats> &wfreq)
{
  // open parse files. the old parse can be stored in a single file or in multiple files
  FILE *old_parse_file = open_aux_file(inputFileName.c_str(),EXTPARS0,"rb");
  FILE *new_parse_file = open_aux_file(inputFileName.c_str(),EXTPARSE,"wb");

  // recompute occ as an extra check
  uint64_t hash;
  while(true) {
    size_t s = fread(&hash,sizeof(hash),1,old_parse_file);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    word_int_t rank = wfreq.at(hash).rank;
    s = fwrite(&rank,sizeof(rank),1,new_parse_file);
    if(s!=1) die("Error writing to new parse file");
  }
  if(fclose(new_parse_file)!=0) die("Error closing new parse file");
  if(fclose(old_parse_file)!=0) die("Error closing old parse file");
}

// =============== pfparser class ===================
class pfparser {
private:
    map<uint64_t, word_stats> wordFreq;
    string file_prefix;
    FILE* tmp_parse_file;

    // intermediate values for parsing
    string word;
    KR_window krw;

    // arguments
    bool probing;
    size_t w;
    size_t p;
    
    // Helper method to save and update words
    void save_update_word(string& window, uint64_t hash);

public:
    pfparser(string file_prefix, size_t w, size_t p, bool probing);
    ~pfparser();
    
    // Process a single string and accumulate results
    void process_string(const string& input_string);
        
    // Run the second pass (dictionary construction and remapping)
    void finish_parse();

};

// Constructor
pfparser::pfparser(string file_prefix, size_t w, size_t p, bool probing) :
      file_prefix(file_prefix),
      w(w), 
      p(p), 
      probing(probing),
      krw(w) {
    // Initialize any necessary state
    
    tmp_parse_file = open_aux_file(file_prefix.c_str(),EXTPARS0,"wb");
    // init first word in the parsing with a Dollar char
    word.append(1,Dollar);
}

// Destructor - clean up any open files
pfparser::~pfparser() {}

// Modified save_update_word as a member function
void pfparser::save_update_word(string& window, uint64_t hash) {
    if(window.size() <= w) return;
        
    // Update frequency table for current hash
    if (probing) {
        auto slot = wordFreq.find(hash);
        while(slot != wordFreq.end() && wordFreq[hash].str != window) { 
            hash += 1;
            slot = wordFreq.find(hash);
        }
        if (slot == wordFreq.end()) {
            wordFreq[hash].occ = 1;
            wordFreq[hash].str = window;
        }
        else {
            wordFreq[hash].occ += 1;
        }
    }
    else {
        if(wordFreq.find(hash) == wordFreq.end()) {
            wordFreq[hash].occ = 1;
            wordFreq[hash].str = window;
        }
        else {
            wordFreq[hash].occ += 1;
            if(wordFreq[hash].occ <= 0) {
                cerr << "Emergency exit! Maximum # of occurrence of dictionary word (";
                cerr << MAX_WORD_OCC << ") exceeded\n";
                exit(1);
            }
            if(wordFreq[hash].str != window) {
                cerr << "Emergency exit! Hash collision for strings:\n";
                cerr << wordFreq[hash].str << "\n  vs\n" <<  window << endl;
                exit(1);
            }
        }
    }

    if(fwrite(&hash, sizeof(hash), 1, tmp_parse_file) != 1) die("parse write error");
      
    // Keep only the overlapping part of the window
    window.erase(0, window.size() - w);
}

// Modified process_string as a member function - takes a string instead of filename
void pfparser::process_string(const string& input_string) {
    uint64_t hash;
    int c;
    // Main processing logic for the input string    
    // Process each character in the input string
    for (size_t i = 0; i < input_string.length(); i++) {
        c = input_string[i];
        if(c <= Dollar) { cerr << "Invalid char found in input string: no additional chars will be read\n"; break;}
        word.append(1, c);
        hash = krw.addchar(c);
        if(hash % p == 0) {
            save_update_word(word, hash);
        }
    }
}

// Run the second pass
void pfparser::finish_parse() {
    // virtually add w null chars at the end of the file and add the last word in the dict
    word.append(w, Dollar);
    uint64_t hash = kr_hash(word);
    save_update_word(word, hash);
    // close input and output files
    if(fclose(tmp_parse_file)!=0) die("Error closing parse file");
    
    // Check # distinct words
    uint64_t totDWord = wordFreq.size();
    if(totDWord > MAX_DISTINCT_WORDS) {
        cerr << "Emergency exit! The number of distinct words (" << totDWord << ")\n";
        cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
        exit(1);
    }
    
    // Create array of dictionary words
    vector<word_stats *> dictArray;
    dictArray.reserve(totDWord);
    
    // Fill array
    for (auto& x: wordFreq) {
        dictArray.push_back(&x.second);
    }
      
    // Sort dictionary
    sort(dictArray.begin(), dictArray.end(), pword_statsCompare);
    
    // Write plain dictionary and occ file, also compute rank for each hash
    writeDictOcc(file_prefix, wordFreq, dictArray);
    dictArray.clear(); // reclaim memory
    
    // Remap parse file
    remapParse(file_prefix, wordFreq);
}


#endif // NEWSCAN_HPP

