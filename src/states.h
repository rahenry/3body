#ifndef STATES_H
#define STATES_H

#include "general.h"

// contains all information about a particular angular basis state
// this includes various indices, which are used to quickly retrieve elements of caches
struct angular_subsystem;
struct state_angular{
  int index, index_rel, index_E;
  int lr, lrho, lcoup, lcom;
  int nr, nrho;

  state_angular(){};
  state_angular(int index_, int index_rel_, int index_E_, int lr_, int lrho_, int lcoup_, int lcom_, int nr_, int nrho_) : index(index_), index_rel(index_rel_), index_E(index_E_), lr(lr_), lrho(lrho_), lcoup(lcoup_), lcom(lcom_), nr(nr_), nrho(nrho_){};

  vector<int> allowed_transformations; // shows which of states_transformed this state can potentially transform to, based on angular momentum rules
  void generate_allowed_transformations(vector<state_angular> &states_transformed);

  friend ostream &operator<<(ostream &os, state_angular const &s);
};

struct state_full{
  int index;
  state_angular *ang;

  double alpha_r, alpha_rho, alpha_com;
  int channel;
  // channel denotes the odd one out when defining coordinates, with positive and negative corresponding to pos/neg permutations
  //
  double norm;

  void calc_norm();
  state_full() {};
  state_full(int index_, state_angular *ang_, double alpha_r_, double alpha_rho_, double alpha_com_, int channel_) : index(index_), ang(ang_), alpha_r(alpha_r_), alpha_rho(alpha_rho_), alpha_com(alpha_com_), channel(channel_) {calc_norm();};
  state_full(string input_string, angular_subsystem *as);

  state_full(state_full &s_, state_angular *ang_) : index(0), ang(ang_), alpha_r(s_.alpha_r), alpha_rho(s_.alpha_rho), alpha_com(s_.alpha_com), channel(s_.channel) {calc_norm();};

  void permute();
  
  friend ostream &operator<<(ostream &os, state_full const &s) { 
    return os << *s.ang << "  ~|~  " << s.index << " " << s.channel << " " << s.alpha_r << setw(5) << " " << s.alpha_rho << setw(5) << " " << s.alpha_com << setw(5) << endl;
  }

};

void generate_angular_states(vector<state_angular> &tar, int lmax_r, int lmax_rho, int lmax_coup, int lmax_com, int ltotal, string parity, bool single_angular_state, bool use_n_powers, int max_quanta);

int find_angular_index(state_angular *s, angular_subsystem &as);
#endif
