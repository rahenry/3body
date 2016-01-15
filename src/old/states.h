#ifndef STATES_H
#define STATES_H

#include "general.h"

// contains all information about a particular angular basis state
// this includes various indices, which are used to quickly retrieve elements of caches
struct state_angular{
  int index, index_rel, index_E;
  int lr, lrho, lcoup, lcom;
  int nr, nrho;

  state_angular(int index_, int index_rel_, int index_E_, int lr_, int lrho_, int lcoup_, int lcom_, int nr_, int nrho_) : index(index_), index_rel(index_rel_), index_E(index_E_), lr(lr_), lrho(lrho_), lcoup(lcoup_), lcom(lcom_), nr(nr_), nrho(nrho_){};

  vector<int> allowed_transformations; // shows which of states_transformed this state can potentially transform to, based on angular momentum rules
  void generate_allowed_transformations(vector<state_angular> &states_transformed);

  friend ostream &operator<<(ostream &os, state_angular const &s);
};

struct state_full{
  int index;
  state_angular *angular;

  double alpha_r, alpha_rho, alpha_R;
  int channel;

  state_full(int index_, state_angular *angular_, double alpha_r_, double alpha_rho_, double alpha_R_, int channel_) : index(index_), angular(angular_), alpha_r(alpha_r_), alpha_rho(alpha_rho_), alpha_R(alpha_R_), channel(channel_) {};
  
  friend ostream &operator<<(ostream &os, state_full const &s) { 
    return os << *s.angular << "  ~|~  " << s.channel << " " << s.alpha_r << setw(5) << " " << s.alpha_rho << setw(5) << " " << s.alpha_R << setw(5) << endl;
  }
};

void generate_angular_states(vector<state_angular> &tar, int lmax_r, int lmax_rho, int lmax_coup, int lmax_com, int ltotal, string parity, bool single_angular_state, bool use_n_powers, int max_quanta);

#endif
