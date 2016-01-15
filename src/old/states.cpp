#include "states.h"
#include "general.h"

ostream& operator<<(ostream &os, state_angular const &s){
  os << s.lr << " " << s.lrho << " " << s.lcoup << " " << s.lcom << " (" << s.nr << " " << s.nrho << ") | " 
    << setw(3) << s.index << " " << setw(3) << s.index_rel << " " << setw(3) << s.index_E << " | ";
  for (auto &a : s.allowed_transformations){
    os << a << setw(3) << " ";
  }
  return os;
}

void generate_angular_states(vector<state_angular> &tar, int lmax_r, int lmax_rho, int lmax_coup, int lmax_com, int ltotal, string parity, bool single_angular_state, bool use_n_powers, int max_quanta){
 tar = vector<state_angular>();
  if (single_angular_state){
    if (ltotal == 0) tar.push_back(state_angular(0,0,0,0,0,0,0,0,0));
    else if (ltotal == 1) tar.push_back(state_angular(0,0,0,1,1,0,0,0,0));
    else{
      cout << "ERROR: single_angular_state is not usable for ltotal >= 2" << endl;
      throw 0;
    }
  }
  else{
    int state_index = 0;
    int state_index_rel = 0;
    int nmax = use_n_powers ? max_quanta / 2 : 0;
    cout << "zzz: " << nmax << endl;
    for (int lr=0; lr <= lmax_r; lr++){
      for (int nr=0; nr <= nmax; nr++){
        for (int lrho=0; lrho <= lmax_rho; lrho++){
          for (int nrho=0; nrho <= nmax; nrho++){
            if (lr + 2*nr + lrho + 2*nrho > max_quanta) continue;
            for (int lcoup=0; lcoup <= lmax_coup; lcoup++){
              bool state_test_rel = false;
              for (int lcom=0; lcom <= lmax_com; lcom++){
                if (ltotal > lcoup + lcom ||
                    ltotal < abs(lcoup - lcom) ||
                    lcoup > lr + lrho ||
                    lcoup < abs(lr - lrho)) continue;

                int ltest = (lr + lrho + lcom + ltotal) % 2;
                if (parity == "natural" && ltest != 0) continue;
                if (parity == "unnatural" && ltest != 1) continue;


                state_test_rel = true;
                int state_index_E = lcoup + (lmax_coup+1) * (lrho + (lmax_rho+1) * lr);
                tar.push_back(state_angular(state_index++, state_index_rel, state_index_E, lr, lrho, lcoup, lcom, nr, nrho));
              }
              if (state_test_rel) state_index_rel++;
            }
          }
        }
      }
    }
  }
}


void state_angular::generate_allowed_transformations(vector<state_angular> &states_transformed){
  for (auto &st : states_transformed){
    if (st.lcom != lcom) continue;
    if (st.lcoup != lcoup) continue;
    if (lr + lrho + 2*nr + 2*nrho != st.lr + st.lrho + 2*st.nr + 2*st.nrho) continue;
    
    allowed_transformations.push_back(st.index);
  }
}

