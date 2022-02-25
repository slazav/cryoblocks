#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <cmath>
#include <cstring>

#include "c_err.h"
#include "c_args.h"
#include "c_blocks.h"
#include "c_links.h"
#include "c_dim.h"

#include <gsl/gsl_linalg.h>

class Calculator {
  std::map<std::string, std::shared_ptr<BlockBase> > blocks; // storage for blocks
  std::map<std::string, std::shared_ptr<LinkBase> >  links;  // links, block thermal connections

  std::map<std::string, double> temps; // block temperatures
  std::map<std::string, std::pair<std::string, std::string> > conn; // connections link -> block1,block2

  typedef std::vector<std::string>::const_iterator arg_cit;

  std::vector<std::string> print_list; // list of parameters to print

  // magnetic field and its sweeping rate [T, T/s]
  double B=0, Bdot=0;

  // current time [s]
  double t=0;


  public:
  /***********************************/

  // get temperature of a block
  double get_block_temp(const std::string & name) const {
    auto b = temps.find(name);
    if (b==temps.end()) throw Err() << "Unknown block: " << name;
    return b->second;
  }

  // get heat flow through a link
  double get_link_flow(const std::string & name) const {
    auto l = links.find(name);
    auto c = conn.find(name);
    if (l==links.end() || c==conn.end())
      throw Err() << "Unknown link: " << name;
    auto T1 = get_block_temp(c->second.first);
    auto T2 = get_block_temp(c->second.second);
    return l->second->get_qdot(T1,T2);
  }

  void set_magn_field(const double v) {B = v;}
  void set_magn_field_rate(const double v) {Bdot = v;}

  /***********************************/
  // Add a block
  void add_block(const std::string & name, const double temp,
                 const arg_cit & b, const arg_cit & e){
    if (blocks.count(name)>0)
      std::cout << "# replacing existing block\n";

    temps.emplace(name, temp);
    blocks.emplace(name, create_block(b,e));
  }

  /***********************************/
  // Add a link
  void add_link(const std::string & name,
                const std::string & bl1, const std::string & bl2,
                const arg_cit & b, const arg_cit & e){
    if (blocks.count(bl1) == 0) throw Err() << "no such block: " << bl1;
    if (blocks.count(bl2) == 0) throw Err() << "no such block: " << bl2;
    if (links.count(name)>0)
      std::cout << "# replacing existing link\n";

    links.emplace(name, create_link(b,e));
    conn.emplace(name, std::make_pair(bl1, bl2));
  }

  /***********************************/
  // Set list of output values
  void set_print_list(const std::vector<std::string> list){
    print_list = list;
  }

  // print heat flows and temperatures according to print_list
  void print_data() {
    size_t n = 0;
    for (const auto v:print_list){
      if (n!=0) std::cout << "\t"; n++;

      // Temperature of a block: T(<name>)
      if (v.size()>3 && v[0]=='T' && v[1]=='(' && v[v.size()-1]==')'){
        auto n = v.substr(2,v.size()-3);
        // extra check to have more understandable error message:
        if (blocks.count(n)==0) throw Err() << "Unknown block name for data output: " << n;
        std::cout << get_block_temp(n); 
        continue;
      }

      // Heat flow through a link: Q(<name>)
      if (v.size()>3 && v[0]=='Q' && v[1]=='(' && v[v.size()-1]==')'){
        auto n = v.substr(2,v.size()-3);
        // extra check to have more understandable error message:
        if (links.count(n)==0) throw Err() << "Unknown link name for data output: " << n;
        std::cout << get_link_flow(n);
        continue;
      }

      // Magnetic field
      if (v == "B"){ std::cout << B; continue; }

      // Time
      if (v == "t"){ std::cout << t; continue;
      }
    }
    std::cout << "\n";
  }

  /***********************************/
  // Run a calculation for some time, print results
  void run(double te, double dt){


    // print table header
    std::cout << "#" << print_list << "\n";

    // time cycle
    te+=t;
    for (; t<=te; t+=dt){

      // find and set temperatures of zero-C blocks
      double max = 1e-3;
      for (int i=0; i<100; i++){
        max = adjust_zero_c_temps(max);
        if (max < 1e-6) break;
        if (i==99) throw Err() << "Can't find temperatures of zero-C blocks";
      }

      // find heat flows to each block
      std::map<std::string, double> bq;
      for (const auto & l : links){
        auto n = l.first;
        auto b1n = conn[n].first;
        auto b2n = conn[n].second;
        auto qdot = get_link_flow(n);
        bq[b1n] = ((bq.count(b1n) == 0)? 0.0 : bq[b1n]) - qdot;
        bq[b2n] = ((bq.count(b2n) == 0)? 0.0 : bq[b2n]) + qdot;
      }

      // find temperature change of each block (except zero-c)
      for (const auto & b : blocks){
        if (b.second->is_zero_c()) continue;
        auto n = b.first;
        auto dQ = bq.count(n)? bq[n] * dt : 0;
        auto T  = get_block_temp(n);
        temps[n] += b.second->get_dt(dQ, T, B, Bdot*dt);
      }

      // repeat calculation of zero-C blocks to have then in the equilibrium after the step
      max = 1e-3;
      for (int i=0; i<100; i++){
        max = adjust_zero_c_temps(max);
        if (max < 1e-6) break;
        if (i==99) throw Err() << "Can't find temperatures of zero-C blocks";
      }

      print_data();

      // change magnetic field
      B = B + Bdot*dt;
    }
  }

  // Do one step of zero-C temperature calculation.
  // Return maximum of relative temperature change dT/T.
  // Parameter rstep is a relative temperature step (dT/T) for
  // estimating derivatives.
  double adjust_zero_c_temps(double rstep) {

    // count zero-c blocks
    std::map<std::string, size_t> num; // name -> number
    size_t nzblocks = 0;
    for (const auto & b : blocks)
      if (b.second->is_zero_c()) num[b.first] = nzblocks++;

    // Q - Heat flow to each zero-C block;
    // dQdT - Derivative dQi/dTj for i,j in zero-c blocks.
    gsl_matrix *dQdT = gsl_matrix_alloc(nzblocks, nzblocks);
    gsl_vector *Q    = gsl_vector_alloc(nzblocks);
    gsl_vector_set_all(Q, 0.0);
    gsl_matrix_set_all(dQdT, 0.0);

    // Find all links which connect zero-c blocks
    // Fill Q and dQdT arrays
    for (const auto & l : links){
      auto n = l.first;
      auto b1n = conn[n].first;
      auto b2n = conn[n].second;
      auto b1 = blocks[b1n];
      auto b2 = blocks[b2n];
      if (!b1->is_zero_c() && !b2->is_zero_c()) continue;

      // temperatures of both blocks and flow between them:
      auto T1 = temps[b1n];
      auto T2 = temps[b2n];
      auto qdot = l.second->get_qdot(T1,T2);

      // If block1 is zero-c block, add to
      // Q(block1), dQ(block1)/dT1, dQ(block2)/dT1
      if (b1->is_zero_c()){
        auto T1a = T1*(1+rstep);
        auto qdota = l.second->get_qdot(T1a,T2);
        auto dqdt1 = (qdota-qdot)/(T1a-T1);
        auto i = num[b1n];
        *gsl_vector_ptr(Q, i) -= -qdot;
        *gsl_matrix_ptr(dQdT,i,i) -= dqdt1;
        if (b2->is_zero_c()){
          auto j = num[b2n];
          *gsl_matrix_ptr(dQdT,j,i) += dqdt1;
        }
      }
      // Same for block2
      if (b2->is_zero_c()){
        auto T2a = T2*(1+rstep);
        auto qdota = l.second->get_qdot(T1,T2a);
        auto dqdt2 = (qdota-qdot)/(T2a-T2);
        auto i = num[b2n];
        *gsl_vector_ptr(Q,i) += -qdot;
        *gsl_matrix_ptr(dQdT,i,i) += dqdt2;
        if (b1->is_zero_c()){
          auto j = num[b1n];
          *gsl_matrix_ptr(dQdT,j,i) -= dqdt2;
        }
      }
    }

    // Newton method: solve linear system  dQ_j/dTj * DTj = -Q_i

    int s;
    gsl_vector *DT = gsl_vector_alloc(nzblocks);
    gsl_permutation *p = gsl_permutation_alloc (nzblocks);
    gsl_linalg_LU_decomp (dQdT, p, &s);
    gsl_linalg_LU_solve (dQdT, p, Q, DT);
    gsl_permutation_free(p);
    gsl_vector_free(Q);
    gsl_matrix_free(dQdT);

    // Adjust temperatures of zero-c blocks
    double max = 0.0;
    for (const auto & n : num){
      auto v = gsl_vector_get(DT, n.second);
      temps[n.first] += v;
      double e = fabs(v/temps[n.first]);
      if (max < e) max = e;
    }
    gsl_vector_free(DT);

    return max;
  }

};

/********************************************************************/
// Main program block: read commands from a file and send them
// to the Calculator class

int
main(int argc, char *argv[]){
try{

  Calculator calc;

  // program should be run with one argument - name of command file
  if (argc!=2){
    std::cerr << "Usage: " << argv[0] << " <command file>\n";
    return 1;
  }

  // open command file
  std::ifstream in_c(argv[1]);
  if (!in_c.good()) throw Err() << "Can't open command file: " << argv[1];

  // find prefix (command file name without extension)
  //const char * pos1 = rindex(argv[1], '/');
  //const char * pos2 = rindex(argv[1], '.');
  //if (!pos1) pos1 = argv[1];
  //auto pref = pos2 && pos2>pos1+1 ? std::string(argv[1], pos2-argv[1]) : argv[1];

  // Main cycle: read command file line by line
  while (1){
    // Read one line, detect EOF
    std::string line;
    if (!getline(in_c, line)){
      std::cout << "# end of command file\n";
      break;
    }

    // Remove comments:
    size_t c = line.find('#');
    if (c != std::string::npos)
      line = line.substr(0, c);

    // Split the line into command and arguments:
    std::istringstream in(line);
    std::string cmd;
    in >> cmd;
    if (!in) continue; // empty command

    std::vector<std::string> args;
    while (1) {
      std::string a;
      in >> a;
      if (!in) break;
      args.push_back(a);
    }
    int narg = args.size();

    // logging to cout
    std::cout << "# " << line << "\n";

    /*****************************************/
    // Process commands

    // Exit
    if (cmd == "exit") { break; }

    // Define a block
    if (cmd == "block") {
      if (args.size() < 2)
        throw Err() << "Not enough arguments, expect: block <name> <temp> [options]";
      calc.add_block(args[0], read_temp(args[1]), args.begin()+2, args.end());
      continue;
    }

    // Define a link
    if (cmd == "link") {
      if (args.size() < 3)
        throw Err() << "Not enough arguments, expect: link <name> <block1> <block2> [options]";
      calc.add_link(args[0], args[1], args[2], args.begin()+3, args.end());
      continue;
    }

    // Wait for some time
    if (cmd == "run") {
      if (args.size() != 2)
        throw Err() << "Wrong number of arguments. Expect: run <time> <time step>";
      calc.run(read_time(args[0]), read_time(args[1]));
      continue;
    }

    // Set list of output values
    if (cmd == "print") {
      calc.set_print_list(args);
      continue;
    }

    // Set magnetic field
    if (cmd == "field") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: field <B>";
      calc.set_magn_field(read_magn_field(args[0]));
      continue;
    }

    // Set magnetic field sweep rate
    if (cmd == "field_rate") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: field_rate <dB/dt>";
      calc.set_magn_field_rate(read_magn_field_rate(args[0]));
      continue;
    }

    /*****************************************/
    throw Err() << "Unknown command: " << cmd;
  }
  return 0;

}
catch (const Err & e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
