#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <memory>
#include <map>
#include <cmath>
#include <cstring>

#include "inc/err.h"
#include "c_args.h"
#include "c_blocks.h"
#include "c_links.h"
#include "c_dim.h"

#include <gsl/gsl_linalg.h>

void my_gsl_handler (const char * reason, const char * file, int line, int gsl_errno){ return; }

class Calculator {
  std::map<std::string, std::shared_ptr<BlockBase> > blocks; // storage for blocks
  std::map<std::string, std::shared_ptr<LinkBase> >  links;  // links, block thermal connections

  std::map<std::string, double> temps0; // block temperatures
  std::map<std::string, std::pair<std::string, std::string> > conn; // connections link -> block1,block2

  std::map<std::string, std::string> defs; // parameter definitions, name -> value

  typedef std::vector<std::string>::const_iterator arg_cit;
  typedef std::map<std::string, double> d_blpars;

  std::vector<std::string> print_list; // list of parameters to print

  std::ostream * out = &std::cout; // pointer to output stream
  std::shared_ptr<std::ostream> out_f; // storage for std::ofstream

  // magnetic field and its sweeping rate [T, T/s]
  double B=0, Bdot=0;

  // current time [s]
  double t=0;

  // Parameters for adaptive steps
  double max_tempstep = 1e-2; // max relative temperature change on each calculation step
  double max_tempacc  = 1e-6; // max relative accuracy on each step

  bool print_substeps = false; // print all calculation steps instead of using

  public:
  /***********************************/

  // get block (shared_ptr), throw error if block does not exist
  std::shared_ptr<BlockBase> get_block(const std::string & name) const {
    auto b = blocks.find(name);
    if (b==blocks.end()) throw Err() << "Unknown block: " << name;
    return b->second;
  }

  // get link (shared_ptr), throw error if link does not exist
  std::shared_ptr<LinkBase> get_link(const std::string & name) const {
    auto l = links.find(name);
    if (l==links.end()) throw Err() << "Unknown link: " << name;
    return l->second;
  }

  // get connection (name1-name2 pair), throw error if link does not exist
  std::pair<std::string, std::string> get_conn(const std::string & name) const {
    auto c = conn.find(name);
    if (c==conn.end()) throw Err() << "Unknown link: " << name;
    return c->second;
  }

  // get temperature of a block
  double get_block_temp(const d_blpars & temps, const std::string & name) const {
    auto b = temps.find(name);
    if (b==temps.end()) throw Err() << "Unknown block: " << name;
    return b->second;
  }

  /***********************************/

  // get heat flow through a link
  double get_link_flow(const d_blpars & temps, const std::string & name) const {
    auto l = get_link(name);
    auto c = get_conn(name);
    auto T1 = get_block_temp(temps, c.first);
    auto T2 = get_block_temp(temps, c.second);
    return l->get_qdot(T1,T2,B);
  }

  void set_magn_field(const double v) {B = v;}
  void set_magn_field_rate(const double v) {Bdot = v;}
  void set_print_substeps(const bool v) {print_substeps = v;}
  void set_max_tempstep(const double v) {max_tempstep = v;}
  void set_max_tempacc(const double v) {max_tempacc = v;}

  void set_print_to_file(const std::string & name){
    if (name == "-") {
      out_f.reset();
      out = & std::cout;
    } else {
      out_f = std::shared_ptr<std::ostream>(new std::ofstream(name));
      out = out_f.get();
    }
    if (!*out) throw Err() << "Can't open file: " << name;
  }

  /***********************************/
  // Define parameters
  void define(const std::string & name, const std::string & value) {
    defs.emplace(name, value);
  }
  void apply_defs(std::string & str){
    while (1) {
      auto n1 = str.find("${");
      if (n1 == std::string::npos) break;
      auto n2 = str.find("}", n1);
      auto v = defs.find(str.substr(n1+2, n2-n1-2));
      if (v == defs.end()) throw Err() << "Unknown parameter: " << str.substr(n1, n2-n1+1);
      str.replace(n1,n2-n1+1, v->second);
    }
  }

  /***********************************/
  // Add a block, return 1 if block was replaced
  int add_block(const std::string & name, const double temp,
                 const arg_cit & b, const arg_cit & e){
    int ret = 0;
    if (blocks.count(name)>0){
      blocks.erase(name);
      ret = 1;
    }
    temps0.emplace(name, temp);
    blocks.emplace(name, create_block(b,e));
    return ret;
  }

  /***********************************/
  // Add a link, return 1 if link was replaced
  bool add_link(const std::string & name,
                const std::string & bl1, const std::string & bl2,
                const arg_cit & b, const arg_cit & e){
    int ret = 0;
    if (blocks.count(bl1) == 0) throw Err() << "no such block: " << bl1;
    if (blocks.count(bl2) == 0) throw Err() << "no such block: " << bl2;
    if (links.count(name)>0){
      links.erase(name);
      ret = 1;
    }
    links.emplace(name, create_link(b,e));
    conn.emplace(name, std::make_pair(bl1, bl2));
    return ret;
  }

  /***********************************/
  // Set list of output values
  void set_print_list(const std::vector<std::string> list){
    print_list = list;
  }

  // print heat flows and temperatures according to print_list
  void print_data(const double t, const double B, const d_blpars & temps) const {
    size_t n = 0;
    for (const auto v:print_list){
      if (n!=0) *out << "\t"; n++;

      // Temperature of a block: T(<name>)
      if (v.size()>3 && v[0]=='T' && v[1]=='(' && v[v.size()-1]==')'){
        auto n = v.substr(2,v.size()-3);
        // extra check to have more understandable error message:
        if (blocks.count(n)==0) throw Err() << "Unknown block name for data output: " << n;
        *out << get_block_temp(temps, n); 
        continue;
      }

      // Heat flow through a link: Q(<name>)
      if (v.size()>3 && v[0]=='Q' && v[1]=='(' && v[v.size()-1]==')'){
        auto n = v.substr(2,v.size()-3);
        // extra check to have more understandable error message:
        if (links.count(n)==0) throw Err() << "Unknown link name for data output: " << n;
        *out << get_link_flow(temps, n);
        continue;
      }

      // Magnetic field
      if (v == "B"){ *out << B; continue; }

      // Time
      if (v == "t"){ *out << t; continue;
      }
    }
    if (print_list.size()>0) *out << "\n";
  }

  /***********************************/
  // Do one step of zero-C temperature calculation.
  // Arguments:
  //   temps - input: block temperatures
  //           return: modified temperatures
  //   rstep - relative temperature step (dT/T) for estimating derivatives
  //
  // Return:
  //   maximum of relative temperature change dT/T
  //
  double do_zeroc_step(d_blpars & temps, const double rstep) const {

    // count zero-c blocks
    std::map<std::string, size_t> num; // name -> number
    size_t nzblocks = 0;
    for (const auto & b : blocks)
      if (b.second->is_zero_c()) num[b.first] = nzblocks++;
    if (nzblocks==0) return 0;

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
      auto c = get_conn(n);
      auto b1 = get_block(c.first);
      auto b2 = get_block(c.second);
      if (!b1->is_zero_c() && !b2->is_zero_c()) continue;

      // temperatures of both blocks and flow between them:
      auto T1 = get_block_temp(temps, c.first);;
      auto T2 = get_block_temp(temps, c.second);;
      auto qdot = l.second->get_qdot(T1,T2,B);

      // If block1 is zero-c block, add to
      // Q(block1), dQ(block1)/dT1, dQ(block2)/dT1
      if (b1->is_zero_c()){
        auto T1a = T1*(1+rstep);
        auto qdota = l.second->get_qdot(T1a,T2,B);
        auto dqdt1 = (qdota-qdot)/(T1a-T1);
        auto i = num[c.first];
        *gsl_vector_ptr(Q, i) -= -qdot;
        *gsl_matrix_ptr(dQdT,i,i) -= dqdt1;
        if (b2->is_zero_c()){
          auto j = num[c.second];
          *gsl_matrix_ptr(dQdT,j,i) += dqdt1;
        }
      }
      // Same for block2
      if (b2->is_zero_c()){
        auto T2a = T2*(1+rstep);
        auto qdota = l.second->get_qdot(T1,T2a,B);
        auto dqdt2 = (qdota-qdot)/(T2a-T2);
        auto i = num[c.second];
        *gsl_vector_ptr(Q,i) += -qdot;
        *gsl_matrix_ptr(dQdT,i,i) += dqdt2;
        if (b1->is_zero_c()){
          auto j = num[c.first];
          *gsl_matrix_ptr(dQdT,j,i) -= dqdt2;
        }
      }
    }

    // Newton method: solve linear system  dQ_j/dTj * DTj = -Q_i

    int s;
    gsl_vector *DT = gsl_vector_alloc(nzblocks);
    gsl_permutation *p = gsl_permutation_alloc (nzblocks);
    try {
      int status = gsl_linalg_LU_decomp (dQdT, p, &s);
      if (status) throw Err() << "Error while solving zero-c blocks: " << gsl_strerror(status);
      status = gsl_linalg_LU_solve (dQdT, p, Q, DT);
      if (status) throw Err() << "Error while solving zero-c blocks: " << gsl_strerror(status);
    }
    catch (Err & e){
      e << "\nMatrix:\n";
      for (const auto & n : num){
        e << n.first << ": ";
        e << "\tQ=" << gsl_vector_get(Q,n.second);
        e << "\tT=" << get_block_temp(temps, n.first);
        e << "\tdQ/dT_j=";
        for (int i=0; i<nzblocks; i++) e << " " << gsl_matrix_get(dQdT, n.second, i);
        e << "\n";
      }
      gsl_permutation_free(p);
      gsl_vector_free(Q);
      gsl_matrix_free(dQdT);
      throw;
    }
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

  /***********************************/
  // Do zero-c calculation, modify array of temperatures
  void do_zeroc_calc(d_blpars & temps) const {
    // find and set temperatures of zero-C blocks
    double max = 1e-3;
    for (int i=0; i<100; i++){
      max = do_zeroc_step(temps, max);
      if (max < 1e-6) break;
      if (i==99) throw Err() << "Can't find temperatures of zero-C blocks";
    }
  }

  /***********************************/
  // Calculate a single step t+dt, B+dB
  // Modify temps array to return new values.
  // Be sure that do_zeroc_calc(temps) was called before.
  void do_step(
        const double t, const double dt,
        const double B, const double dB,
        std::map<std::string, double> & temps) const {

    // find heat flows to each block
    std::map<std::string, double> bq;
    for (const auto & l : links){
      auto nl = l.first;
      auto c = get_conn(nl);
      auto qdot = get_link_flow(temps, nl);
      auto n1 = c.first;
      auto n2 = c.second;
      bq[n1] = ((bq.count(n1) == 0)? 0.0 : bq[n1]) - qdot;
      bq[n2] = ((bq.count(n2) == 0)? 0.0 : bq[n2]) + qdot;
    }

    // find temperature change of each block (except zero-c)
    for (const auto & b : blocks){
      if (b.second->is_zero_c()) continue;
      auto n = b.first;
      auto dQ = bq.count(n)? bq[n] * dt : 0;
      auto T  = get_block_temp(temps, n);
      temps[n] += b.second->get_dt(dQ, T, B, Bdot*dt);
    }

    // repeat calculation of zero-C blocks to have them in the equilibrium after the step
    do_zeroc_calc(temps);
  }

  // Calculate a step t+dt, B+dB.
  // Split it into smaller parts if needed to keep accuracy
  void do_adaptive_step(
        const double t, const double dt,
        const double B, const double dB,
        std::map<std::string, double> & temps) const {
    std::list<double> brkpts;
    brkpts.push_back(t);
    brkpts.push_back(t+dt);
    auto t1=brkpts.begin();
    auto t2=t1; t2++;

    do_zeroc_calc(temps);

    while (1) {
      auto temps1 = temps, temps2=temps;
      auto dt1 = *t2-*t1;
      // dt1 step
      do_step(*t1, dt1, B+(*t1-t)*dB/dt, dt1*dB/dt, temps1);

      // two dt1/2 steps
      do_step(*t1, dt1/2, B+(*t1-t)*dB/dt, dt1/2*dB/dt, temps2);
      do_step(*t1+dt1/2, dt1/2, B+(*t1-t+dt1/2)*dB/dt, dt1/2*dB/dt, temps2);
      // Calculate calculation accuracy:
      // - relative temperature change in temps1,
      // - relative temperature difference between temps2 and temps1
      double m1=0,m2=0;
      for (auto i0 = temps.begin(), i1 = temps1.begin(), i2=temps2.begin();
           i0!=temps.end() && i1!=temps1.end() && i2!=temps2.end(); i0++,i1++,i2++){
        double v1 = fabs((i1->second - i0->second)/i0->second);
        double v2 = fabs((i1->second - i2->second)/i0->second);
        if (m1 < v1) m1 = v1;
        if (m2 < v2) m2 = v2;
      }
      // Split calculation interval, throw error if it's too small
      if (m1>max_tempstep || m2>max_tempacc){
        t2 = brkpts.insert(t2, *t1+dt1/2);
        if (*t1 + dt1/4.0 <= *t1) throw Err() << "timestep limit reached: " << dt1/2;
        continue;
      }
      // Finish or move to the next interval
      else {
        temps.swap(temps1);
        t1 = t2; t2++; // move to the next interval
        if (t2 == brkpts.end()) break;
        // print intermediate data if needed
        if (print_substeps) print_data(*t1, B+(*t1-t)*dB/dt, temps);
      }
    }
  }

  /***********************************/
  // Run a calculation for some time, print results
  void run(double te, double dt){

    // print table header
    *out << "#" << print_list << "\n";

    // time cycle
    te+=t;
    auto dB = Bdot*dt;
    do_zeroc_calc(temps0);
    print_data(t, B, temps0);
    for (; t<te; t+=dt, B+=dB){
      do_adaptive_step(t,dt,B,dB,temps0);
      print_data(t+dt, B+dB, temps0);
    }
  }


};

/********************************************************************/
// Main program block: read commands from a file and send them
// to the Calculator class

int
main(int argc, char *argv[]){
try{

  Calculator calc;
  gsl_set_error_handler (&my_gsl_handler);

  bool print_cmd = 1;                  // show commands on output

  // program should be run with one argument - name of command file
  if (argc!=2){
    std::cerr << "Usage: " << argv[0] << " <command file>\n";
    return 1;
  }

  // open command file
  std::istream * inf = &std::cin;
  std::shared_ptr<std::ifstream> inff;
  if (strcmp(argv[1],"-")!=0) {
    inff.reset(new std::ifstream(argv[1]));
    inf = inff.get();
  }
  if (!*inf) throw Err() << "Can't open command file: " << argv[1];

  // find prefix (command file name without extension)
  //const char * pos1 = rindex(argv[1], '/');
  //const char * pos2 = rindex(argv[1], '.');
  //if (!pos1) pos1 = argv[1];
  //auto pref = pos2 && pos2>pos1+1 ? std::string(argv[1], pos2-argv[1]) : argv[1];

  // Main cycle: read command file line by line
  while (1){
    // Read one line, detect EOF
    std::string line;
    getline(*inf, line);
    if (!*inf){
      if (print_cmd) std::cerr << "# end of command file\n";
      break;
    }

    // Remove comments:
    size_t c = line.find('#');
    if (c != std::string::npos)
      line = line.substr(0, c);

    // apply definitions
    calc.apply_defs(line);

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

    // logging command
    if (print_cmd) std::cerr << "# " << line << "\n";

    /*****************************************/
    // Process commands

    // Exit
    if (cmd == "exit") { break; }

    // Define a block
    if (cmd == "block") {
      if (args.size() < 2)
        throw Err() << "Not enough arguments, expect: block <name> <temp> [options]";
      auto ret = calc.add_block(args[0], read_value(args[1], "K"), args.begin()+2, args.end());
      if (ret && print_cmd) std::cerr << "# replacing existing block: " << args[0] << "\n";
      continue;
    }

    // Define a link
    if (cmd == "link") {
      if (args.size() < 3)
        throw Err() << "Not enough arguments, expect: link <name> <block1> <block2> [options]";
      auto ret = calc.add_link(args[0], args[1], args[2], args.begin()+3, args.end());
      if (ret && print_cmd) std::cerr << "# replacing existing link: " << args[0] << "\n";
      continue;
    }

    // Wait for some time
    if (cmd == "run") {
      if (args.size() != 2)
        throw Err() << "Wrong number of arguments. Expect: run <time> <time step>";
      calc.run(read_value(args[0], "s"), read_value(args[1], "s"));
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
      calc.set_magn_field(read_value(args[0], "T"));
      continue;
    }

    // Set magnetic field sweep rate
    if (cmd == "field_rate") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: field_rate <dB/dt>";
      calc.set_magn_field_rate(read_value(args[0], "T/s"));
      continue;
    }

    // Set print_substeps parameter
    if (cmd == "print_substeps") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: print_substeps 0|1";
      calc.set_print_substeps(atoi(args[0].c_str()));
      continue;
    }

    // Set max_tempstep parameter
    if (cmd == "max_tempstep") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: max_tempstep <value>";
      calc.set_max_tempstep(read_value(args[0], ""));
      continue;
    }

    // Set max_tempacc parameter
    if (cmd == "max_tempacc") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: max_tempacc <value>";
      calc.set_max_tempacc(read_value(args[0], ""));
      continue;
    }

    // Define a parameter
    if (cmd == "define") {
      if (args.size() != 2)
        throw Err() << "Wrong number of arguments. Expect: define <name> <value>";
      calc.define(args[0], args[1]);
      continue;
    }

    // Set output file
    if (cmd == "print_to_file") {
      if (args.size() != 1)
        throw Err() << "Wrong number of arguments. Expect: print_to_file <name>";
      calc.set_print_to_file(args[0]);
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
