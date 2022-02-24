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


class Calculator {
  std::map<std::string, std::shared_ptr<BlockBase> > blocks; // storage for blocks
  std::map<std::string, std::shared_ptr<LinkBase> >  links;  // links, block thermal connections
  std::map<std::string, double> temps; // block temperatures

  typedef std::vector<std::string>::const_iterator arg_cit;

  std::vector<std::string> print_list; // list of parameters to print

  // magnetic field and its sweeping rate [T, T/s]
  double B=0, Bdot=0;

  // current time [s]
  double t=0;

  public:

  // get temperature of a block
  double get_block_temp(const std::string & name){
    auto b = temps.find(name);
    if (b==temps.end()) throw Err() << "Accessing non-existing block: " << name;
    return b->second;
  }

  // get heat capacity of a block
  double get_block_heat_cap(const std::string & name){
    auto b = blocks.find(name);
    if (b==blocks.end()) throw Err() << "Accessing non-existing block: " << name;
    return b->second->get_heat_cap(get_block_temp(name), B);
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

    // extract type parameter
    auto type = get_key_val(b,e, "type");
    if (type == "") throw Err() << "block type is not set";

    // thermal bath, infinite heat capacity
    if (type == "bath") {
      auto opts = get_key_val_args(b,e, {"type=", "print="});
      blocks.emplace(name, new BlockSimple(INFINITY));
      return;
    }

    if (type == "simple") {
      auto opts = get_key_val_args(b,e, {"type=", "print=", "C=1J/K"});
      blocks.emplace(name, new BlockSimple(read_heat_cap(opts["C"])));
      return;
    }
    throw Err() << "unknown block type: " << type;


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

    // extract type
    auto type = get_key_val(b,e, "type");
    if (type == "")   throw Err() << "link type is not set";

    if (type=="const"){
      auto opts = get_key_val_args(b,e, {"type=", "Qdot=0W"});
      auto Q = read_power(opts["Qdot"]);
      links.emplace(name, new LinkConst(bl1, bl2, Q));
      return;
    }

    throw Err() << "unknown link type: " << type;
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
        if (blocks.count(n)==0) throw Err() << "Unknown block name for data output: " << n;
        std::cout << get_block_temp(n); 
        continue;
      }

      // Heat flow through a link: Q(<name>)
      if (v.size()>3 && v[0]=='Q' && v[1]=='(' && v[v.size()-1]==')'){
        auto n = v.substr(2,v.size()-3);
        if (links.count(n)==0) throw Err() << "Unknown link name for data output: " << n;
        auto b1n = links[n]->get_block1();
        auto b2n = links[n]->get_block2();
        auto qdot = links[n]->get_qdot(get_block_temp(b1n), get_block_temp(b2n));
        std::cout << qdot;
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
      // TODO: find and set temperatures of zero-C blocks

      // find heat flows to each block
      std::map<std::string, double> bq;
      for (const auto & l : links){
        auto b1n = l.second->get_block1();
        auto b2n = l.second->get_block2();
        auto qdot = l.second->get_qdot(get_block_temp(b1n), get_block_temp(b2n));
        bq[b1n] = ((bq.count(b1n) == 0)? 0.0 : bq[b1n]) - qdot;
        bq[b2n] = ((bq.count(b2n) == 0)? 0.0 : bq[b2n]) + qdot;
      }

      // find temperature change of each block affected by heat flows
      for (const auto & b : bq){
        auto C = get_block_heat_cap(b.first);
        temps[b.first] += b.second * dt / C;
      }
      print_data();

      // change magnetic field
      B = B + Bdot*dt;
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

  // program should be run with one argument - name of command file
  if (argc!=2){
    std::cerr << "Usage: " << argv[0] << " <command file>\n";
    return 1;
  }

  // open command file
  std::ifstream in_c(argv[1]);
  if (!in_c.good()) throw Err() << "Can't open command file: " << argv[1];

  // find prefix (command file name without extension)
  const char * pos1 = rindex(argv[1], '/');
  const char * pos2 = rindex(argv[1], '.');
  if (!pos1) pos1 = argv[1];
  auto pref = pos2 && pos2>pos1+1 ? std::string(argv[1], pos2-argv[1]) : argv[1];

  std::ofstream out_l((pref + ".run.log").c_str()); // log commands

  // Main cycle: read command file line by line
  while (1){
    // Read one line, detect EOF
    std::string line;
    if (!getline(in_c, line)){
      std::cout << "End of command file, stop calculations\n";
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
