#include <string>
#include <vector>

#include <iostream>

#include "inc/err.h"
#include "inc/assert_err.h"
#include "c_args.h"
#include "c_blocks.h"
#include "c_links.h"
#include "c_dim.h"


/********************************************************************/

int
main(int argc, char *argv[]) {
try{

  { // simple/zero-c/bath blocks
    auto b1 = create_block({"type=simple", "C=1mJ/K"}); // simple block with 1mJ/K heat capacity
    assert_eq(b1->is_zero_c(), false);
    assert_feq(b1->get_dt(1e-6, 1, 2, 3), 1e-3, 1e-12); // Apply 1uW, temperature change is 1e-3K
    assert_feq(b1->get_dt(1e-9, 3, 4, 5), 1e-6, 1e-12);

    auto b2 = create_block({"type=simple", "C=0J/K"}); // zero heat capacity in the simple block
    assert_eq(b2->is_zero_c(), true);

    auto b3 = create_block({"type=simple"}); // default heat capacity is zero
    assert_eq(b3->is_zero_c(), true);

    auto b4 = create_block({"type=zero-c"}); // zero-c block
    assert_eq(b4->is_zero_c(), true);

    auto b5 = create_block({"type=bath"});
    assert_eq(b5->is_zero_c(), false);
    assert_feq(b5->get_dt(1.1, 2, 3, 4), 0, 1e-12); // zero change of temperature
  }

  {
    auto l1 = create_link({"type=const", "Qdot=1mW"}); // link with a constant heat flow
    assert_feq(l1->get_qdot(1.1,2.2,3.3,0), 1e-3, 1e-12); // always 1mW

    auto l2 = create_link({"type=const", "Qdot=-1mW"}); // link with a constant heat flow
    assert_feq(l2->get_qdot(1.1,2.2,3.3,0), -1e-3, 1e-12); // always -1mW

    auto l3 = create_link({"type=const"}); // link with a constant heat flow
    assert_feq(l3->get_qdot(1.1,2.2,3.3,0), 0, 1e-12); // always zero
  }

  {
    auto l1 = create_link({"type=simple_bar", "Ka=2", "Kb=2", "S=1mm^2", "L=1m"}); // K = 2*T^2
    // Q = S*K*dT/dL  -->  Q = S/L int(K dT) = S/L 2/3 (T1^3-T2^3)
    assert_feq(l1->get_qdot(0.1,0.05, 0,0), 5.83333e-10, 0.00001e-10);

    auto l2 = create_link({"type=simple_bar", "material=CuNi", "S=1mm^2", "L=1cm"});
    assert_feq(l2->get_qdot(0.1,0.05, 0,0), 1.88514e-8, 0.00001e-8);

    auto l3 = create_link({"type=metal_bar", "R=1mOhm"});
    assert_feq(l3->get_qdot(0.1,0.05, 0,0), 9.15e-8, 0.0001e-8);

    auto l4 = create_link({"type=metal_bar", "R=1uOhm"});
    assert_feq(l4->get_qdot(0.01,0.0, 0,0), 1.22e-6, 0.0001e-6);

    auto l5 = create_link({"type=korringa", "material=copper", "moles=1"});
    assert_feq(l5->get_qdot(0.002,0.001,0.01,0), 1.282e-10, 0.001e-10);

    auto l6 = create_link({"type=el_ph", "material=copper", "moles=1"});
    assert_feq(l6->get_qdot(0.002,0.001,0.01,0), 4.408e-10, 0.001e-10);

    auto l7 = create_link({"type=kap_res_he3", "area=2m^2"});
    assert_feq(l7->get_qdot(0.0021,0.002,0.01,0), 4.555e-10, 0.001e-10);
  }

  return 0;

}
catch (const Err & e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
