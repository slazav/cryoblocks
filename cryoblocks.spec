Name:         cryoblocks
Version:      1.3
Release:      alt1

Summary:      Cryoblocks -- thermal flow calculator for cryogenic (and any other) systems
Group:        System
URL:          https://github.com/slazav/cryoblocks
License:      GPL

Packager:     Vladislav Zavjalov <slazav@altlinux.org>

Source:       %name-%version.tar
BuildRequires: gcc-c++
BuildRequires: libgsl-devel

%description
Cryoblocks -- thermal flow calculator for cryogenic (and any other) systems.

%prep
%setup -q

%build
%make

%install
install -m 755 -D cryoblocks %buildroot%_bindir/cryoblocks
install -m 755 -D cryoblocks_fit %buildroot%_bindir/cryoblocks_fit

%files
%_bindir/cryoblocks

%changelog
* Wed Mar 22 2023 Vladislav Zavjalov <slazav@altlinux.org> 1.3-alt1
- v1.3
  - fix error: update temperature when replacing a block
  - fix error in paramagnet heat capacity
  - add Curie-Weiss magnet
  - add Bfactor parameter for magnetic blocks
  - examples: update for new "run" command interface
  - in run command reduce the last step if needed to reach exact final time
  - link type=field_heat_leak: dB2 parameter, heat leak proportional to (dB/dt)^2

* Thu May 05 2022 Vladislav Zavjalov <slazav@altlinux.org> 1.2-alt1
- v1.2:
 - Change in `run` command: no second argument, options "abs",
   "to_field", "step". Also add `run_to` command which is now depreciated.
 - Improve calculation of a mixture of big and small heat capacity blocks..
   (dump tiny temperature oscillations).
 - New commands: "read", "read_factors", "time_shift", "delete", "include"
 - Allow reading commands from stdin (use "-" as filename)
 - Add temperature dependence to simple blocks and links. Now heat capacity
   of a simple block and thermal conductivity of a simple link  can be set
   as <factor>*T^<power>.
 - New link type: field_heat_leak: heat flow proportional to
   B, B^2, (dB/dt)^2
 - kap_res_he3 link: add C parameter to override resistance scaling.
 - Add `cryoblocks_fit` - a python script for fitting data to a cryoblocks
   model. All fit parameters are defined in cryoblocks command file,
   same file can be passed to cryoblocks and cryoblocks_fit programs.

* Sat Apr 16 2022 Vladislav Zavjalov <slazav@altlinux.org> 1.1-alt1
- v1.1:
  - adaptive steps with configurable accuracy;
  - printing data to a file (print_to_file <name>);
  - basic variables (define <name> <value>);
  - support for blocks with zero heat capacity;
  - block types: bath, zero-c, simple, paramagnet, liquid_he3;
  - link types: const, simple, simple_bar, metal_bar,
    korringa, el_ph, kap_res_he3, dilution_cooling, dilution_circ;

* Tue Apr 12 2022 Vladislav Zavjalov <slazav@altlinux.org> 1.0-alt1
- First version

