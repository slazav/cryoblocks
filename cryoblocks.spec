Name:         cryoblocks
Version:      1.1
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

%files
%_bindir/cryoblocks

%changelog
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

