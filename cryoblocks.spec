Name:         cryoblocks
Version:      1.0
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
* Tue Apr 12 2022 Vladislav Zavjalov <slazav@altlinux.org> 1.0-alt1
- First version

