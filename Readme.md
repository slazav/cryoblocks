## Cryoblocks -- thermal flow calculator for cryogenic (and any other) systems

The system is defined as a set of `blocks` and `links` with certain properties. Program read a command file
with the following structure:

* Comments start with `#`.

* Empty lines are skipped.

* `block <name> <temperature> [parameters]` -- Define a block with a given name and startig temperature.

* `link  <name> <block1> <block2> [parameters]` -- Define a named link which connects block1 and block2.

* `run  <time> <time step>` -- Do calculations for a given time, print table with temperatures and heat flows.

* `print <par1> ...` -- Set list of parameters for output during `run` command.
  Use `T(<name>)` for a block temperature, `Q(name)` for heat flow
  through a link, `t` for time, `B` for magnetic field. Example: `print t
  H T(block1) T(block2) Q(link1)`.

* `field <value>` -- Set magnetic field.

* `field_rate <value>` -- Set magnetic field rate.

* `exit` -- Stop processing file and exit.

#### Dimention parameter reading

All parameters usually have some dimensions, it should be specified. For example,
time in `run` commabs can be written as `1s`, `0.1m`, `1e-2ms`, etc.
