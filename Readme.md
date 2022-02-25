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

#### Physics

Consider a block with temperature `T` with external heating power `dQ`.
Entropy of the block `S` depends on temperature and maybe other
parameters (magnetic field `B`).


Then `dQ = T dS = C dT + D dB`, where `C =  T \partial S / \partial T)` is
heat capacity, and `D = T \partial S / \partial B` - cooling power of
demagnetization refrigerator (or how to call it better?). Using this
formula each block object can calculate change of temperature `dT` as a
function of `dQ`, `T`, `B`, and `dB`.

A thermal link object can calculate heat transfer depending on
temperatures of two blocks connected by the link.

Calculation is done using a fixed time step (TODO: it would be good to
implement adaptive steps). On each step initial tempetures of all blocks
are known. For each link heat flows are calculated, Total power applied
to each block is found. Using these values new temperatures of each block
are found.

Currently blocks with zero heat capacity are not supported. To implement
this one should solve equation `Q(T) = 0` for such blocks at every step.


#### Dimention parameter reading

All parameters usually have some dimensions, it should be specified. For example,
time in `run` commabs can be written as `1s`, `0.1m`, `1e-2ms`, etc.

#### Blocks

Following types of blocks are supported:

* `block <name> <temperature> type=bath` -- A heat bath with intinite heat capacity and
constant temperature.

* `block <name> <temperature> type=bath C=<value>` -- A block with constant heat capacity.

* `block <name> <temperature> type=paramagnet [parameters]` -- paramagnetic material
(e.g. copper nuclei). Parameters:
  * `material=<v>` -- `copper` (nuclei), `he3` (solid)
  * `Bint=<v>` -- internal field (overrides value set by `material` parameter)
  * `gyro=<v>` -- gyromagnetic ratio (overrides value set by `material` parameter)
  * `spin=<v>` -- spin (overrides value set by `material` parameter)
  * `nmol=<v>` -- number of mols

#### Links

Following types of links are supported:

* `link <name> <block1> <block2> type=const` -- A constant heat transfer. Not a physical process,
but good for tests.

