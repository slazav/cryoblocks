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

If blocks have zero heat capacity they are calculated differently:
before and after each step a zero of total heat flow on every such block
is found as a function of temperaures of these blocks.

#### Dimention parameter reading

All parameters usually have some dimensions, it should be specified. For example,
time in `run` commabs can be written as `1s`, `0.1m`, `1e-2ms`, etc.

#### Blocks

Following types of blocks are supported:

* `block <name> <temperature> type=bath` -- A heat bath with intinite heat capacity and
constant temperature.

* `block <name> <temperature> type=zero-c` -- A block with zero heat capacity.

* `block <name> <temperature> type=simple C=<value>` -- A block with constant heat capacity.

* `block <name> <temperature> type=paramagnet [parameters]` -- paramagnetic material
(e.g. copper nuclei). Parameters:
  * `material=<v>` -- Set parameters `Bint`, `gyro`, `spin` for a certain material: `copper` (nuclei), `he3` (solid).
  * `Bint=<v>` -- Internal field (overrides value set by `material` parameter), default 0.
  * `gyro=<v>` -- Gyromagnetic ratio (overrides value set by `material` parameter).
  * `spin=<v>` -- Spin (overrides value set by `material` parameter).
  * `moles=<v>` -- Number of moles.
  * `mass=<v>` -- If material is set then mass can be used instead of moles.

* `block <name> <temperature> type=liquid_he3 [parameters]` -- liquid He3. Parameters:
  * `P=<v>` -- Pressure
  * `moles=<v>` -- Number of moles
  * `mass=<v>` -- Can be used instead of moles.
  * `volume=<v>` -- Can be used instead of moles (pressure-dependent molar volume is used).

#### Links

Following types of links are supported:

* `link <name> <block1> <block2> type=const` -- A constant heat transfer. Not a physical process,
but good for tests.

* `link <name> <block1> <block2> type=metal_bar R=<v>` -- A bar made of metal.
Total resistance R, Wiedemann-Franz low is used to calculate heat conductivity.

* `link <name> <block1> <block2> type=simple_bar [paramters]` -- A bar made of some material
with Ka*T^Kb heat conductivity [W/m/K]. Length L, cross-section area S. Parameters:
  * `material=<v>` -- torlon4203, GRP, nylon, G10-CR, macor, stycast1266, stycast2850ft,
     araldite_ct200, CuNi, Manganin. 
  * `Ka=<v>` -- Ka, overrides value set by `material` parameter
  * `Kb=<v>` -- Kb, overrides value set by `material` parameter
  * `S=<v>`  -- cross-section area of the bar
  * `L=<v>`  -- bar length

* `link <name> <block1> <block2> type=korringa [parameters]` --
Spin-lattice coupling (heat flow between nuclear spin system and electrons,
Korringa low). Parameters:

  * `material=<v>` -- predefined values for. Variants: copper
  * `Bint=<v>` -- internal field (overrides value set by `material` parameter), default 0
  * `gyro=<v>` -- gyromagnetic ratio (overrides value set by `material` parameter)
  * `spin=<v>` -- spin (overrides value set by `material` parameter)
  * `kappa=<v>` -- high-field value of Karringa constant (overrides value set by `material` parameter)
  * `alpha=<v>` -- parameter alpha in field dependence of Karringa constant (overrides value set by `material` parameter), default 1
  * `moles=<v>` -- number of moles
  * `mass=<v>` -- if material is set then mass can be used instead of moles

* `link <name> <block1> <block2> type=el-ph [parameters]` --
Electron-phonon coupling (Qdot = C*(T1^5-T2^5)*nmoles, see Pobell book f.10.9).
Parameters:

  * `material=<v>` -- predefined values for. Variants: copper
  * `C=<v>` -- set parameter C (overrides value set by `material` parameter), default 0
  * `moles=<v>` -- number of moles
  * `mass=<v>` -- if material is set then mass can be used instead of moles

* `link <name> <block1> <block2> type=kap_res_he3 area=<S> power=<N>` --
Kapitza resistance between He3 and solids. `R = 900/T/S [W]` for temperature `T` and area `S`.
Power should be 1, 2, or 3.

