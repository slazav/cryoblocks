## Cryoblocks -- thermal flow calculator for cryogenic (and any other) systems

### Basic example:

```
#!../cryoblocks

# Define initial temperature of our system:
define T0 20mK

# Describe copper block with mass 100g
# in a contact with 1 mole of liquid He3 at 0 bar:
block  Cu ${T0} type=paramagnet material=copper mass=100g
block  He ${T0} type=liquid_he3 volume=36cm^3 P=0bar

# They are linked by Kapitza thernal resistance which
# is proportional to contact surface area. Let's use
# area 1m^2 (it's not that big for sintered-powder heat exchangers)
link   k1 Cu He type=kap_res_he3 area=1m^2

# We will cool the system by mixing chamber (MC) of a simple dilution
# fridge with cooling power proportional to T^2, circulation 50 umol/s
# and external heat leak 0.2uW.

block  MC ${T0} type=simple C=0.2J/K # use some reasonable heat capacity
block  bath 300K type=bath   # this is needed as a sink for heat leaks
link MC_cooling MC bath type=dilution_cooling ndot=50umol/s
link MC_heatleak bath MC type=const Qdot=0.2uW

# The copper block is connected to the mixing chamber with a heat switch
# with electrical resistance 10 uOhm:
link precool MC Cu type=metal_bar R=10uOhm

# Magnetic field is 8T
field 8T

# We will be calculating with some large steps (30min). To keep the accuracy
# intermediate steps will be always done by the program. We can show/hide them
# by setting print_substeps to 1/0
print_substeps 0

# Now define what do we want to print: time, magnetic field,
# temperatures of a few blocks and heat flow through a link
print t B T(MC) T(Cu) T(He) Q(precool)

# Where we print result (use '-' for stdout):
print_to_file doc_example.dat

# Precool for 72h.
run 72h 30m

####
# Now we precooled to about 9mK. Let's disconnect HS and
# demagnetize to 0.4T (20 times) during 10h:

link precool MC Cu type=metal_bar R=1MOhm
field_rate -0.76T/h
run 10h 30m

# And now just wait
field_rate 0T/h
run 18h 30m

```

![Result](https://github.com/slazav/cryoblocks/blob/main/examples/doc_example.png)

### Command file

The system is defined as a set of `blocks` and `links` with certain
properties. Program reads a command file with the following structure:

* Comments start with `#`.

* Empty lines are skipped.

* `block <name> <temperature> [parameters]` -- Define a block with a given name and starting temperature.

* `link  <name> <block1> <block2> [parameters]` -- Define a named link which connects block1 and block2.

* `delete (block|link) <name>` -- Delete block or link.

* `run  <time> [step=<time step>] [abs=<0|1>] [to_field=<value>]` --
Do calculations for a given time period, print results. If `step` step
parameter is missing then a single step is used (internally it will be
divided by substeps to keep accuracy, substep results can be printed by
setting `print_substeps 1`). If `abs` parameter is missing or 0 then
`<time>` is a calculation period length. otherwise it's final time at the
end of the calculation (including shift set by `time_shift` command). If
`to_field` parameter is used, then field rate is set to reach the
specified field at the end of the calculation. Thern field rate is reset
to a previous value.

* `run_to  <final time> [<time step>]` -- DO NOT USE. Equivalent
to `run <final time> abs=1`.

* `include <file name>` -- Read commands from a file.

* `exit` -- Stop processing file and exit.


* `define <name> <value>` -- Define a parameter which can be used later.
`${<name>}` will be substituted by `<value>`.

* `field <value>` -- Set magnetic field.

* `field_rate <value>` -- Set magnetic field rate.

* `time_shift <value>` -- When printing result add the value to time.

* `max_tempstep <value>` -- Use smaller calculation steps to keep
relative change of temperatures on each step less than given value.
(Default: 1e-2).

* `max_tempacc <value>` -- Use smaller calculation steps to keep
relative temperature accuracy (difference between single dt and double
dt/2 steps) less than given value. (Default: 1e-6).


* `print <par1> ...` -- Set list of parameters for output during `run` command.
  Use `T(<name>)` for a block temperature, `Q(name)` for heat flow
  through a link, `t` for time, `B` for magnetic field. Example: `print t
  H T(block1) T(block2) Q(link1)`.

* `print_to_file <name>` -- Set filename for printing data. If name is `-` then
use `stdout` (default).

* `print_substeps 0|1` -- Print values for requested time steps or
for actual calculation steps (they could be smaller to achieve reasonable
accuracy).


* `read_data <file name> <par1> ...` -- Read time grid and additional parameters and
do calculation. Parameters determin data columns: `T(<name>)` for a block
temperature, `B` for magnetic field, `-` for skipping a column (numeric
or non-numeric).  If `<file name>` is `-` then data is read from the current
command file. Reading stops on `EOF` line or end of file.

* `read_factors <time factor> <par1 factor> ...` -- Set factors to
multiply values in `read_data` command. First value is factor for the time column.
Columns defined as `-` should be skipped here.
Example: to read time in hours and temperature of block `MC` in mK from 1st and 3rd
columns of a file use commands:
```
read_factors 3600 1e-3
read_data file.dat - T(MC)
```

#### Physics

Consider a block with temperature `T` with external heating power `dQ`.
Entropy of the block `S` depends on temperature and maybe other
parameters (magnetic field `B`).

Then `dQ = T dS = C dT + D dB`, where `C =  T \partial S / \partial T` is
heat capacity, and `D = T \partial S / \partial B` - cooling power of
demagnetization process. Using this
formula each block object can calculate change of temperature `dT` as a
function of `dQ`, `T`, `B`, and `dB`.

A thermal link object can calculate heat transfer depending on
temperatures of two blocks connected by the link.

Calculation is done using adaptive time steps. On each step initial
tempetures of all blocks are known. For each link heat flows are
calculated, Total power applied to each block is found. Using these
values new temperatures of each block are found. Blocks with zero heat
capacity are calculated differently: before and after each step a
zero of total heat flow on every such block is found as a function of
temperatures of these blocks.

#### Parameter reading

All parameters usually have some dimensions, it should be specified. For example,
time in `run` command can be written as `1s`, `0.1m`, `1e-2ms`, etc.

#### Blocks

Following types of blocks are supported:

* `block <name> <temperature> type=bath` -- A heat bath with intinite heat capacity and
constant temperature.

* `block <name> <temperature> type=zero-c` -- A block with zero heat capacity.

* `block <name> <temperature> type=simple C=<value> power=<p> factor=<f>` --
One can use `C` parameter to specify constant heat capacity in units
which can be converted to `J/K`. Alternatively one can give `power` and
`factor` dimensionless values, to have `C=<factor>*T^<power> [J/K]`. Zero
heat capacity is allowed.

* `block <name> <temperature> type=paramagnet [parameters]` -- paramagnetic material
(e.g. copper nuclei). Parameters:
  * `material=<v>` -- Set parameters `Bint`, `gyro`, `spin` for a certain material: `copper` (nuclei), `he3` (solid).
  * `Bint=<v>` -- Internal field (overrides value set by `material` parameter), default 0.
  * `gyro=<v>` -- Gyromagnetic ratio (overrides value set by `material` parameter).
  * `spin=<v>` -- Spin (overrides value set by `material` parameter).
  * `moles=<v>` -- Number of moles.
  * `mass=<v>` -- If material is set then mass can be used instead of moles.

* `block <name> <temperature> type=liquid_he3 [parameters]` -- liquid He3 (uses he3lib).
Parameters:
  * `P=<v>` -- Pressure
  * `moles=<v>` -- Number of moles
  * `mass=<v>` -- Can be used instead of moles.
  * `volume=<v>` -- Can be used instead of moles (pressure-dependent molar volume is used).


#### Links

Following types of links are supported:

* `link <name> <block1> <block2> type=const Qdot=<value>` -- A constant
heat transfer. Not a physical process, but good for tests.

* `link <name> <block1> <block2> type=simple K=<v> power=<p> factor=<f>` --
One can use `K` parameter to specify constant heat conductivity in units
which can be converted to `W/K`. Alternatively one can give `power` and
`factor` dimensionless values, to have `K=<factor>*T^<power> [W/K]`.
Heat transfer is calculated as `Qdot = <factor>*(T1^P - T2^P)/P`, where `P=<power>+1`.

* `link <name> <block1> <block2> type=metal_bar R=<v>` -- A bar made of metal.
Total resistance R, Wiedemann-Franz low is used to calculate heat conductivity.

* `link <name> <block1> <block2> type=field_heat_leak B=<v> B2=<v>` -- A field-dependent
heat leak, `Qdot = <B>*B + <B2>*B^2`. By default `B=0W/T, B2=0W/T^2`.

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
Electron-phonon coupling (`Qdot = C*(T1^5-T2^5)*nmoles`, see Pobell book f.10.9).
Parameters:

  * `material=<v>` -- predefined values for. Variants: `copper`
  * `C=<v>` -- set parameter C (overrides value set by `material` parameter), default 0
  * `moles=<v>` -- number of moles
  * `mass=<v>` -- if material is set then mass can be used instead of moles

* `link <name> <block1> <block2> type=kap_res_he3 area=<S> power=<N> C=<v>` --
Kapitza resistance between He3 and solids. `R = 900/T/S [W]`, or `R=41.0/T^2/S`,
or `R = 0.1/T^3/S` for temperature `T` and area `S`. Power should be 1, 2, or 3.
Parameter `C` is used to override numerical factor.

* `link <name> <block1> <block2> type=dilution_cooling ndot=<V>`
Cooling power of a dilution refrigerator.  block1 should be mixing
chamber, block2 - thermal bath at any temperature.

* `link <name> <block1> <block2> type=dilution_circ ndot=<V> phase=<C|D>`
Heat transfer by circulation in a dilution refrigerator. Phase parameter
is "C" or "D". T2 is not used in the calculation.
