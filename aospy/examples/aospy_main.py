#! /usr/bin/env python
"""Main script for executing calculations using aospy.

Before using this script
------------------------

It is best to copy this template into a separate directory before populating it
with the objects from your own object library.  You can also always get a fresh
copy from https://github.com/spencerahill/aospy/examples/aospy_main.py

How to use this script
----------------------

On the example library and data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script comes pre-populated with objects taken from the example aospy
object library that is included in this directory in the `example_obj_lib.py`
module. So you can try it out on the sample data without modifying anything at
all.

This simple example library includes only one Proj, one Model, and one Run, but
it also includes multiple Var and Region objects over which you can automate
computations.  The date range, sub-annual averaging period, and output data
types can all also be modified.

On your own library and data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a user, there are only two places you need to modify the code:

1. (Only done once) Replace `example_obj_lib` with the name of your object
   library.  Consult the documentation for instructions on how to make your
   object library findable by Python if this statement generates errors.

2. (Done as often as desired) Replace the dummy project, model, run, var, and
   region objects with your objects that you want to perform
   calculations with.

   Also alter the other parameters -- date_range, etc. -- to your liking.

Running the script
------------------

Once the parameters are all set as desired, execute the script from the command
line ::

  ./aospy_main.py  # after `cd` to the directory where you've made your copy

"""
# C.f. instructions above, replace `example_obj_lib` with your own object
# library if you wish to use this on your own data.
from aospy import submit_mult_calcs
import example_obj_lib as lib


# This dictionary contains all of the specifications of calculations that you
# wish to permute over.
calc_suite_specs = dict(
    # Consult `Calc` API reference for further explanation of each option and
    # accepted values.

    # The desired library of aospy objects.
    library=lib,
    # List of the Proj objects, or 'default', or 'all'.
    projects=[lib.example_proj],
    # List of the Model objects, or 'default', or 'all'.
    models=[lib.example_model],
    # List of the Run objects, or 'default', or 'all'.
    runs=[lib.example_run],
    # List of the Var objects, or 'default', or 'all'.
    variables=[lib.precip_largescale, lib.precip_convective, lib.precip_total,
               lib.precip_conv_frac],
    # List of the Region objects, or 'default', or 'all'.
    regions='all',

    # Start and end dates (inclusive).  Either 'default' or a list comprising
    # tuples of the form (start_date, end_date), where start_date and end_date
    # are datetime.datetime objects.  Be sure to add `import datetime` above if
    # using `datetime.datetime` objects.
    date_ranges='default',

    # Sub-annual time-interval to average over.  List of 'ann', seasonal
    # string (e.g. 'djf'), or month integer (1 for Jan, 2 for Feb, etc).
    output_time_intervals=['ann'],
    # List of strings indicating the desired spatiotemporal reductions.
    output_time_regional_reductions=['av', 'std', 'reg.av', 'reg.std',
                                     'reg.ts'],
    # List of desired vertical reductions to perform.
    output_vertical_reductions=[None],

    # List of time spacing of input data.
    input_time_intervals=['monthly'],
    # List of time type of input data.
    input_time_datatypes=['ts'],
    # List the time offset dictionaries (if desired) to apply to the input
    # data (e.g. [{'days': -15}]).
    input_time_offsets=[None],
    # List of vertical data type of input data.
    input_vertical_datatypes=[False],

)


# This dictionary contains options regarding how the calculations are displayed
# to you, submitted for execution, and saved upon execution.
calc_exec_options = dict(
    # List calculations to be performed and prompt for your verification before
    # submitting them for execution.
    prompt_verify=True,

    # Submit all calculations in parallel.  If parallelize is True and client
    # is None, a LocalCluster will be started; the client argument can be used
    # to specify an external dask.distributed Client for use in parallelizing
    # computations
    parallelize=False,
    client=None,

    # Save results of calculations to .tar files, one for each Run object.
    # These tar files are placed using the same directory structure as the
    # standard output relative to their root directory, which is specified via
    # the `tar_direc_out` argument of each Proj object's instantiation.
    write_to_tar=True,
)


# Don't modify this statement.
if __name__ == '__main__':
    calcs = submit_mult_calcs(calc_suite_specs, calc_exec_options)
