########
Examples
########

aospy comes with some `sample data files
<https://github.com/spencerahill/aospy/tree/develop/aospy/test/data/netcdf>`_,
which can be used to illustrate some of its basic features.  These
files contain monthly mean time series output of two variables (the
large-scale and convective components of the total precipitation rate)
from an idealized aquaplanet climate model.  A simple computation one
could seek to do from this model output would be to compute some
statistics of the total precipitation rate (large-scale plus
convective).

Here's a quick summary of the included data:

.. ipython:: python
   
   import xarray as xr
   xr.open_mfdataset('../aospy/test/data/netcdf/000[4-6]0101.precip_monthly.nc',
                     decode_times=False)

In this particular model, the large-scale component of the precipitation rate
is called "condensation_rain" and the convective component is called
"convection_rain."

Defining the simulation metadata
================================

The first step is to create an :py:class:`aospy.Run` object that
stores metadata about this simulation.  This includes giving it a
name, a description, and specifying where files are located through a
DataLoader.

.. ipython:: python

    from datetime import datetime
    
    from aospy import Run
    from aospy.data_loader import DictDataLoader
    file_map = {'monthly': '../aospy/test/data/netcdf/000[4-6]0101.precip_monthly.nc'}
    example_run = Run(
        name='example_run',
        description=(
            'Control simulation of the idealized moist model'
        ),
        data_loader=DictDataLoader(file_map)
    )
    
We then need to associate this ``Run`` with an :py:class:`aospy.Model` object:

.. ipython:: python

    from aospy import Model
    example_model = Model(
        name='example_model',
        grid_file_paths=(
            '../aospy/test/data/netcdf/00040101.precip_monthly.nc',
            '../aospy/test/data/netcdf/im.landmask.nc'
        ),
        runs=[example_run]
    )

Finally, we need to associate the ``Model`` object with an
:py:class:`aospy.Proj` object.  Here we can specify the location that
aospy will save its output files.

.. ipython:: python

    from aospy import Proj
    example_proj = Proj(
        'example_proj',
        direc_out='example-output',
        tar_direc_out='example-tar-output',
        models=(example_model,)
    )

Now the metadata associated with this simulation is fully defined.  We
can move on to computing the total precipitation.

Computing the annual mean total precipitation rate
==================================================

We can start by defining a simple
python function that computes the total precipitation from condensation and
convection rain arguments:

.. ipython:: python

    def total_precipitation(condensation_rain, convection_rain):
        return condensation_rain + convection_rain

To hook this function into the aospy framework, we need to connect it
to an :py:class:`aospy.Var` object, as well as define the ``Var``
objects it depends on (variables that are natively stored in model
output files).

.. ipython:: python

   from aospy import Var
   condensation_rain = Var(
       name='condensation_rain',
       alt_names=('prec_ls',),
       def_time=True,
       description=('condensation rain'),
   )

   convection_rain = Var(
       name='convection_rain',
       alt_names=('prec_conv',),
       def_time=True,
       description=('convection rain'),
   )

   precip = Var(
       name='total_precipitation',
       def_time=True,
       description=('total precipitation rate'),
       func=total_precipitation,
       variables=(condensation_rain, convection_rain)
   )

Here the func attribute of the precip ``Var`` object is the function
we defined, and the variables attribute is a tuple containing the
``Var`` objects the function depends on, in the order of the
function's call signature.

If we'd like to compute the time-mean total precipitation rate from
year four to year six using aospy, we can create an
:py:class:`aospy.Calc` object.  This is currently done through passing
an :py:class:`aospy.CalcInterface` object to a ``Calc`` object; once
created, the computation can be submitted by simply calling the
compute function of ``Calc``.

.. ipython:: python

    from aospy import CalcInterface, Calc
    calc_int = CalcInterface(
        proj=example_proj,
        model=example_model,
        run=example_run,
        var=precip,
        date_range=(datetime(4, 1, 1), datetime(6, 12, 31)),
        intvl_in='monthly',
        dtype_in_time='ts',
        intvl_out='ann',
        dtype_out_time='av'
    )
    Calc(calc_int).compute()

The result is stored in a netcdf file, whose path and filename
contains metadata about where it came from:

.. ipython:: python
             
    calc_int.path_out['av']

Using xarray we can open and plot the results of the calculation:

.. ipython:: python

    @savefig plot_ann_total_precipitation.png width=80%
    xr.open_dataset(calc_int.path_out['av']).total_precipitation.plot()

Computing the global annual mean total precipitation rate
=========================================================

Not only does aospy enable reductions along the time dimension, it
also enables area weighted regional averages.  As a simple
introduction, we'll show how to compute the global mean total
precipitation rate from this ``Run``.  To do so, we'll make use of the
infrastructure defined above, and also define an
:py:class:`aospy.Region` object:

.. ipython:: python

    from aospy import Region
    globe = Region(
        name='globe',
        description='Entire globe',
        lat_bounds=(-90, 90),
        lon_bounds=(0, 360),
        do_land_mask=False
    )

To compute the global annual mean total precipitation rate, we can now create
another ``Calc`` object:

.. ipython:: python

   calc_int = CalcInterface(
       proj=example_proj,
       model=example_model,
       run=example_run,
       var=precip,
       date_range=(datetime(4, 1, 1), datetime(6, 12, 31)),
       intvl_in='monthly',
       dtype_in_time='ts',
       intvl_out='ann',
       dtype_out_time='reg.av',
       region={'globe': globe}
   )
   Calc(calc_int).compute()

This produces a new file, located in:

.. ipython:: python

   calc_int.path_out['reg.av']

We find that the global annual mean total precipitation rate for this
run (converting to units of mm per day) is:

.. ipython:: python

    xr.open_dataset(calc_int.path_out['reg.av']).globe * 86400.


.. ipython:: python
    :suppress:
       
    from shutil import rmtree
    rmtree('example-output')
    rmtree('example-tar-output')
