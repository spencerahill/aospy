Climate data analysis, database, and visualization package for Python

*Installation*

User-defined aospy objects (Proj, Calc, Var, Region, Units) must be defined in `~/aospy`, with projects in that directory as their own modules (e.g. `aero_3agcm.py` for a project named `aero_3agcm`), while the others must be subdirectories and stored in `__init__.py` files (e.g. `calcs/__init__.py`, `regions/__init__.py`, etc.).  This is an imperfect and temporary solution.

With this setup, aospy will import all of the users objects and be fully functional.

------

This version of aospy gives the user the option to use the xray package internally to do calculations. Currently this is supported on a run-by-run basis; really this should be an option in the main script. That way one could easily toggle it based on the calculations one wanted to do (rather than have to change each run object individually). 
