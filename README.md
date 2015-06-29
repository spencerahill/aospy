Climate data analysis, database, and visualization package for Python

*Installation*

User-defined aospy objects (Proj, Calc, Var, Region, Units) must be defined in `~/aospy`, with projects in that directory as their own modules (e.g. `aero_3agcm.py` for a project named `aero_3agcm`), while the others must be subdirectories and stored in `__init__.py` files (e.g. `calcs/__init__.py`, `regions/__init__.py`, etc.).  This is an imperfect and temporary solution.

With this setup, aospy will import all of the users objects and be fully functional.
