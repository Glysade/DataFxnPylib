
Data Function Background
========================

This document contains information on building Glysade Python data functions for Spotfire.

To quickly understand how to construct a data function script that can be included in the data function definition
YAML, see the examples in :ref:`getting-started`.

A core API and utility functions for handling chemical and biological data is described in :ref:`core-api`.  This core
API should be sufficient to build most data functions.

The library contains a number of builtin data functions described in :ref:`df-modules`.  These data functions
are executed by dynamically importing modules and executing a class method, rather than running a script.

The remainder of the documentation concerns the *ruse* packages.  These modules are required by the builtin data
functions.  If possible, avoid using these modules as they may change in the future.
