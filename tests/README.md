Directory for unittests for jannasutils

## Running the unit tests

from the jannas-utils git repository run:
$ python3 -m unittest

## Conventions for creating unit tests
Create a separate testSOMETHING.py file for each module SOMETHING.py Each class represents a seperate unit test; i.e. testing a specific functionality in the module. Within the unit test, different methods can exist (for example to test for different inputs).
