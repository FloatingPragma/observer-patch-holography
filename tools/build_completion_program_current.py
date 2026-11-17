#!/usr/bin/env python3
"""Public entry point for the official common-source acceptance programme."""
try:
    from tools.build_completion_program import *
except ModuleNotFoundError:
    from build_completion_program import *

if __name__=='__main__':main()
