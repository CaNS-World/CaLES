#!/usr/bin/python
# values should be consistent with dns.in
h     = 1.
ub    = 1.
visci = 100.
#
uconv = 0. # if we solve on a convective reference frame; else = 0.
#
# parameters for averaging
#
tbeg   = 600.
tend   = 1000.
fldstp = 50
#
# case name (e.g., the Retau)
#
casename = '100'.zfill(5)
