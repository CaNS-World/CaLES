#!/usr/bin/python
# values should be consistent with dns.in
h     = 1.
ub    = 1.
visci = 180.
#
uconv = 0. # if we solve on a convective reference frame; else = 0.
#
# parameters for averaging
#
tbeg   = 60000.
tend   = 100000.
fldstp = 5000
#
# case name (e.g., the Retau)
#
casename = '180'.zfill(5)
