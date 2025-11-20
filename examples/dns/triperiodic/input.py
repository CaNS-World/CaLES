#!/usr/bin/python
# values should be consistent with dns.in
h     = 1.
ub    = 1.
visci = 1000.
#
uconv = 0. # if we solve on a convective reference frame; else = 0.
#
# parameters for averaging
#
tbeg   = 6000.
tend   = 10000.
fldstp = 500
#
# case name (e.g., the Retau)
#
casename = '1000'.zfill(5)
