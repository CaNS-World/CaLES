#!/usr/bin/python
# values should be consistent with dns.in
h     = 1.
ub    = 1.
visci = 10000.
#
uconv = 0. # if we solve on a convective reference frame; else = 0.
#
# parameters for averaging
#
tbeg   = 200.
tend   = 800.
fldstp = 100
#
# case name (e.g., the Retau)
#
casename = '550'.zfill(5)
