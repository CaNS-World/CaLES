#!/usr/bin/python
# values should be consistent with dns.in
h     = 6.283185307179586
ub    = 1.
visci = 1600.
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
casename = '1600'.zfill(5)
