import numpy as np
#
# profile folding function
#
def fold(var,cf='C',isym=1):
   if(abs(isym) != 1): exit()
   n = np.size(var[:])
   s = 0
   if(cf == 'F'): s = 1
   var[0:n//2-s:+1] = 0.5*(var[0:n//2-s:+1]+isym*var[n-1-s:n//2-1:-1])
   var[n-1-s:n//2-1:-1] = isym*var[0:n//2-s:+1]
   return var
def interp(var,bc='D',cf='F'):
   n = np.size(var[:])
   #
   var_aux = np.zeros(n+2)
   var_aux[1:n+1] = var[:]
   if(  bc + cf == 'DF'):
       var_aux[0] = 0.
       var_aux[n] = 0.
       var_aux[n+1] = var_aux[n-1] # not needed
   elif(bc + cf == 'NF'):
       var_aux[0] = var_aux[1  ]
       var_aux[n] = var_aux[n-1]
   elif(bc + cf == 'DC'):
       var_aux[0  ] = -var_aux[1]
       var_aux[n+1] = -var_aux[n]
   elif(bc + cf == 'NC'):
       var_aux[0  ] = var_aux[1]
       var_aux[n+1] = var_aux[n]
   #
   res = np.zeros(n)
   if(  cf == 'F'):
       res[:] = (var_aux[2:n+2] + var_aux[1:n+1])/2
   elif(cf == 'C'):
       res[:] = (var_aux[1:n+1] + var_aux[0:n])/2
   return res
def ddz(var,dzc,dzf,bc='N',cf='F',lz=1.,order=1):
   n = np.size(var[:])
   dzf_aux  = np.zeros(n+2)
   dzf_aux[1:n+1] = dzf[:]
   dzf_aux[0] = dzf_aux[1]
   dzf_aux[n+1] = dzf_aux[n]
   dzc_aux  = np.zeros(n+2)
   dzc_aux[1:n+1] = dzc[:]
   dzc_aux[0] = dzf_aux[0]
   dzc_aux[n+1] = dzf_aux[n+1]
   #
   var_aux = np.zeros(n+2)
   var_aux[1:n+1] = var[:]
   if(  bc + cf == 'DF'):
       var_aux[0] = 0.
       var_aux[n] = 0.
       var_aux[n+1] = var_aux[n-1] # not needed
   elif(bc + cf == 'NF'):
       var_aux[0] = var_aux[1  ]
       var_aux[n] = var_aux[n-1]
   elif(bc + cf == 'DC'):
       var_aux[0  ] = -var_aux[1]
       var_aux[n+1] = -var_aux[n]
   elif(bc + cf == 'NC'):
       var_aux[0  ] = var_aux[1]
       var_aux[n+1] = var_aux[n]
   #
   res = np.zeros(n)
   if(order == 1):
       if(  cf == 'F'):
           res[:] = (var_aux[1:n+1] - var_aux[0:n])/dzf_aux[1:n+1]
       elif(cf == 'C'):
           res[0:n-1] = (var_aux[2:n+1] - var_aux[1:n])/dzc_aux[1:n]
           res[n-1] = res[n-2] # not needed
   return res
def estim_std_var(signal,m_max=100):
    #
    # see Appendix of Hoyas & Jimenez, PoF 20 101511 (2008)
    #
    n = np.size(signal)
    m_max = min(m_max,n)
    m_min = min(2    ,n) # 2 as the minimum signal subsampling
    mean = np.average(signal[:])
    var  = np.average(signal[:]**2)-mean**2
    std = (var/n)**.5
    #print('stdev/mean, naive calculation: ',std,std/abs(mean))
    #arr = np.zeros([m_max,3])
    for m in range(m_max,m_min-1,-1):
        nm = n//m
        mean = 0.
        var  = 0.
        for j in range(m):
            lbound = j*nm
            ubound = min(n,(j+1)*nm)
            mean += np.average(signal[lbound:ubound])
            var  += np.average(signal[lbound:ubound]**2)
            ##
            ## alternatively, use midpoint value
            ##
            #mbound = (lbound+ubound)//2
            #mean += signal[mbound]
            #var  += signal[mbound]**2
        mean = mean/m
        std = ((var/m-mean**2)/m)**.5
        #print('Stdev/mean, m = ', m,':',std,std/abs(mean))
        #arr[m-1,0:3] = np.array([m,mean,std])
#
# compute mean pressure gradient
#
use_new_arrays = False
shuffle_data   = False
import sys
args = sys.argv[1:]
casedir = args[3] if len(args) > 3 else ''
datadir = casedir+''
print('Data directory: ' + datadir+'\n')
resultsdir = casedir+'results/'
import os
os.makedirs(resultsdir,exist_ok=True)
#
# import parameters
#
sys.path.append(casedir) # put input.py in the case directory
from input import *
tbeg   = float(args[0]) if len(args) > 0 else tbeg
tend   = float(args[1]) if len(args) > 1 else tend
fldstp = int(  args[2]) if len(args) > 2 else fldstp
fname_dpdx = datadir+'forcing.out'
data = np.loadtxt(fname_dpdx)
tarr, ind = np.unique(data[:,0], return_index=True) # this requires the same number of tasks used for restart
if(shuffle_data):
    print('Note: dataset will be randomly shuffled before fixing the time window\n')
    nlen = np.size(data[ind,1][np.where((tarr>tbeg) & (tarr<tend))])
    dpdx_arr = data[ind,1][np.where(tarr>tbeg)]
    np.random.shuffle(dpdx_arr[:])
    ub_arr = data[ind,4][np.where(tarr>tbeg)]
    np.random.shuffle(ub_arr[:])
else:
    dpdx_arr = data[ind,1][np.where((tarr>tbeg) & (tarr<tend))]
    nlen = np.size(dpdx_arr[:])
    ub_arr = data[ind,4][np.where((tarr>tbeg) & (tarr<tend))]
dpdx = -np.average(dpdx_arr[0:nlen])
ub   =  np.average(ub_arr  [0:nlen])
print("ub = ", ub)
#
visc = visci**(-1)
reb  = 2*h*ub/visc
#
utau = np.sqrt(dpdx*h)
retau = utau*h/visc
dnu   = visc/utau
cf    = utau**2/(ub**2/2.)
np.savetxt(resultsdir+'stats.txt',np.c_[retau,utau,dnu])
#
print("Pressure gradient = ", dpdx)
print("u_tau/u_bulk = ", np.sqrt(dpdx*h)/ub)
print("dnu/h        = ", dnu/h)
print("Friction Reynolds number = ", np.sqrt(dpdx*h)*h/visc)
print("Bulk Reynolds number = ", reb)
#
# compute single point statistics
#
saves_log = np.loadtxt(datadir+'time.out')
step_log, ind = np.unique(saves_log[:,0], return_index=True)
time_log = saves_log[ind,2]
def find_closest_save(time_log,step_log,time,stp):
    index = np.searchsorted(time_log, time, side="right")
    if(index < np.size(step_log[:])):
        isave = int(step_log[index])-stp//2
        isave = isave - isave%stp
    else:
        isave = -1
    return isave
fldbeg = find_closest_save(time_log,step_log,tbeg,fldstp)
fldend = find_closest_save(time_log,step_log,tend,fldstp)
import glob
import re
import os
files = np.sort(glob.glob(datadir+'velstats_fld_???????.out'))
lastfile = files[-1]
fldmax = list(map(int,re.findall(r'([0-9]{7,})', lastfile)))[-1]
if(fldend == -1):
    fldend = fldmax
    index = np.searchsorted(step_log[:], fldmax, side="left")
    print("\nN.B.: Maximum time exceeds number of saves, overwriting end time...")
    tend = time_log[index]
firstfile= files[0]
fldmin = list(map(int,re.findall(r'([0-9]{7,})',firstfile)))[-1]
if(fldmin > fldbeg):
    fldbeg = fldmin
    index = np.searchsorted(step_log[:], fldmin, side="left")
    print("\nN.B.: Minimum larger than the first found save, overwriting  time...")
    tbeg = time_log[index]
if(fldbeg > fldend):
    print("\nError in the selected time interval. Aborting...")
    exit()
t_dns = 1.
tb    = (2.*h)/ub
te    = h/utau
print("Max, min, and interval averaging time in half bulk times (h/ub) = ", tbeg         , tend         , (tend-tbeg)         ,"\n")
print("Max, min, and interval averaging time in turnover times         = ", tbeg*t_dns/te, tend*t_dns/te, (tend-tbeg)*t_dns/te,"\n")
#
fname_beg = datadir+'velstats_fld_'
fname_ext = '.out'

nflds = (fldend-fldbeg)//fldstp + 1
flds = np.array(range(fldbeg,fldend+fldstp,fldstp))
if(shuffle_data):
    np.random.shuffle(flds)
flds = flds[0:nflds]
norm = 1/nflds
grid = np.loadtxt(datadir+'geometry.out')
nz   = int(grid[0,2])
#
# single point statistics
#
izc    = 0
izf    = 1
iu1    = 2
iv1    = 3
iw1    = 4
iu2    = 5
iv2    = 6
iw2    = 7
iuw    = 8
iu3    = 9
iv3    = 10
iw3    = 11
iu4    = 12
iv4    = 13
iw4    = 14
ip1    = 15
ip2    = 16
iomx1  = 17
iomy1  = 18
iomz1  = 19
iomx2  = 20
iomy2  = 21
iomz2  = 22
iu2m   = 23
iv2m   = 24
iw2m   = 25
iuwm   = 26
ivisct = 27
iuwv   = 28

idzc  = -2
idzf  = -1

u1    = np.zeros(nz)
v1    = np.zeros(nz)
w1    = np.zeros(nz)
u2    = np.zeros(nz)
v2    = np.zeros(nz)
w2    = np.zeros(nz)
uw    = np.zeros(nz)
u3    = np.zeros(nz)
v3    = np.zeros(nz)
w3    = np.zeros(nz)
u4    = np.zeros(nz)
v4    = np.zeros(nz)
w4    = np.zeros(nz)
p1    = np.zeros(nz)
p2    = np.zeros(nz)
omx1  = np.zeros(nz)
omy1  = np.zeros(nz)
omz1  = np.zeros(nz)
omx2  = np.zeros(nz)
omy2  = np.zeros(nz)
omz2  = np.zeros(nz)
u2m   = np.zeros(nz)
v2m   = np.zeros(nz)
w2m   = np.zeros(nz)
uwm   = np.zeros(nz)
visct = np.zeros(nz)
uwv   = np.zeros(nz)
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_ext
    print(fname)
    data = np.loadtxt(fname)
    #
    zc        = data[:,izc]
    zf        = data[:,izf]
    dzc       = data[:,idzc]
    dzf       = data[:,idzf]
    u1[:]    += data[:,iu1]
    v1[:]    += data[:,iv1]
    w1[:]    += data[:,iw1]
    u2[:]    += data[:,iu2]
    v2[:]    += data[:,iv2]
    w2[:]    += data[:,iw2]
    uw[:]    += data[:,iuw]
    u3[:]    += data[:,iu3]
    v3[:]    += data[:,iv3]
    w3[:]    += data[:,iw3]
    u4[:]    += data[:,iu4]
    v4[:]    += data[:,iv4]
    w4[:]    += data[:,iw4]
    p1[:]    += data[:,ip1]
    p2[:]    += data[:,ip2]
    omx1[:]  += data[:,iomx1]
    omy1[:]  += data[:,iomy1]
    omz1[:]  += data[:,iomz1]
    omx2[:]  += data[:,iomx2]
    omy2[:]  += data[:,iomy2]
    omz2[:]  += data[:,iomz2]
    u2m[:]   += data[:,iu2m]
    v2m[:]   += data[:,iv2m]
    w2m[:]   += data[:,iw2m]
    uwm[:]   += data[:,iuwm]
    visct[:] += data[:,ivisct]
    uwv[:]   += data[:,iuwv]
u1[:]    = fold(u1[:]  ,'C',+1)*norm
v1[:]    = fold(v1[:]  ,'C',+1)*norm
w1[:]    = fold(w1[:]  ,'F',-1)*norm
u2[:]    = fold(u2[:]  ,'C',+1)*norm
v2[:]    = fold(v2[:]  ,'C',+1)*norm
w2[:]    = fold(w2[:]  ,'F',+1)*norm
uw[:]    = fold(uw[:]  ,'F',-1)*norm
u3[:]    = fold(u3[:]  ,'C',+1)*norm
v3[:]    = fold(v3[:]  ,'C',+1)*norm
w3[:]    = fold(w3[:]  ,'F',-1)*norm
u4[:]    = fold(u4[:]  ,'C',+1)*norm
v4[:]    = fold(v4[:]  ,'C',+1)*norm
w4[:]    = fold(w4[:]  ,'F',+1)*norm
p1[:]    = fold(p1[:]  ,'C',+1)*norm
p2[:]    = fold(p2[:]  ,'C',+1)*norm
omx1[:]  = fold(omx1[:],'F',+1)*norm
omy1[:]  = fold(omy1[:],'F',-1)*norm
omz1[:]  = fold(omz1[:],'C',+1)*norm
omx2[:]  = fold(omx2[:],'F',+1)*norm
omy2[:]  = fold(omy2[:],'F',+1)*norm
omz2[:]  = fold(omz2[:],'C',+1)*norm
u2m[:]   = fold(u2m[:] ,'C',+1)*norm
v2m[:]   = fold(v2m[:] ,'C',+1)*norm
w2m[:]   = fold(w2m[:] ,'F',+1)*norm
uwm[:]   = fold(uwm[:] ,'F',-1)*norm
visct[:] = fold(visct[:],'C',+1)*norm
uwv[:]   = fold(uwv[:] ,'F',-1)*norm
#
u2[:] -= u1[:]**2
v2[:] -= v1[:]**2
w2[:] -= w1[:]**2
p2[:] -= p1[:]**2
omx2[:] -= omx1[:]**2
omy2[:] -= omy1[:]**2
omz2[:] -= omz1[:]**2
#uw[:] -= interp(u1[:],'D',C')*w1[:] # unecessary -> w1[:] is zero (to machine precision, actually)
u3[:] -= 3*u2[:]*u1[:] + u1[:]**3
v3[:] -= 3*v2[:]*v1[:] + v1[:]**3
w3[:] -= 3*w2[:]*w1[:] + w1[:]**3
u4[:] -= 6*u2[:]*u1[:]**2 + 4*u3[:]*u1[:] + u1[:]**4
v4[:] -= 6*v2[:]*v1[:]**2 + 4*v3[:]*v1[:] + v1[:]**4
w4[:] -= 6*w2[:]*w1[:]**2 + 4*w3[:]*w1[:] + w1[:]**4
uwv[:] = -visc*uwv[:]
#
fname = resultsdir + 'stats-single-point-chan-' + casename + fname_ext
np.savetxt(fname,np.c_[zc,zf,u1,v1,w1,u2,v2,w2,uw,u3,v3,w3,u4,v4,w4,p1,p2,omx1,omy1,omz1,omx2,omy2,omz2,u2m,v2m,w2m,uwm,visct,uwv])
utau_s = ((u1[0]+uconv)/zc[0]*visc)**.5
uc     = u1[nz//2-1]+ub # convective reference frame
uu_max = np.max(u2[:])
uw_max = np.max(uw[:])
print("Friction Reynolds number obtained from the velocity gradient at the wall: ", utau_s*h/visc)
print("Centreline velocity:                  ", uc)
print("Maximum streamwise velocity variance: ", uu_max)
print("Maximum Reynolds shear stresses:      ", uw_max)
#
# MKE and shear stress budgets
#
fname_beg = datadir+'velstats_fld_'
fname_mid = '_reystr_budget'
fname_ext = '.out'
#
iz_c      = 0
iz_f      = 1
iu1_c     = 2
iu1_f     = 3
idu1dz1_f = 4
idu2dz1_f = 5
iuw_f     = 6
iuw_c     = 7
idu1dz1_c = 8
#
z_c      = np.zeros(nz)
z_f      = np.zeros(nz)
u1_c     = np.zeros(nz)
u1_f     = np.zeros(nz)
du1dz1_f = np.zeros(nz)
du2dz1_f = np.zeros(nz)
uw_f     = np.zeros(nz)
uw_c     = np.zeros(nz)
du1dz1_c = np.zeros(nz)
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_mid+fname_ext
    data = np.loadtxt(fname)
    #
    z_c          = data[:,iz_c]
    z_f          = data[:,iz_f]
    u1_c[:]     += data[:,iu1_c]
    u1_f[:]     += data[:,iu1_f]
    du1dz1_f[:] += data[:,idu1dz1_f]
    du2dz1_f[:] += data[:,idu2dz1_f]
    uw_f[:]     += data[:,iuw_f]
    uw_c[:]     += data[:,iuw_c]
    du1dz1_c[:] += data[:,idu1dz1_c]
u1_c[:]     = u1_c[:]    *norm
u1_f[:]     = u1_f[:]    *norm
du1dz1_f[:] = du1dz1_f[:]*norm
du2dz1_f[:] = du2dz1_f[:]*norm
uw_f[:]     = uw_f[:]    *norm
uw_c[:]     = uw_c[:]    *norm
du1dz1_c[:] = du1dz1_c[:]*norm
#
# shear stresses budget
#
term1 = fold(-uw_c[:],'C',-1)
term2 = fold(du1dz1_c[:]*visc,'C',-1)
tot   = term1 + term2
fname = resultsdir + 'stats-single-point-chan-shear-stress-balance-' + casename + fname_ext
np.savetxt(fname,np.c_[zc,tot,term1,term2])
#
# MKE budget
#
prod = fold(dpdx*u1_c[:],'C',+1)
diss = fold(-visc*du1dz1_c[:]**2,'C',+1)
work = fold(uw_c[:]*du1dz1_c[:],'C',+1)
trans_visc = fold(visc*ddz(ddz(u1_c**2/2.,dzc,dzf,'D','C'),dzc,dzf,'N','F'),'C',+1)
trans_turb = fold(ddz(-u1_f[:]*uw_f[:],dzc,dzf,'D','F'),'C',+1)
trans = trans_visc[:]+trans_turb[:]
tot = prod[:]+diss[:]+work[:]+trans[:]
fname = resultsdir + 'stats-single-point-chan-mke-budget-' + casename + fname_ext
np.savetxt(fname,np.c_[zc[:],tot,prod,diss,work,trans_visc,trans_turb])
#
# uu budget
#
iuuw_f     = 9
ip_c       = 10
ipdudx_c   = 11
idiss_uu_c = 12
idiss_uu_x_c = 31
idiss_uu_y_c = 32
idiss_uu_z_f = 33
#
uuw_f     = np.zeros(nz)
p_c       = np.zeros(nz)
pdudx_c   = np.zeros(nz)
diss_uu_c = np.zeros(nz)
if(use_new_arrays):
    diss_uu_x_c = np.zeros(nz)
    diss_uu_y_c = np.zeros(nz)
    diss_uu_z_f = np.zeros(nz)
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_mid+fname_ext
    data = np.loadtxt(fname)
    #
    uuw_f[:]     += data[:,iuuw_f]
    p_c[:]       += data[:,ip_c]
    pdudx_c[:]   += data[:,ipdudx_c]
    diss_uu_c[:] += data[:,idiss_uu_c]
    if(use_new_arrays):
        diss_uu_x_c[:] += data[:,idiss_uu_x_c]
        diss_uu_y_c[:] += data[:,idiss_uu_y_c]
        diss_uu_z_f[:] += data[:,idiss_uu_z_f]
uuw_f[:]     = uuw_f[:]    *norm
p_c[:]       = p_c[:]      *norm
pdudx_c[:]   = pdudx_c[:]  *norm
diss_uu_c[:] = diss_uu_c[:]*norm
if(use_new_arrays):
    diss_uu_x_c[:] = diss_uu_x_c[:]*norm
    diss_uu_y_c[:] = diss_uu_y_c[:]*norm
    diss_uu_z_f[:] = diss_uu_z_f[:]*norm
    diss_uu_z_c = interp(diss_uu_z_f[:]-du1dz1_f[:]**2,'N','F')
#
prod  = -work[:]
diss  = fold(-visc*(diss_uu_c[:]-du1dz1_c[:]**2),'C',+1)
if(use_new_arrays):
    diss  = fold(-visc*(diss_uu_x_c[:]+diss_uu_y_c[:]+diss_uu_z_c[:]),'C',+1)
dist  = fold(pdudx_c[:],'C',+1)
trans_visc_uu = fold(visc*ddz(du2dz1_f[:]-ddz(u1_c**2,dzc,dzf,'D','C'),dzc,dzf,'N','F')/2,'C',+1)
trans_turb_uu = fold(-ddz(uuw_f[:]-2.*uw_f[:]*u1_f[:],dzc,dzf,'D')/2.,'C',+1)
prod_uu  = prod[:]
diss_uu  = diss[:]
dist_uu  = dist[:]
trans_uu = trans_visc_uu[:]+trans_turb_uu[:]
tot_uu = prod_uu[:]+diss_uu[:]+dist_uu[:]+trans_uu[:]
fname = resultsdir + 'stats-single-point-chan-uu-budget-' + casename + fname_ext
np.savetxt(fname,np.c_[zc[:],tot_uu,prod_uu,diss_uu,dist_uu,trans_visc_uu,trans_turb_uu])
#
# vv budget
#
idv2dz1_f  = 13
ivvw_f     = 14
ipdvdy_c   = 15
idiss_vv_c = 16
idiss_vv_x_c = 34
idiss_vv_y_c = 35
idiss_vv_z_f = 36
#
dv2dz1_f  = np.zeros(nz)
vvw_f     = np.zeros(nz)
pdvdy_c   = np.zeros(nz)
diss_vv_c = np.zeros(nz)
if(use_new_arrays):
    diss_vv_x_c = np.zeros(nz)
    diss_vv_y_c = np.zeros(nz)
    diss_vv_z_f = np.zeros(nz)
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_mid+fname_ext
    data = np.loadtxt(fname)
    #
    dv2dz1_f[:]  += data[:,idv2dz1_f ]
    vvw_f[:]     += data[:,ivvw_f    ]
    pdvdy_c[:]   += data[:,ipdvdy_c  ]
    diss_vv_c[:] += data[:,idiss_vv_c]
    if(use_new_arrays):
        diss_vv_x_c[:] += data[:,idiss_vv_x_c]
        diss_vv_y_c[:] += data[:,idiss_vv_y_c]
        diss_vv_z_f[:] += data[:,idiss_vv_z_f]
dv2dz1_f[:]  = dv2dz1_f[:] *norm
vvw_f[:]     = vvw_f[:]    *norm
pdvdy_c[:]   = pdvdy_c[:]  *norm
diss_vv_c[:] = diss_vv_c[:]*norm
if(use_new_arrays):
    diss_vv_x_c[:] = diss_vv_x_c[:]*norm
    diss_vv_y_c[:] = diss_vv_y_c[:]*norm
    diss_vv_z_f[:] = diss_vv_z_f[:]*norm
    diss_vv_z_c = interp(diss_vv_z_f[:],'N','F')
#
diss_vv  = fold(-visc*(diss_vv_c[:]),'C',+1)
if(use_new_arrays):
    diss_vv  = fold(-visc*(diss_vv_x_c[:]+diss_vv_y_c[:]+diss_vv_z_c[:]),'C',+1)
dist_vv  = fold(pdvdy_c[:],'C',+1)
trans_turb_vv = fold(-ddz(vvw_f[:],dzc,dzf,'D','F')/2.,'C',+1)
trans_visc_vv = fold(visc*ddz(dv2dz1_f[:],dzc,dzf,'D','F')/2.,'C',+1)
trans_vv = trans_visc_vv[:]+trans_turb_vv[:]
tot_vv = diss_vv[:]+dist_vv[:]+trans_vv[:]
fname = resultsdir + 'stats-single-point-chan-vv-budget-' + casename + fname_ext
np.savetxt(fname,np.c_[zc,tot_vv,diss_vv,dist_vv,trans_turb_vv,trans_visc_vv])
#
# ww budget
#
idw2dz1_f  = 17
iwww_f     = 18
iwp_f      = 19
ipdwdz_c   = 20
idiss_ww_c = 21
idiss_ww_x_f = 37
idiss_ww_y_f = 38
idiss_ww_z_c = 39
#
dw2dz1_f  = np.zeros(nz)
www_f     = np.zeros(nz)
p_f       = np.zeros(nz)
wp_f      = np.zeros(nz)
pdwdz_c   = np.zeros(nz)
diss_ww_c = np.zeros(nz)
if(use_new_arrays):
    diss_ww_x_f = np.zeros(nz)
    diss_ww_y_f = np.zeros(nz)
    diss_ww_z_c = np.zeros(nz)
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_mid+fname_ext
    data = np.loadtxt(fname)
    #
    dw2dz1_f [:]  += data[:,idw2dz1_f ]
    www_f    [:]  += data[:,iwww_f    ]
    wp_f     [:]  += data[:,iwp_f     ]
    pdwdz_c  [:]  += data[:,ipdwdz_c  ]
    diss_ww_c[:]  += data[:,idiss_ww_c]
    if(use_new_arrays):
        diss_ww_x_f[:] += data[:,idiss_ww_x_f]
        diss_ww_y_f[:] += data[:,idiss_ww_y_f]
        diss_ww_z_c[:] += data[:,idiss_ww_z_c]
dw2dz1_f [:]  = dw2dz1_f[:] *norm
www_f    [:]  = www_f[:]    *norm
p_f      [:]  = p_f[:]      *norm
wp_f     [:]  = wp_f[:]     *norm
pdwdz_c  [:]  = pdwdz_c[:]  *norm
diss_ww_c[:]  = diss_ww_c[:]*norm
if(use_new_arrays):
    diss_ww_x_f[:] = diss_ww_x_f[:]*norm
    diss_ww_y_f[:] = diss_ww_y_f[:]*norm
    diss_ww_z_c[:] = diss_ww_z_c[:]*norm
    diss_ww_x_c = interp(diss_ww_x_f[:],'N','F')
    diss_ww_y_c = interp(diss_ww_y_f[:],'N','F')
#
diss_ww  = fold(-visc*(diss_ww_c[:]),'C',+1)
if(use_new_arrays):
    diss_ww  = fold(-visc*(diss_ww_x_c[:]+diss_ww_y_c[:]+diss_ww_z_c[:]),'C',+1)
dist_ww  = fold(pdwdz_c[:],'C',+1)
trans_turb_ww = fold(-ddz(www_f[:],dzc,dzf,'D','F')/2.,'C',+1)
trans_pres_ww = fold(-ddz(wp_f[:],dzc,dzf,'D','F'),'C',+1)
trans_visc_ww = fold(visc*ddz(dw2dz1_f[:],dzc,dzf,'N','F')/2.,'C',+1)
trans_ww = trans_visc_ww[:]+trans_turb_ww[:]+trans_pres_ww
tot_ww = diss_ww[:]+dist_ww[:]+trans_ww[:]
fname = resultsdir+'stats-single-point-chan-ww-budget-' + casename + fname_ext
np.savetxt(fname,np.c_[zc,tot_ww,diss_ww,dist_ww,trans_turb_ww,trans_visc_ww,trans_pres_ww])
#
# TKE budget
#
prod = prod_uu[:]
diss = diss_uu[:]+diss_vv[:]+diss_ww[:]
trans_visc = trans_visc_uu[:]+trans_visc_vv[:]+trans_visc_ww[:]
trans_turb = trans_turb_uu[:]+trans_turb_vv[:]+trans_turb_ww[:]
trans_pres =                                   trans_pres_ww[:]
tot = prod[:]+diss[:]+trans_visc[:]+trans_turb[:]+trans_pres[:]
fname = resultsdir+'stats-single-point-chan-tke-budget-' + casename + fname_ext
np.savetxt(fname,np.c_[zc,tot,prod,diss,trans_turb,trans_visc,trans_pres])
#
# uw budget
#
iww_c      = 22
iduwdz1_f  = 23
iww_f      = 24
iuww_f     = 25
ip_f       = 26
iup_f      = 27
ips_c      = 28
idiss_uw_c = 29
#
ww_c      = np.zeros(nz)
duwdz1_f  = np.zeros(nz)
ww_f      = np.zeros(nz)
uww_f     = np.zeros(nz)
p_f       = np.zeros(nz)
up_f      = np.zeros(nz)
ps_c      = np.zeros(nz)
diss_uw_c = np.zeros(nz)
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_mid+fname_ext
    data = np.loadtxt(fname)
    #
    ww_c     [:]  += data[:,iww_c     ]
    duwdz1_f [:]  += data[:,iduwdz1_f ]
    ww_f     [:]  += data[:,iww_f     ]
    uww_f    [:]  += data[:,iuww_f    ]
    p_f      [:]  += data[:,ip_f      ]
    up_f     [:]  += data[:,iup_f     ]
    ps_c     [:]  += data[:,ips_c     ]
    diss_uw_c[:]  += data[:,idiss_uw_c]
ww_c     [:]  = ww_c[:]     *norm
duwdz1_f [:]  = duwdz1_f[:] *norm
ww_f     [:]  = ww_f[:]     *norm
uww_f    [:]  = uww_f[:]    *norm
p_f      [:]  = p_f[:]      *norm
up_f     [:]  = up_f[:]     *norm
ps_c     [:]  = ps_c[:]     *norm
diss_uw_c[:]  = diss_uw_c[:]*norm
#
prod_uw  = fold(-ww_c[:]*du1dz1_c[:],'C',-1)
diss_uw  = fold(-visc*(diss_uw_c[:])*2,'C',-1)
dist_uw  = fold(ps_c[:] - du1dz1_c[:]*p_c[:],'C',-1)
trans_visc_uw = fold(visc*ddz(duwdz1_f[:],dzc,dzf,'N'),'C',-1)
trans_turb_uw = fold(-ddz(uww_f[:]-u1_f[:]*ww_f[:],dzc,dzf,'D'),'C',-1)
trans_pres_uw = fold(-ddz(up_f[:]-u1_f[:]*p_f[:],dzc,dzf,'D'),'C',-1)
trans_uw = trans_visc_uw[:]+trans_turb_uw[:]+trans_pres_uw[:]
tot = prod_uw[:]+diss_uw[:]+dist_uw[:]+trans_uw[:]
fname = resultsdir+'stats-single-point-chan-uw-budget-' + casename + fname_ext
np.savetxt(fname,np.c_[zc,tot,prod_uw,diss_uw,dist_uw,trans_turb_uw,trans_visc_uw,trans_pres_uw])
#
# save histories
#
with open(resultsdir+'histories.out', 'ab') as f:
    np.savetxt(f, np.c_[(tend-tbeg)*te/t_dns,utau,utau_s,retau,cf,uc,uu_max,uw_max])
def cmpt_std_err(arr):
    return np.std(arr[:])/np.size(arr[:])**.5
dpdx_avg = np.mean(-dpdx_arr[:])
dpdx_std = np.std(dpdx_arr[:])
dpdx_err = cmpt_std_err(dpdx_arr[:])
utau_arr = (-dpdx_arr[:]*h)**.5 # utau = sqrt(-dpdx*h)
utau_avg = np.mean(utau_arr[:])
utau_err = cmpt_std_err(utau_arr[:])
retau_arr = utau_arr*h/visc
retau_avg = np.mean(retau_arr[:])
retau_err = cmpt_std_err(retau_arr[:])
cf_arr = utau_arr**2/(ub**2/2.)
cf_avg = np.mean(cf_arr[:])
cf_err =  cmpt_std_err(cf_arr[:])
with open(resultsdir+'history_utau.out', 'ab') as f:
    np.savetxt(f, np.c_[(tend-tbeg)*te/t_dns,dpdx_avg,dpdx_err,utau_avg,utau_err,retau_avg,retau_err,cf_avg,cf_err])
