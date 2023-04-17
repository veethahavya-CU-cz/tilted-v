#!/usr/bin/python

import sys, os, shutil
import subprocess
import numpy as np
from datetime import datetime, timedelta
import pickle

sys.path.append('/opt/gwswex_multilay/lib')
from gwswex_wrapper import gwswex

from modflowapi import ModflowApi

import bmi.wrapper as bw

# subprocess.run("bash -c 'source /opt/intel/oneapi/setvars.sh --config=/opt/intel/oneapi/sv_exclude-py'", shell=True)


### set-up and initialize gwswex
with open('./common/data/mdef.pkl', 'rb') as pkl:
    mdef = pickle.load(pkl)
pkl.close()
with open('./common/data/ugrid.pkl', 'rb') as pkl:
    ugrid_dict = pickle.load(pkl)
pkl.close()

tstart = datetime.strptime(mdef['model']['domain']['tstart'], '%Y%m%d %H%M%S')
tstop = datetime.strptime(mdef['model']['domain']['tstop'], '%Y%m%d %H%M%S')
ne = int(mdef['model']['domain']['nelements'])
nlay = mdef['model']['domain']['nlay']
Gdt = mdef['model']['domain']['dt']
Gnts = mdef['auxiliary']['Gnts']

ic = np.load('./common/data/ic.npz')
gw_ini = np.array(ic['gw'], dtype=np.float64, order='F')
sw_ini = np.array(ic['sw'], dtype=np.float64, order='F')

gwswex.init('./gwswex/tilted-v.yml', gw_ini, sw_ini) #NOTE: also solves 1 ts



### set-up and initialize MF6
lib_name = '/opt/mf6/lib/libmf6.so'
mf6 = ModflowApi(lib_name)
mf6.working_directory = os.path.abspath('./mf6/tilted-v/')
mf6.initialize()
mf6_ts = 0



### set-up and initialize DFM
dfm = bw.BMIWrapper(engine='dflowfm', configfile='./dfm/tilted-v.mdu')
dfm.initialize()
Lnts = gwswex.get_lnts()
assert Lnts == 1, 'Lnts != 1'


### first coupled run
curr_time = datetime.utcfromtimestamp(gwswex.get_curr_time_unix())
print("Solving: ", curr_time - timedelta(seconds=Gdt))

gw_dis = np.empty((ne, Lnts), dtype=np.float64, order='F')
sw_dis = np.empty((ne, Lnts), dtype=np.float64, order='F')
gwswex.grab_result('gw_dis_l', gw_dis)
gwswex.grab_result('sw_dis_l', sw_dis)
gw_dis = gw_dis/Gdt
sw_dis = sw_dis/Gdt
sw_dis = np.multiply(sw_dis[:,0], ugrid_dict['ea'])


# check ts fidelity
gwswex_ts = (gwswex.get_gts() - 2) * Gdt
mf6_ts = mf6.get_current_time()
dfm_ts = dfm.get_current_time()
assert int(mf6_ts) == int(dfm_ts), "Infedility b/w mf6_ts ({}) and dfm_ts ({}).".format(int(mf6_ts), int(dfm_ts))
assert int(mf6_ts) == int(gwswex_ts), "Infedility b/w mf6_ts ({}) and gwswex_time ({}).".format(int(mf6_ts), int(gwswex_ts))


# update lat
dfm.get_var('qplat')[:] = sw_dis.ravel()
# solve dfm ts
dfm.update()


# update rch
address = ['BOUND', 'TILTED-V', 'RCH_0']
rch_tag = mf6.get_var_address(*address)
mf6.get_value_ptr(rch_tag)[:] = gw_dis
# solve mf6 ts
mf6.prepare_time_step(mf6_ts)
mf6.prepare_solve(1)
address = ["MXITER", "SLN_1"]
mxittag = mf6.get_var_address(*address)
mxit = mf6.get_value_ptr(mxittag)
mxit = 100
kiter = 0
while kiter < mxit:
    has_converged = mf6.solve(1)
    kiter += 1
    if has_converged:
        break
if not has_converged:
    print("Solver did not converge")
mf6.finalize_solve(1)
mf6.finalize_time_step()



### iterative coupling
while curr_time < tstop:
    curr_time = datetime.utcfromtimestamp(gwswex.get_curr_time_unix())
    print("Solving: ", curr_time - timedelta(seconds=Gdt), "ts = ", gwswex.get_gts())

    ## update gwswex ini, solve ts, and fetch (gw_dis, sw_dis)
    gwswex.update_ini(gw_ini=gw_ini, sw_ini=sw_ini, auto_advance=True)
    Lnts = gwswex.get_lnts()
    assert Lnts == 1, 'Lnts != 1'
    # fetch gwswex prescribed discharges
    gw_dis = np.empty((ne, Lnts), dtype=np.float64, order='F')
    sw_dis = np.empty((ne, Lnts), dtype=np.float64, order='F')
    gwswex.grab_result('gw_dis_l', gw_dis)
    gwswex.grab_result('sw_dis_l', sw_dis)
    gw_dis = gw_dis/Gdt
    sw_dis = sw_dis/Gdt
    sw_dis = np.multiply(sw_dis[:,0], ugrid_dict['ea'])


    # check ts fidelity
    gwswex_ts = (gwswex.get_gts() - 2) * Gdt
    mf6_ts = mf6.get_current_time()
    dfm_ts = dfm.get_current_time()
    assert int(mf6_ts) == int(dfm_ts), "Infedility b/w mf6_ts ({}) and dfm_ts ({}).".format(int(mf6_ts), int(dfm_ts))
    assert int(mf6_ts) == int(gwswex_ts), "Infedility b/w mf6_ts ({}) and gwswex_time ({}).".format(int(mf6_ts), int(gwswex_ts))

    ## dfm
    # update lat
    dfm.get_var('qplat')[:] = sw_dis.ravel()
    # solve dfm ts
    dfm.update()


    ## mf6
    # update rch
    address = ['BOUND', 'TILTED-V', 'RCH_0']
    rch_tag = mf6.get_var_address(*address)
    mf6.get_value_ptr(rch_tag)[:] = gw_dis
    # solve mf6 ts
    mf6.prepare_time_step(mf6_ts)
    mf6.prepare_solve(1)
    address = ["MXITER", "SLN_1"]
    mxittag = mf6.get_var_address(*address)
    mxit = mf6.get_value_ptr(mxittag)
    mxit = 100
    kiter = 0
    while kiter < mxit:
        has_converged = mf6.solve(1)
        kiter += 1
        if has_converged:
            break
    if not has_converged:
        print("Solver did not converge")
    mf6.finalize_solve(1)
    mf6.finalize_time_step()

    if gwswex.get_gts() == Gnts+1:
        break



gwswex_out = {}
gwswex_out['Gstorages'] = {}
gwswex_out['Gdischarges'] = {}

gws = np.empty((ne, Gnts+1), dtype=np.float64, order='F')
sws = np.empty((ne, Gnts+1), dtype=np.float64, order='F')
sms = np.empty((nlay, ne, Gnts+1), dtype=np.float64, order='F')
epv_sm = np.empty((nlay, ne, Gnts+1), dtype=np.float64, order='F')
gwswex.pass_vars_nlay(gws, sws, sms, epv_sm)
gwswex_out['Gstorages']['gw'], gwswex_out['Gstorages']['sw'], gwswex_out['Gstorages']['sm'], gwswex_out['Gstorages']['epv_sm'] = gws, sws, sms, epv_sm

uzs = np.empty((ne, Gnts+1), dtype=np.float64, order='F')
epv_uz = np.empty((ne, Gnts+1), dtype=np.float64, order='F')
gwswex.pass_vars(gws, sws, uzs, epv_uz)
gwswex_out['Gstorages']['uzs'], gwswex_out['Gstorages']['epv_uz'] = uzs, epv_uz

gw_dis = np.empty(gws.shape, dtype=np.float64, order='F')
sw_dis = np.empty(gws.shape, dtype=np.float64, order='F')
uz_dis = np.empty(gws.shape, dtype=np.float64, order='F')
qdiff = np.empty(gws.shape, dtype=np.float64, order='F')
gwswex.pass_dis(gw_dis, uz_dis, sw_dis, qdiff)
gwswex_out['Gdischarges']['gw'], gwswex_out['Gdischarges']['sw'], gwswex_out['Gdischarges']['uz'], gwswex_out['Gdischarges']['qdiff'] = gw_dis, uz_dis, sw_dis, qdiff

with open('./gwswex/output/coupled_run.pkl', 'wb') as file:
    pickle.dump(gwswex_out, file, pickle.HIGHEST_PROTOCOL)


gwswex.finalize()
mf6.finalize()
dfm.finalize()