#!/usr/bin/python
import os
import numpy as np
import flopy as fp
import pickle
from scipy.io import FortranFile



def fwrite(fpath, val):
	Ffile = FortranFile(fpath, 'w')
	Ffile.write_record(val.T)
	Ffile.close()

def fread(fpath, shape):
	Ffile = FortranFile(fpath, 'r')
	val = Ffile.read_reals().reshape(shape, order='F')
	Ffile.close()
	return val



with open('../../common/data/mdef.pkl', 'rb') as pkl:
    mdef = pickle.load(pkl)


ne = mdef['model']['domain']['nelements']
nlay = 1 # purposeful infedility b/w GWSWEX and MF6 as GWSWEX layers are virtual (i.e. for numerical presision gain alone)

tstart = mdef['auxiliary']['tstart_dt']
tstop = mdef['auxiliary']['tstop_dt']
Gnts = int(mdef['auxiliary']['Gnts'])
Gdt = mdef['model']['domain']['dt']

ip_path = mdef['paths']['dirs']['input']
ks_path = os.path.join(ip_path, mdef['model']['domain']['layer1']['ks'])
ks = fread(ks_path, (1,ne))
gw_ini_path = os.path.join(ip_path, mdef['model']['initial conditions']['GW'])
gw_ini = fread(gw_ini_path, (1,ne))

sim = fp.mf6.MFSimulation.load(sim_name='tilted-v', version='mf6', sim_ws='../tilted-v/')
tdis = sim.tdis
ims = sim.ims
gwf = sim.get_model('tilted-v')
disv = gwf.disv

sim.remove_package('tdis')
period_data = [[Gdt, 1, 1.0]]*Gnts
tdis = fp.mf6.ModflowTdis(sim, time_units="seconds", nper=Gnts, perioddata=period_data, start_date_time=tstart.isoformat())

ict = np.ones((nlay,ne))
npf = fp.mf6.ModflowGwfnpf(gwf, icelltype=ict, k=ks, save_flows=True)
ic = fp.mf6.ModflowGwfic(gwf, strt=gw_ini)

ss = np.full((nlay,ne), 1e-5)
icv = np.ones((nlay,ne))
sto = fp.mf6.ModflowGwfsto(gwf, save_flows=True, iconvert=icv, ss=ss, sy=ss)

rch = fp.mf6.ModflowGwfrch(gwf, maxbound=ne, stress_period_data=0, save_flows=True, print_input=True)

saverec = [('HEAD', 'LAST'), ('BUDGET', 'LAST')]
printrec = [('HEAD', 'LAST')]
oc = fp.mf6.ModflowGwfoc(gwf, saverecord=saverec, printrecord=printrec, head_filerecord='heads.hds', budget_filerecord='budget.cbb')

sim.write_simulation() #silent=True