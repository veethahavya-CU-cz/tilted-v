#!/usr/bin/python
import os
import numpy as np
import pickle
import yaml
from scipy.io import FortranFile
from datetime import datetime, timedelta



def fwrite(fname, val):
	Ffile = FortranFile(os.path.join(ip_path,fname), 'w')
	Ffile.write_record(val.T)
	Ffile.close()

def fread(fname, shape):
	Ffile = FortranFile(os.path.join(op_path,fname), 'r')
	val = Ffile.read_reals().reshape(shape, order='F')
	Ffile.close()
	return val



with open('../../common/data/ugrid.pkl', 'rb') as ugrid_file:
    ugrid_dict = pickle.load(ugrid_file)
ugrid_file.close()



### ACTIVE MODEL DEFINITIONS
tstart = datetime.strptime("01.02.1995 00:00:00", '%d.%m.%Y %H:%M:%S')
tstop = tstart + timedelta(hours=120)
Gdt = timedelta(seconds=600)
Gnts = int((tstop - tstart)/Gdt)
assert Gnts%10 == 0, "choose Gdt in such a way that it produces whole number of timesteps!"


root_path = os.path.abspath('../')
ip_path = os.path.join(root_path, 'input')
if not os.path.exists(ip_path): os.makedirs(ip_path)
op_path = os.path.join(root_path, 'output')
if not os.path.exists(op_path): os.makedirs(op_path)


ne = ugrid_dict['ne']
nlay = 5
uz_thickness = 5
lay_thickness_ratios = np.array([0.1, 0.1, 0.2, 0.2, 0.4])
if not lay_thickness_ratios.sum() == 1: lay_thickness_ratios = lay_thickness_ratios/lay_thickness_ratios.sum()

top = ugrid_dict['ez']
fwrite(os.path.join(ip_path, 'top.ip'), top)

bot = np.empty((nlay,ne))
for lay in range(nlay):
    if lay == 0:
        bot[lay] = top - (lay_thickness_ratios[lay] * uz_thickness)
    else:
        bot[lay] = bot[lay-1] - (lay_thickness_ratios[lay] * uz_thickness)
fwrite(os.path.join(ip_path, 'bot.ip'), bot)

isactive = np.full(ne, True, dtype=bool, order='F')
for lay in range(nlay):
    fwrite(os.path.join(ip_path, 'l{}_active.ip'.format(lay+1)), isactive)

ks = np.full(ne, 10/3600, dtype=np.float64, order='F') # 10 m/h
for lay in range(nlay):
    fwrite(os.path.join(ip_path, 'l{}_ks.ip'.format(lay+1)), ks)

porosity = np.full(ne, 0.4, dtype=np.float64, order='F')
for lay in range(nlay):
    fwrite(os.path.join(ip_path, 'l{}_porosity.ip'.format(lay+1)), porosity)


chd = np.full(ne, 0, dtype=int, order='F')
fwrite(os.path.join(ip_path, 'GW_chd.ip'), chd)

gw_ini = np.full(ne, (top - 2), dtype=np.float64, order='F')
fwrite(os.path.join(ip_path, 'GW_ini.ip'), gw_ini)
sw_ini = np.full(ne, 0, dtype=np.float64, order='F')
fwrite(os.path.join(ip_path, 'SW_ini.ip'), sw_ini)
np.savez('../../common/data/ic.npz', gw=gw_ini, sw=sw_ini)

p = np.zeros((ne,Gnts+1), dtype=np.float64, order='F')
p[:,:120] = 0.1/3600 # 0.1 m/h
fwrite(os.path.join(ip_path, 'p.ip'), p)
et = np.full((ne,Gnts+1), 0, dtype=np.float64, order='F')
fwrite(os.path.join(ip_path, 'et.ip'), et)

vG = [6.0, 2.0, 0.08, 0.4]


mdef_tformat = '%Y%m%d %H%M%S'
mdef = {
    'model': {
        'domain': {
            'nelements': ne,
            'dt': int(Gdt.total_seconds()),
            'tstart': tstart.strftime(mdef_tformat),
            'tstop': tstop.strftime(mdef_tformat),
            'nlay': nlay,
            'top': 'top.ip',
            'bot': 'bot.ip',
            'layer1': {
                'name': 'l1',
                'isactive': 'l1_active.ip',
                'vanG': vG[:],
                'ks': 'l1_ks.ip',
                'porosity': 'l1_porosity.ip'
            },
            'layer2': {
                'name': 'l2',
                'isactive': 'l2_active.ip',
                'vanG': vG[:],
                'ks': 'l2_ks.ip',
                'porosity': 'l2_porosity.ip'
            },
            'layer3': {
                'name': 'l3',
                'isactive': 'l3_active.ip',
                'vanG': vG[:],
                'ks': 'l3_ks.ip',
                'porosity': 'l3_porosity.ip'
            },
            'layer4': {
                'name': 'l4',
                'isactive': 'l4_active.ip',
                'vanG': vG[:],
                'ks': 'l4_ks.ip',
                'porosity': 'l4_porosity.ip'
            },
            'layer5': {
                'name': 'l5',
                'isactive': 'l5_active.ip',
                'vanG': vG[:],
                'ks': 'l5_ks.ip',
                'porosity': 'l5_porosity.ip'
            }
        },
        'boundary conditions': {
            'GW CHD': 'GW_chd.ip'
        },
        'initial conditions': {
            'GW': 'GW_ini.ip',
            'SW': 'SW_ini.ip'
        },
        'external forcings': {
            'p': 'p.ip',
            'et': 'et.ip'
        },
        'solver settings': {
            'pet_intensities': [0.05/3600, 0.11/3600, 0.21/3600],
            'pet_nts': [1, 1, 1]
        }
    },
    'paths': {
        'dirs': {
            'root': root_path,
            'input': ip_path,
            'output': op_path
        },
        'files': {
            'DMN.TOP': 1,
            'DMN.BOT': 1,
            'DMN.LAY.ACT': 1,
            'DMN.LAY.KS': 1,
            'DMN.LAY.POR': 1,
            'BND.CHD': 1,
            'IC.GW': 1,
            'IC.SW': 1,
            'EXTF.p': 1,
            'EXTF.et': 1
        }
    },
    'utils': {
        'logger': {
            'level': 1,
            'fname': 'GWSWEX.log'
        }
    }
}



def represent_list(dumper, data):
    if len(data) == 1:
        return dumper.represent_scalar(u'tag:yaml.org,2002:seq', str(data[0]))
    else:
        return dumper.represent_sequence(u'tag:yaml.org,2002:seq', data, flow_style=True)

yaml.add_representer(list, represent_list)

with open('../tilted-v.yml', 'w') as outfile:
    yaml.dump(mdef, outfile, sort_keys=False)

mdef['auxiliary'] = {}
mdef['auxiliary']['mdef_tformat'] = mdef_tformat
mdef['auxiliary']['Gnts'] = Gnts
mdef['auxiliary']['tstart_dt'] = tstart
mdef['auxiliary']['tstop_dt'] = tstop

with open('../../common/data/mdef.pkl', 'wb') as file:
    pickle.dump(mdef, file)