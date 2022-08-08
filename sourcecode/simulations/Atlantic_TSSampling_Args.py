import netCDF4 as nc
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
import pandas as pd
from parcels.tools.statuscodes import ErrorCode
import numpy as np
from kernels.samplefields import SampleTSFields
import sys


def delete_particle(particle, fieldset, time):
    particle.delete()

simulation_dt = 30 #(days: for how long to run the simulation)
    
# verify input parameters: arg1= Year, arg2= Month, arg3=Startday,arg4= simulation depth
args = sys.argv
assert len(args) == 5

start_year = np.int32(args[1])
start_mon = np.int32(args[2])
start_day = np.int32(args[3])

# verify depth argument and assign the indices to load from the depth dimension
r_depth = np.int32(args[4])
if r_depth == 0:
    min_ind, max_ind = 0, 1
elif r_depth == 50:
    min_ind, max_ind = 22, 23
elif r_depth == 100:
    min_ind, max_ind = 29, 30
elif r_depth == 200:
    min_ind, max_ind = 36, 37
elif r_depth == 500:
    min_ind, max_ind = 47, 48
else:
    raise ValueError('Depth indices have not been setup.')

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(start_year, start_mon, start_day, 12, 0, 0)
days=[simulation_start+timedelta(days=i) for i in range(simulation_dt+1)]


ufiles = [data_path + 'ROMEO.01_1d_uo_{0}{1}{2}_grid_U.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
vfiles = [data_path + 'ROMEO.01_1d_vo_{0}{1}{2}_grid_V.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
tfiles = [data_path + 'ROMEO.01_1d_thetao_{0}{1}{2}_grid_T.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
sfiles = [data_path + 'ROMEO.01_1d_so_{0}{1}{2}_grid_T.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
assert len(ufiles)==simulation_dt+1
print(ufiles[0],ufiles[-1])

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': vfiles},
             'T': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': tfiles},
             'S': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': sfiles}
             }

variables = {'U': 'uo',
             'V': 'vo',
             'T': 'thetao',
             'S': 'so'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'}

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': [min_ind, max_ind]}, chunksize=False)

coords = np.load('/nethome/manra003/analysis/paper01/H3_Res5_release_points.npz')


class Particle(JITParticle):
    min_temp = Variable('min_temp', initial=999.0)
    max_temp = Variable('max_temp', initial=-999.0)
    min_sal = Variable('min_sal', initial=999.0)
    max_sal = Variable('max_sal', initial=-999.0)


if r_depth == 0:
    depth_arg = None
else:
    depth_arg = [r_depth for i in range(len(coords['Longitude']))]

pset = ParticleSet.from_list(fieldset=fieldset,
                             pclass=Particle,
                             lon=coords['Longitude'],
                             lat=coords['Latitude'],
                             depth=depth_arg,
                             time=simulation_start)

mon = simulation_start.strftime("%b")
output_file = pset.ParticleFile(
    name="/nethome/manra003/sim_out/tara{3}m/FullTara_Res5_TS_{0}{1}{2}_dt600_z{3}.zarr".format(start_day, mon,
                                                                                              start_year,
                                                                                              r_depth),
#     outputdt=timedelta((simulation_end - simulation_start).days))
    outputdt=timedelta(days=5))

sample_kernel = pset.Kernel(SampleTSFields)

pset.execute(AdvectionRK4 + sample_kernel,
             runtime=timedelta(days=simulation_dt),
             dt=600,
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: delete_particle})

output_file.close()
