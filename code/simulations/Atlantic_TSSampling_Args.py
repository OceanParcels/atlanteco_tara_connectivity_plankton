import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
import pandas as pd
from parcels.tools.statuscodes import ErrorCode
import numpy as np
from kernels.samplefields import SampleTSFields
import sys


def delete_particle(particle, fieldset, time):
    particle.delete()


# verify input parameters: arg1= Year, arg2= Month, arg3=Startday
args = sys.argv
assert len(args) == 4

start_year = np.int32(args[1])
assert 2009 <= start_year <= 2018

start_mon = np.int32(args[2])
assert 1 <= start_mon <= 12

start_day = np.int32(args[3])
assert 1 <= start_day <= 31

# by default- changes only in case od 1-12-2018
end_day = start_day

if start_mon == 12:
    if start_year == 2018:
        raise ValueError('data unavailable for complete simulation.')
    end_mon = 1
    end_year = start_year + 1
else:
    end_mon = start_mon + 1
    end_year = start_year

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(start_year, start_mon, start_day, 12, 0, 0)
simulation_end = datetime(end_year, end_mon, end_day, 12, 0, 0)
r_depth = 100

ufiles = sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'. \
                     format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'. \
                     format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

vfiles = sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'. \
                     format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'. \
                     format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

tfiles = sorted(glob(data_path + 'ROMEO.01_1d_thetao_{0}{1}*_T.nc'. \
                     format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                glob(data_path + 'ROMEO.01_1d_thetao_{0}{1}*_T.nc'. \
                     format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

sfiles = sorted(glob(data_path + 'ROMEO.01_1d_so_{0}{1}*_T.nc'. \
                     format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                glob(data_path + 'ROMEO.01_1d_so_{0}{1}*_T.nc'. \
                     format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

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

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles) - 1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': [29, 30]}, chunksize=False)

coords = pd.read_csv(r'/nethome/manra003/data/Nemo_H3Release_LatLon_Res5.csv')


class Particle(JITParticle):
    min_temp = Variable('min_temp', initial=999.0)
    max_temp = Variable('max_temp', initial=-999.0)
    min_sal = Variable('min_sal', initial=999.0)
    max_sal = Variable('max_sal', initial=-999.0)


pset = ParticleSet.from_list(fieldset=fieldset,
                             pclass=Particle,
                             lon=coords['Longitudes'],
                             lat=coords['Latitudes'],
                             time=simulation_start,
                             depth=[r_depth for i in range(len(coords))])

mon = simulation_start.strftime("%b")
output_file = pset.ParticleFile(
    name="/nethome/manra003/sim_out/tara{3}m/FullTara_Res5_TS_{0}{1}{2}_dt600_z{3}.nc".format(start_day, mon,
                                                                                              start_year,
                                                                                              r_depth),
    outputdt=timedelta((simulation_end - simulation_start).days))

sample_kernel = pset.Kernel(SampleTSFields)

pset.execute(AdvectionRK4 + sample_kernel,
             endtime=simulation_end,
             dt=600,
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: delete_particle})

output_file.close()
