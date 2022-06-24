import numpy as np
from glob import glob
import xarray as xr

data_folder = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
output_folder = '/nethome/manra003/analysis/paper01/avg_CMCC/'


def month_average(mon, field, file_pattern):
    # get all the files for a given month

    files = glob(data_folder + file_pattern.format(mon))

    temp_ds = xr.open_mfdataset(files)
    if field == 't':
        month_avg = temp_ds.thetao.isel(deptht=0).mean('time_counter')
    elif field == 'u':
        month_avg = temp_ds.uo.isel(deptht=0).mean('time_counter')
    elif field == 'v':
        month_avg = temp_ds.vo.isel(deptht=0).mean('time_counter')
    else:
        raise ValueError('incorrect field value')
    month_avg.to_netcdf(output_folder + '{0}_{1}.nc'.format(field, mon))


# ONLY WORKING WITH SURFACE CURRENTS AT THE MOMENT- deptht=0

# region: Monthly average of fields
months = np.array(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
file_pattern = 'ROMEO.01_1d_thetao_20*{0}*_T.nc'
[month_average(mon, 't', 'ROMEO.01_1d_thetao_20*{0}*_T.nc') for mon in months]

[month_average(mon, 'u', 'ROMEO.01_1d_uo_20*{0}*_U.nc') for mon in months]

[month_average(mon, 'v', 'ROMEO.01_1d_vo_20*{0}*_V.nc') for mon in months]
# endregion

# region: Average over all the data
