import cdsapi

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "mean_sea_level_pressure"
    ],
    "year": ["2024"],
    "month": ["06", "07"],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "25", "26", "27",
        "28", "29", "30"
    ],
    "time": [
        "00:00", "01:00", "02:00",
        "03:00", "04:00", "05:00",
        "06:00", "07:00", "08:00",
        "09:00", "10:00", "11:00",
        "12:00", "13:00", "14:00",
        "15:00", "16:00", "17:00",
        "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [35, -100, 5, -10]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download("/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/ERA5/era5_sfc_wind_press_2024.nc")

