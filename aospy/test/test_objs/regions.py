from aospy.region import Region


# No land mask, simple lat/lon bounds
nh = Region(
    name='nh',
    description='Northern Hemisphere',
    lat_bounds=(0, 90),
    lon_bounds=(0, 360),
    do_land_mask=False
)

# Land mask, more complicated lat/lon bounds
sahel = Region(
    name='sahel',
    description='African Sahel',
    mask_bounds=[((10, 20), (0, 40)), ((10, 20), (342, 360))],
    do_land_mask=True
)

# Ocean mask, simple lat/lon bounds
nh_ocean = Region(
    name='nh_ocean',
    description='Northern Hemisphere Ocean',
    lat_bounds=(0, 90),
    lon_bounds=(0, 360),
    do_land_mask='ocean'
)
