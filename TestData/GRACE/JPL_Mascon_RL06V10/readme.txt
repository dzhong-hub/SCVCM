Driver: netCDF/Network Common Data Format
Files: GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc
Size is 512, 512
Coordinate System is `'
Metadata:
  NC_GLOBAL#acknowledgement=GRACE is a joint mission of NASA (USA) and DLR (Germany).  Use the digital object identifier provided in the id attribute when citing this data.  See https://podaac.jpl.nasa.gov/CitingPODAAC
  NC_GLOBAL#C_20_substitution=Cheng, M., Ries, and Tapley (2011), J. Geophys. Res., 116, B01409.
  NC_GLOBAL#Conventions=CF-1.6, ACDD-1.3, ISO 8601
  NC_GLOBAL#creator_email=grace@podaac.jpl.nasa.gov
  NC_GLOBAL#creator_institution=NASA/JPL
  NC_GLOBAL#creator_name=David Wiese
  NC_GLOBAL#creator_type=group
  NC_GLOBAL#creator_url=https://grace.jpl.nasa.gov
  NC_GLOBAL#date_created=2018-10-01T12:10:07Z
  NC_GLOBAL#geocenter_correction=Swenson, Chambers, and Wahr (2008), J. Geophys. Res., 113, 8410.
  NC_GLOBAL#geospatial_lat_max=89.75
  NC_GLOBAL#geospatial_lat_min=-89.75
  NC_GLOBAL#geospatial_lat_resolution=0.5 degree grid; however the native resolution of the data is 3-degree equal-area mascons
  NC_GLOBAL#geospatial_lat_units=degrees_north
  NC_GLOBAL#geospatial_lon_max=359.75
  NC_GLOBAL#geospatial_lon_min=0.25
  NC_GLOBAL#geospatial_lon_resolution=0.5 degree grid; however the native resolution of the data is 3-degree equal-area mascons
  NC_GLOBAL#geospatial_lon_units=degrees_east
  NC_GLOBAL#GIA_removed=ICE6G-D; Peltier, W. R., D. F. Argus, and R. Drummond (2018) Comment on the paper by Purcell et al. 2016 entitled An assessment of ICE-6G_C (VM5a) glacial isostatic adjustment model, J. Geophys. Res. Solid Earth, 122.
  NC_GLOBAL#id=10.5067/TEMSC-3MJ06
  NC_GLOBAL#institution=NASA/JPL
  NC_GLOBAL#journal_reference=Watkins, M. M., D. N. Wiese, D.-N. Yuan, C. Boening, and F. W. Landerer (2015) Improved methods for observing Earth's time variable mass distribution with GRACE using spherical cap mascons, J. Geophys. Res., 120, 2648-2671. 
  NC_GLOBAL#keywords=Solid Earth, Geodetics/Gravity, Gravity, liquid_water_equivalent_thickness
  NC_GLOBAL#keywords_vocabulary=NASA Global Change Master Directory (GCMD) Science Keywords
  NC_GLOBAL#license=https://science.nasa.gov/earth-science/earth-science-data/data-information-policy
  NC_GLOBAL#Metadata_Conventions=Unidata Dataset Discovery v1.0
  NC_GLOBAL#months_missing=2002-06;2002-07;2003-06;2011-01;2011-06;2012-05;2012-10;2013-03;2013-08;2013-09;2014-02;2014-07;2014-12;2015-06;2015-10;2015-11;2016-04;2016-09;2016-10;2017-02
  NC_GLOBAL#naming_authority=org.doi.dx
  NC_GLOBAL#platform=GRACE
  NC_GLOBAL#postprocess_1= A MASCON-AVERAGED REPRESENTATION OF THE OCEAN_ATMOSPHERE_DEALIAS_MODEL (GAD), MONTHLY_AVE, ADDED BACK
  NC_GLOBAL#postprocess_2=Water density used to convert to equivalent water height: 1000 kg/m^3
  NC_GLOBAL#processing_level=2 and 3
  NC_GLOBAL#product_version=v1.0
  NC_GLOBAL#program=NASA Earth Science System Pathfinder
  NC_GLOBAL#project=NASA Gravity Recovery and Climate Experiment (GRACE)
  NC_GLOBAL#publisher_email=podaac@jpl.nasa.gov
  NC_GLOBAL#publisher_institution=NASA/JPL
  NC_GLOBAL#publisher_name=Physical Oceanography Distributed Active Archive Center
  NC_GLOBAL#publisher_type=group
  NC_GLOBAL#publisher_url=https://podaac.jpl.nasa.gov
  NC_GLOBAL#source=GRACE JPL RL06M
  NC_GLOBAL#standard_name_vocabulary=NetCDF Climate and Forecast (CF) Metadata Convention-1.6
  NC_GLOBAL#summary=Monthly gravity solutions from GRACE as determined from the JPL RL06M mascon solution - CRI filter is NOT applied
  NC_GLOBAL#time_coverage_end=2017-06-11T23:59:59Z
  NC_GLOBAL#time_coverage_start=2002-04-16T00:00:00Z
  NC_GLOBAL#time_epoch=2002-01-01T00:00:00Z
  NC_GLOBAL#time_mean_removed=2004.000 to 2009.999
  NC_GLOBAL#title=JPL GRACE MASCON RL06M
  NC_GLOBAL#user_note_1=The accelerometer on the GRACE-B spacecraft was turned off after August 2016.  After this date, the accelerometer on GRACE-A was used to derive the non-gravitational accelerations acting on GRACE-B using a transplant procedure.  This has led to a subsequent degradation in the quality of the gravity fields derived.  The uncertainties in this file have been scaled to accomodate this degradation; however, we still recommend caution when using data past August 2016.
Subdatasets:
  SUBDATASET_1_NAME=NETCDF:"GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01.nc":lwe_thickness
  SUBDATASET_1_DESC=[163x360x720] Liquid_Water_Equivalent_Thickness (64-bit floating-point)
  SUBDATASET_2_NAME=NETCDF:"GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01.nc":uncertainty
  SUBDATASET_2_DESC=[163x360x720] uncertainty (64-bit floating-point)
  SUBDATASET_3_NAME=NETCDF:"GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01.nc":lat_bounds
  SUBDATASET_3_DESC=[360x2] lat_bounds (64-bit floating-point)
  SUBDATASET_4_NAME=NETCDF:"GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01.nc":lon_bounds
  SUBDATASET_4_DESC=[720x2] lon_bounds (64-bit floating-point)
  SUBDATASET_5_NAME=NETCDF:"GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01.nc":time_bounds
  SUBDATASET_5_DESC=[163x2] time_bounds (64-bit floating-point)
Corner Coordinates:
Upper Left  (    0.0,    0.0)
Lower Left  (    0.0,  512.0)
Upper Right (  512.0,    0.0)
Lower Right (  512.0,  512.0)
Center      (  256.0,  256.0)
