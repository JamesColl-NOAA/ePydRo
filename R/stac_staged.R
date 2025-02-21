# # https://owp-spatial.r-universe.dev/elevationSources/data/catalog_table/json
# format_for_owpstac <- function(path_to_domain) {
#   # path_to_domain <- "c:/Users/james.coll/Desktop/"
#
#   ## Start
#   if(!file.exists(file.path(path_to_domain,"ePydRo_tracker.csv",fsep = .Platform$file.sep))) {
#     message("Alert, that's not a formatted domain folder")
#     stop()
#   }
#   core_table <- data.table::fread(file.path(path_to_domain,"ePydRo_tracker.csv",fsep = .Platform$file.sep))
#
#   # Rename to match owpstac_names
#   data.table::setnames(core_table, c("sdsfeature", "V2"), c("region","mouksou"))
#
#   # Add hardocde colimns
#   core_table[,domain:="domain"]
#   core_table[,region:=core_table$sdsfeature]
#   core_table[,resolution:="10"]
#   core_table[,has_topo:="False"]
#   core_table[,has_bathymetry:="True"]
#   core_table[,priority=3]
#
#   # Remove extra columns by name
#   core_table[,surveyjobi:=NULL]
#   core_table[,sdsid:=NULL]
#   core_table[,sdsmetadat:=NULL]
#   core_table[,surveytype:=NULL]
#   core_table[,channelare:=NULL]
#   core_table[,dateupload:=NULL]
#   core_table[,usacedistr:=NULL]
#   core_table[,surveydate:=NULL]
#   core_table[,surveyda_1:=NULL]
#   core_table[,sourcedata:=NULL]
#   core_table[,sourceproj:=NULL]
#   core_table[,mediaidfk:=NULL]
#   core_table[,projecteda:=NULL]
#   core_table[,sdsfeatu_1:=NULL]
#   core_table[,dateloaded:=NULL]
#   core_table[,datenotifi:=NULL]
#   core_table[,sourceda_1:=NULL]
#   core_table[,plotsheetl:=NULL]
#   core_table[,sourceagen:=NULL]
#   core_table[,globalid:=NULL]
#   core_table[,horizontal_crs:=NULL] # NAD83(2011)
#   core_table[,vertical_datum:=NULL] # NAVD88 height
#   core_table[,vertical_datum_offset:=NULL]
#   core_table[,pos_z_convention:=NULL]
#   core_table[,vertical_datum_conversion:=NULL]
#   core_table[,epydro_source:=NULL]
#   core_table[,epydro_entwine:=NULL]
#   core_table[,epydro_domain:=NULL]
#
#   core_table[,epydro_url:=NULL]
#   core_table[,notes:=NULL]
#   core_table[,region:=NULL]
#
#   data.table::setcolorder(core_table, c("domain","region","source","resolution","has_topo","has_bathymetry","horizontal_crs","vertical_datum","vertical_datum_conversion","priority","source_url","asset_urls"))
#
#   data.table::fwrite(core_table,file.path(path_to_domain,"OWP_stac_catalog_form.csv",fsep = .Platform$file.sep))
#   return(TRUE)
# }
#
#
# ### STAC
# ```{python}
# #| warning: false
# #| echo: false
# #| eval: false
#
# import os
# import datetime
# import argparse
# import re
# import pandas
# import geopandas
# import shapely
# import pystac
#
# feature_collection = geopandas.read_file("/media/sf_G_DRIVE/data/raw/eHydro/eHydro_Survey_Data_-8591331731584774235/SurveyJob.shp")
# feature_collection = feature_collection[feature_collection['surveytype'] == 'CS']
# feature_collection = feature_collection[feature_collection.geometry.geom_type == 'Polygon']
# feature_collection = feature_collection[feature_collection.geometry.geom_type != None]
# feature_collection.dropna(subset=['geometry'], inplace=True)
#
# norm_path = "s3://lynker-spatial/bathymetry/eHydro/ePydRo/"
# catalog = pystac.Catalog(id="ePydRo", description="eHydro, nomalized and ARCO'd, in STAC form")
# for x in range(len(feature_collection.index)):
#   item = pystac.Item(
#     id=feature_collection.iloc[x]["surveyjobi"],
#     geometry=shapely.geometry.mapping(feature_collection.iloc[x]["geometry"]),
#     bbox=[feature_collection.iloc[x]["geometry"].bounds[2],
#           feature_collection.iloc[x]["geometry"].bounds[1],
#           feature_collection.iloc[x]["geometry"].bounds[0],
#           feature_collection.iloc[x]["geometry"].bounds[3]],
#     datetime=pandas.to_datetime(feature_collection.iloc[x]['dateloaded'], utc=True),
#     properties={
#       "reach name": feature_collection.iloc[x]["sdsfeature"],
#       "survey type": feature_collection.iloc[x]["surveytype"],
#       "USACE District": feature_collection.iloc[x]["usacedistr"],
#       "source data": feature_collection.iloc[x]["sourcedata"],
#       "source project": feature_collection.iloc[x]["sourceproj"],
#       "plot sheet": feature_collection.iloc[x]["plotsheetl"],
#       "globalid": feature_collection.iloc[x]["globalid"],
#       "data owner": "eHydro",
#       "datum applied": "FALSE",
#       "datescrape": "10-12-2024",
#     },
#   )
#
# item.add_asset(
#   key=feature_collection.iloc[x]["surveyjobi"],
#   asset=pystac.Asset(href=os.path.join(norm_path,feature_collection.iloc[x]["surveyjobi"],"ePydRo.fgb"), media_type=pystac.MediaType.FLATGEOBUF))
# catalog.add_item(item)
#
# #if norm_path != 'NONE':
# #    exit
#
# catalog.normalize_hrefs("/media/sf_G_DRIVE/data/raw/eHydro/eHydro_Survey_Data_-8591331731584774235/stac")
# catalog.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)
#
# ### Push to AWS
# `aws s3 sync G:\data\raw\eHydro\eHydro_Survey_Data_-8591331731584774235\ePydRo\ s3://lynker-spatial/bathymetry/eHydro/ePydRo/`
