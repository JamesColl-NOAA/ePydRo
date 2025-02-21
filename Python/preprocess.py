#!/usr/bin/env python3

import re
import shutil
from os.path import basename, dirname
from pathlib import Path
import pandas as pd
import pyproj
import numpy as np
import geopandas as gpd
USE_PYGEOS=1
import h5py
from shapely.geometry import LineString

class Preprocess:
    
    
    def __init__(self):

        self.label = "Clean and preprocess HEC-RAS data"
    

    def check_for_missing_hdf_files(self,working_dir,model_folder,models_missing_proj):

        ''' Check HEC-RAS model folder for models without HDF files '''
        
        hdf_folder = working_dir / 'hdf_data'
        hdf_folder.mkdir(parents=True,exist_ok=True)

        # Find all hdf files in base directory
        hdf_files = []
        for path in model_folder.rglob("*.g[0-9][0-9].hdf"):
            hdf_files.append(path.as_posix())
        
        # Find all g## files in base directory
        g_files = []
        for path in model_folder.rglob("*.g[0-9][0-9]"):
            print(path.name)
            g_files.append(path.as_posix())

        hdf_matches = []
        g_matches = []
        no_match = []
        for file in g_files:

            # Ignore models with known missing projections
            if not any(x in file for x in models_missing_proj):
                hdf_file = Path(file + ".hdf")
                
                if hdf_file.exists():
                    hdf_matches.append(hdf_file.as_posix())
                    g_matches.append(file)
                else:
                    no_match.append(file)

        # Save table of models that do not have HDF files
        no_match_basename = [basename(i).split(".") for i in no_match]
        no_match_modelname = [i for i,j,k in no_match_basename]
        no_match_file_ext = [f"{j}.{k}" for i,j,k in no_match_basename]
        
        if len(no_match_basename) > 0:
            status = f"missing {len(no_match_basename)} hdf files"
            no_match_table = pd.DataFrame({'model_directory':no_match_basename,'model_name':no_match_modelname,'model_extension':no_match_file_ext})
            no_match_table.to_csv(working_dir / "models_that_need_hdf_files.csv",index=False)
        else: 
            status = 'No missing HDF files'
        
        return hdf_files, status

    
    def create_initial_survey_db(self,hdf_files,model_folder,survey_db_pre_path):
        
        ''' Populate table of models that have hdf files but still need projection file paths. '''
        
        # Parse HDF file paths to get RFC, model folder, and extention
        hdf_files_dirname = [dirname(i) for i in hdf_files]
        hdf_files_no_root = [i.replace(model_folder.as_posix(),"") for i in hdf_files_dirname]
        hdf_files_split = [i.split("/") for i in hdf_files_no_root]
        
        # Get RFC
        hdf_files_rfc = [i[0] for i in hdf_files_split]
        
        # Get model folder
        hdf_files_model_folder = [i[1] for i in hdf_files_split]
        
        # Get file extention
        hdf_files_basename = [basename(i).split(".") for i in hdf_files]
        hdf_files_modelname = [i for i,_j,_k in hdf_files_basename]
        hdf_files_ext = [f"{j}.{k}" for _,j,k in hdf_files_basename]
        
        # Create dataframe and write to file
        hdf_table = pd.DataFrame({'rfc':hdf_files_rfc,'model_folder':hdf_files_model_folder,'model_name':hdf_files_modelname,'model_extension':hdf_files_ext,'shp_prj_path':'','model_directory':hdf_files_dirname})
        hdf_table.to_csv(survey_db_pre_path,index=False)
    

    def create_survey_db(self,survey_db,hdf_folder):

        ''' Generate survey database from list of HDF files '''

        # Initiate dataframe
        proj_table = pd.DataFrame(columns={'rfc': pd.Series(dtype='str'),
                                    'proj': pd.Series(dtype='str'),
                                    'proj_type': pd.Series(dtype='str'),
                                    'model_name': pd.Series(dtype='str'),
                                    'hdf_path': pd.Series(dtype='str'),
                                    'model_folder': pd.Series(dtype='str'),
                                    'data_type': pd.Series(dtype='str')
                                    })      
        
        # Concatenate hdf file name
        survey_db['hdf_file'] = survey_db[['model_name','model_extension']].apply(lambda row: '.'.join(row.values.astype(str)), axis=1)
        
        # Convert hdf file name to lowercase
        survey_db['hdf_file'] = survey_db['hdf_file'].str.lower()
        
        # Filter hdf files with no duplicate names
        no_duplicates = survey_db.drop_duplicates('hdf_file')
        
        # Add non-duplicate files to survey db
        for _,model in no_duplicates.iterrows():
            
            # Copy hdf file (keep name the same)
            model_directory = Path(model.model_directory)
            hdf_path = model_directory / model.hdf_file
            shutil.copy(hdf_path,hdf_folder)
            
            # Open shapefile and get proj
            shp_ = open(Path(model.shp_prj_path),"r").readline()
            proj_string = pyproj.CRS(shp_)
            proj_type = 'string'

            # Append model to survey db
            proj_table = proj_table.append({'rfc': model.rfc,
                                        'proj': proj_string,
                                        'proj_type': proj_type,
                                        'model_name': model.model_name,
                                        'hdf_path': model.hdf_file,
                                        'model_folder': model.model_folder,
                                        'data_type': model.data_type,
                                        },ignore_index=True)
        
        # Add duplicate files to survey db
        duplicates = survey_db[survey_db.duplicated('hdf_file')]
        counter_dict = {}
        for _,model in duplicates.iterrows():
            model_directory = Path(model.model_directory)
            hdf_path = model_directory / model.hdf_file
            
            # Keep track of how many duplicates
            key, _ = re.split(".g[0-9]",model.hdf_file)
            if key in counter_dict.keys():
                counter_dict[key] += 1
            else:
                counter_dict[key] = 1
            
            # Copy hdf file (change name in new dir)
            new_ext = f"{key}_{counter_dict[key]}"
            new_hdf_file = model.hdf_file.replace(key,new_ext)
            shutil.copy(hdf_path,hdf_folder / new_hdf_file)
            
            # Open shapefile and get proj
            shp_ = open(Path(model.shp_prj_path),"r").readline()
            proj_string = pyproj.CRS(shp_)
            proj_type = 'string'

            # Append model to survey db
            proj_table = proj_table.append({'rfc': model.rfc,
                                        'proj': proj_string,
                                        'proj_type': proj_type,
                                        'model_name': model.model_name,
                                        'hdf_path': new_hdf_file,
                                        'model_folder': model.model_folder,
                                        'data_type': model.data_type,
                                        },ignore_index=True)
        
        return proj_table

    def reindex_survey_ids(self,surveys):

        # Copy old ids
        surveys['old_xid'] = surveys['xid']
        
        # Get unique model - id pairs
        combined = surveys['xid'].astype(str) + surveys['model_name']
        combined_zip = zip(combined[1:].values,combined[:-1].values)
        
        # Calculate new unique ids
        new_xid = 1
        xid_list = [new_xid]
        for i in combined_zip:
            if i[0] != i[1]:
                new_xid += 1
            xid_list.append(new_xid)
        
        # Update id column
        surveys['xid'] = xid_list

        return surveys

    
    def generate_py_xyz_from_hdf(self,args):

        hdf_path  = args[0]
        unit      = args[1]
        from_proj = args[2]
        set_proj  = args[3]
        
        print(f"processing: {hdf_path}")
        try:
            hf = h5py.File(hdf_path, 'r')
        except:
            print(f"{hdf_path} not found.")
            return pd.DataFrame()

        if unit == 'foot':
            print(f"Units in feet.. converting to meters")
            conversion_factor = 0.3048 
            elev_units = 'meters'
        elif unit == 'meter':
            print(f"Units in meters, no conversion required")
            conversion_factor = 1.0
            elev_units = 'meters'
        else:
            print(f"Unknown unit {unit}.. skipping conversion")
            conversion_factor = 1.0
            elev_units = 'unknown'
        
        # Get xy points of the plan view of the cross section
        arr_xs_points = hf.get('Geometry/Cross Sections/Polyline Points')
        arr_xs_points = np.array(arr_xs_points)

        # Get number of points per plan view cross section
        arr_pnts_per_xs = hf.get('Geometry/Cross Sections/Polyline Parts')
        arr_pnts_per_xs = np.array(arr_pnts_per_xs)
        arr_pnts_per_xs = [i[1] for i in arr_pnts_per_xs]

        # Get attribute data of the cross section (reach, river, etc...)
        arr_xs_attrib = hf.get('Geometry/Cross Sections/Attributes')
        arr_xs_attrib = np.array(arr_xs_attrib)

        # Get number of points per cross section profile
        arr_xs_profile_num_points = hf.get('Geometry/Cross Sections/Station Elevation Info')
        arr_xs_profile_num_points = np.array(arr_xs_profile_num_points)

        # Get cross section station/ elevation values
        arr_xs_station_elev = hf.get('Geometry/Cross Sections/Station Elevation Values')
        arr_xs_station_elev = np.array(arr_xs_station_elev)

        # Get mannings n values
        arr_xs_manning_n = hf.get("Geometry/Cross Sections/Manning's n Values")
        arr_xs_manning_n = np.array(arr_xs_manning_n)
        if type(arr_xs_manning_n.tolist()) == type(None):
            arr_xs_manning_n = hf.get("Geometry/Cross Sections/Station Manning's n Values")
            arr_xs_manning_n = np.array(arr_xs_manning_n)
            if type(arr_xs_manning_n.tolist()) == type(None):
                print(f"Error: no mannings n values")

        # Get mannings n info
        arr_xs_manning_info = hf.get("Geometry/Cross Sections/Manning's n Info")
        arr_xs_manning_info = np.array(arr_xs_manning_info)
        if type(arr_xs_manning_info.tolist()) == type(None):
            arr_xs_manning_info = hf.get("Geometry/Cross Sections/Station Manning's n Info")
            arr_xs_manning_info = np.array(arr_xs_manning_info)
            if type(arr_xs_manning_info.tolist()) == type(None):
                print(f"Error: no mannings n info")

        hf.close()
        
        int_start_xs_pnt = 0
        xs_count = 1
        prof_pnts_in_xs = 0
        prof_n_pnts_in_xs = 0
        
        geom_interp_pnt = []
        flt_elev = []
        station_dist = []
        xid_list = []
        mannings_list = []
        model_name_list = []
        # For each survey
        for i,v in enumerate(arr_pnts_per_xs):
            
            # Get a list of the plan cross section points
            int_end_xs_pnt = int_start_xs_pnt + v - 1
            
            # Extract survey line coords
            list_xs_points = list(map(tuple, arr_xs_points[int_start_xs_pnt:(int_end_xs_pnt + 1)]))
            
            # Convert survey line coords to linestring
            geom_xs_linestring = LineString(list_xs_points)
            
            # Update index
            int_start_xs_pnt = int_end_xs_pnt + 1
            
            # Get a list of the station - elevation points
            prof_xs_start_pnt = arr_xs_profile_num_points[i][0]
            prof_pnts_in_xs += arr_xs_profile_num_points[i][1]
            list_xs_station = arr_xs_station_elev[prof_xs_start_pnt:prof_pnts_in_xs,0]
            list_xs_elevation = arr_xs_station_elev[prof_xs_start_pnt:prof_pnts_in_xs,1]

            # Get a list of Mannings n values and station distances
            prof_xs_n_start_pnt = arr_xs_manning_info[i][0]
            prof_n_pnts_in_xs += arr_xs_manning_info[i][1]
            list_xs_n_station = arr_xs_manning_n[prof_xs_n_start_pnt:prof_n_pnts_in_xs,0]
            list_xs_n = arr_xs_manning_n[prof_xs_n_start_pnt:prof_n_pnts_in_xs,1]

            # If Manning's n values is unrealistic, skip survey
            if (list_xs_n > 1.0).sum() >= 1:
                continue

            # Bin station points by manning's n 
            xs_n_bins = np.digitize(list_xs_station, bins=list_xs_n_station[1:])
            mannings_list.extend(list_xs_n[xs_n_bins])
            
            # Rebase point distances if start is < 0 
            if list_xs_station[0] < 0.0:
                list_xs_station = [item - list_xs_station[0] for item in list_xs_station]

            # Convert units
            list_xs_station = [conversion_factor * x for x in list_xs_station]
            list_xs_elevation = [conversion_factor * x for x in list_xs_elevation]
            
            # Collect elevations
            flt_elev.extend(list_xs_elevation)

            # Collect station number (xid)
            xid_list.extend([xs_count]*len(list_xs_station))

            # Collect model name
            model_name_list.extend([hdf_path.name]*len(list_xs_station))
            
            # Normalize xs_station distances
            max_sta = np.max(list_xs_station)
            min_sta = np.min(list_xs_station)
            norm_list_xs_station = []

            norm_list_xs_station = [(i - min_sta)*geom_xs_linestring.length / (max_sta - min_sta) for i in list_xs_station]

            # Collect point distance (xid_d)
            list_xs_station_rel = np.divide(norm_list_xs_station, geom_xs_linestring.length)
            station_dist.extend(list_xs_station_rel)

            # Collect point geometries
            geom_interp_pnt.extend([geom_xs_linestring.interpolate(sta) for sta in norm_list_xs_station])
            
            xs_count = xs_count + 1
                
        # Create GeoDataFrame
        sta_elev_pnts = gpd.GeoDataFrame({
                                        'model_name': model_name_list,
                                        'xid':xid_list,
                                        'xid_d':station_dist,
                                        'geometry':geom_interp_pnt,
                                        'z':flt_elev,
                                        'n':mannings_list,
                                        },crs=from_proj)

        # set projection from input value
        sta_elev_pnts = sta_elev_pnts.to_crs(set_proj)
        sta_elev_pnts['x'] = sta_elev_pnts.geometry.x
        sta_elev_pnts['y'] = sta_elev_pnts.geometry.y
        sta_elev_pnts['source'] = 2
        sta_elev_pnts['elev_units'] = elev_units
        
        sta_elev_pnts = pd.DataFrame(sta_elev_pnts.drop(columns='geometry'))

        return sta_elev_pnts


