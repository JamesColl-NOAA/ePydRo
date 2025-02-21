#!/usr/bin/env python3

import warnings
from datetime import date
from os.path import splitext

import geopandas as gpd
import numpy as np
import pandas as pd
from pygeos import from_shapely, linear,from_wkb,to_wkb

import xarray as xr
from shapely import wkb
from shapely.geometry import LineString, MultiPoint, Point
warnings.simplefilter('ignore')


class Spatial:
    
    
    def __init__(self):

        self.label = "Convert survey/stream data to geospatial objects"


    def survey_table_to_points(self, survey_table, set_proj):

        '''Converts survey dataframe to Geopandas points layer'''

        # Convert lat lon to points
        survey_geom = [Point(xy) for xy in zip(survey_table.x, survey_table.y)]

        # Add geometry attribute to survey data
        survey_pts = gpd.GeoDataFrame(survey_table, geometry=survey_geom, crs=set_proj)

        return survey_pts


    def survey_points_to_line(self, survey_pts, id_col, set_proj):

        '''Converts survey Geopandas points layer to lines layer'''

        # Group points by id an convert to lines
        survey_lines = survey_pts.groupby([id_col])['geometry'].apply(
            lambda x: LineString(x.tolist())
        )

        survey_lines = survey_lines.reset_index()

        survey_lines = gpd.GeoDataFrame(survey_lines, crs=set_proj, geometry='geometry')

        return survey_lines

    
    def survey_lines_to_points(self, survey_lines, survey_table, id_col, set_proj):

        '''Converts survey Geopandas lines layer to points layer'''
        
        # Subset survey table
        survey_lines = survey_lines.drop(
            columns=['model_name', 'rfc', 'geometry', 'source', 'thalweg']
        )
        
        survey_table_subset = pd.merge(survey_lines, survey_table, how='left', on=id_col)
        
        # Convert survey data to points
        survey_pts = self.survey_table_to_points(survey_table_subset, set_proj)
        
        # Reset index
        survey_pts = survey_pts.reset_index(drop=True)
        
        return survey_pts

    
    def survey_table_to_lines(self, survey_table, set_proj, id_col):

        '''Converts survey dataframe to Geopandas lines layer'''

        # Convert survey data to points
        survey_pts = self.survey_table_to_points(survey_table, set_proj)

        # Convert survey points to lines
        survey_lines = self.survey_points_to_line(survey_pts, id_col, set_proj)

        return survey_lines

    
    def create_convex_hull(self, survey_lines):

        '''Creates poygon mask of survey lines using convex hull method'''

        linestrings = survey_lines.geometry
        aoi = linestrings.unary_union
        aoi = aoi.convex_hull

        return aoi

    
    def write_spatial(self, spatial_, fileName, layer=None, index=False, verbose=False):

        '''Gets driver name from file extension for Geopandas writing'''

        if verbose:
            print(f"Writing to {fileName}")

        # sets driver
        driverDictionary = {
            '.gpkg': 'GPKG',
            '.geojson': 'GeoJSON',
            '.shp': 'ESRI Shapefile',
            '.feather': 'feather'
        }
        driver = driverDictionary[splitext(fileName)[1]]
        if driver == 'feather':
            spatial_.to_feather(fileName)
        else:
            spatial_.to_file(fileName, driver=driver, layer=layer, index=index)

    
    def calculate_relative_survey_pt_dist(
        self, survey_lines, id_col, survey_table, set_proj
    ):

        '''Calculates relatve distance from thalweg (xid_d) using stream intersection'''

        survey_pts = self.survey_lines_to_points(
            survey_lines, survey_table, id_col, set_proj
        )

        if not survey_pts.crs == survey_lines.crs:
            survey_pts = survey_pts.to_crs(survey_lines.crs)     

        # Calculate thalweg_distance for survey_lines
        params = list(zip(survey_lines.geometry.to_wkb().to_numpy(),survey_lines.thalweg.to_numpy()))
        survey_lines['thalweg_distance'] = np.array([linear.line_locate_point(from_wkb(geoms),from_wkb(thalweg)) for geoms,thalweg in params])

        # Add survey geom to survey_pts
        survey_lines['line_geometry'] = survey_lines.geometry.to_wkb()
        survey_pts = survey_pts.merge(survey_lines[[id_col,'thalweg_distance','line_geometry']],on=[id_col],how='right')

        # Convert geom objects to WKB byte objects
        survey_pts['point_geometry'] = survey_pts.geometry.to_wkb()
        params = survey_pts[['point_geometry','line_geometry','thalweg_distance']].to_numpy()
        
        # Calculate xid_d_post for each point in survey
        survey_pts['xid_d_post'] = [np.round(linear.line_locate_point(from_wkb(l),from_wkb(x)) - t,2) for x,l,t in params]
        
        survey_pts.drop(
            columns=['line_geometry','thalweg_distance','point_geometry'],
            inplace=True
        )

        return survey_pts

    
    def check_data_types(self, data_frame):

        '''Enforce data types for selected columns'''

        dtype_dict = {
            'xid_length': float,
            'comid': int,
            'old_comid': int,
            'to': int,
            'old_to': int,
            'model_name': str,
            'lengthMap': str,
            'rfc': str,
            'comid_rel_dist_ds': float,
            'xid_d': float,
            'xid_d_post': float,
            'x': float,
            'y': float,
            'z': float,
            'n': float,
            'elev_units': str,
            'source': int,
            'survey_ids': str,
        }
        
        for col in data_frame.columns.unique():
            if col in dtype_dict.keys():
                data_frame[col] = data_frame[col].astype(dtype_dict[col])

        return data_frame


    def subset_streams_with_survey_intersections(self, streams, survey_lines, refactor,crosswalk_col, junction_to):

        '''Find streams with survey data (by intersection)'''

        stream_list = []
        # Find stream with cross sections and add to list
        for _, stream in streams.iterrows():

            crosses = survey_lines.crosses(stream.geometry)
            survey_lines_subset = survey_lines[crosses]

            if not survey_lines_subset.empty:

                # Collect survey line ids and stream attributes
                stream_att = {
                    'comid': stream.comid,
                    'to': stream.to,
                    'NHDWaterbodyComID': stream.NHDWaterbodyComID,
                    'order': stream.order,
                    'Length': stream.Length,
                    'gages': stream.gages,
                    'survey_ids': survey_lines_subset.xid.to_list(),
                    'geometry': stream.geometry
                }

                if refactor == 'True': 
                    # Add refactor crosswalk
                    stream_att[crosswalk_col]= stream[crosswalk_col]
                    stream_att[junction_to]= stream[junction_to]
                    stream_att['old_to'] = stream['old_to']
                    stream_att['old_comid'] = stream['old_comid']

                stream_list.append(stream_att)

        relevant_streams = gpd.GeoDataFrame(
            stream_list, geometry='geometry', crs=survey_lines.crs
        )

        return relevant_streams

    
    def remove_adjacent_trib_surveys(self, surveys, id_col, streams):

        '''Remove stream segments/surveys that intersect surveys from adjacent streams'''

        for xid_ in surveys[id_col].unique():
            survey_sel = surveys.loc[surveys.xid == xid_]
            survey_sel = survey_sel.reset_index(drop=True)
            if len(survey_sel) > 1:

                # Identify which stream is closest to the centerline
                stream_segments = survey_sel.comid.to_list()
                mult_streams = streams.loc[streams.comid.isin(stream_segments)]
                mult_streams = mult_streams.reset_index(drop=True)
                largest_order = max(mult_streams.order)
                hiorder_streams = mult_streams.loc[
                    mult_streams.order == largest_order
                ].comid

                # Remove lower order streams
                if len(hiorder_streams) == 1:
                    surveys = surveys.loc[
                        ~(
                            (surveys.xid == xid_)
                            & ~(surveys.comid == hiorder_streams.item())
                        ),
                        :,
                    ]
                else:
                    surveys = surveys.loc[
                        ~(
                            (surveys.xid == xid_)
                            & ~(surveys.comid.isin(hiorder_streams.to_list()))
                        ),
                        :,
                    ]
                    print(
                        f"{len(hiorder_streams.to_list())} streams with {largest_order} order: {hiorder_streams.to_list()}"
                    )

        return surveys

    
    def get_xs_closest_to_stream_point(self, streams, survey_lines, id_col, distance, refactor,crosswalk_col,junction_to):

        '''Find survey closest to specified distance along stream segment'''

        # Get stream point as specified distance
        stream_list = []
        for _,stream in streams.iterrows():

            relevant_surveys = survey_lines.loc[survey_lines[id_col].isin(stream.survey_ids)]
            relevant_surveys = relevant_surveys.reset_index(drop=True)

            stream_point = stream.geometry.interpolate(distance, normalized=True)

            # Distance from point on stream to each survey
            distances = relevant_surveys.distance(stream_point)

            # Find minimum distance
            min_index = np.argmin(distances)
            
            # Get closest survey
            closest_survey = relevant_surveys.loc[min_index].copy()
           
            # Intersect closest survey with stream
            intersect = stream.geometry.intersection(closest_survey.geometry)

            # if multiple intersections are found
            if type(intersect) is MultiPoint:
                print(
                    f"multiple intersections for comid {stream.comid}, RFC {closest_survey.rfc}, \
                        model {closest_survey.model_name}, xid {closest_survey.xid}"
                )
                multipoints = gpd.GeoDataFrame({'geometry': intersect}, crs=streams.crs)
                instersection_distances = multipoints.distance(stream_point)

                # Find minimum distance
                min_index = np.argmin(instersection_distances)
                intersect = multipoints.loc[min_index].geometry

            # Calculate intersection distance along stream segment
            relative_distance = linear.line_locate_point(
                from_shapely(stream.geometry), from_shapely(intersect), normalized=True
            )
            
            # Get closest survey
            stream_att = {
                'comid': stream.comid,
                'to': stream.to,
                'xid': closest_survey.xid,
                'NHDWaterbodyComID': stream.NHDWaterbodyComID,
                'order': stream.order,
                'model_name': closest_survey.model_name,
                'rfc': closest_survey.rfc,
                'source': closest_survey.source,
                'Length': stream.Length,
                'gages': stream.gages,
                'comid_rel_dist_ds': str(relative_distance),
                'thalweg': wkb.dumps(intersect),
                'geometry': closest_survey.geometry
            }

            if refactor == 'True': 
                # Add refactor crosswalk
                stream_att[crosswalk_col]= stream[crosswalk_col]
                stream_att[junction_to]= stream[junction_to]
                stream_att['old_to'] = stream['old_to']
                stream_att['old_comid'] = stream['old_comid']
            
            stream_list.append(stream_att)

        # Convert table to spatial layer 
        selected_surveys = gpd.GeoDataFrame(
            stream_list, geometry='geometry', crs=survey_lines.crs
        )

        return selected_surveys

    
    def export_surveys(self, survey_pts, data_type, set_proj, version, netcdf_path):

        '''Convert survey data to routelink netcdf'''

        # Reset index and name 'feature_id'
        survey_pts.reset_index(inplace = True,drop = True)
        survey_pts.index.name = 'feature_id'

        survey_netcdf = xr.Dataset.from_dataframe(survey_pts)

        # Add global attributes
        survey_netcdf.attrs['featureType'] = 'survey points'
        survey_netcdf.attrs['title'] = 'Natural cross section data'
        survey_netcdf.attrs['data_type'] = f"{data_type}"
        survey_netcdf.attrs['date_generated'] = f"database generated on {date.today()}"
        survey_netcdf.attrs['repository'] = "https://github.com/NOAA-OWP/inland_bathy"
        survey_netcdf.attrs['projection'] = f"{set_proj}"
        survey_netcdf.attrs['version'] = version

        # Save as final output
        survey_netcdf.to_netcdf(netcdf_path)
    

    def subset_refactored_network(self,comid_list,refactored_network_path,crosswalk_col,origin_streams,junction_to):

        '''Subset refactored network based on NWM COMIDs'''

        # Refactored network
        refactored_network = gpd.read_file(refactored_network_path)
        
        # Get list of NWM COMIDs in crosswalk
        crosswalk = list(zip(refactored_network.ID.tolist(),refactored_network[crosswalk_col].str.split(',').tolist()))
        
        # Using a list of comids, find matches in crosswalk and collect the corresponding refactored id
        refactored_ids = []
        junc_dict = {}
        for stream in crosswalk:
            
            refactored_ids += [stream[0] for s in stream[1] if int(s.split('.')[0]) in comid_list]
            origin_segs = [int(s.split('.')[0]) for s in stream[1] if int(s.split('.')[0]) in comid_list]
            
            if len(origin_segs) > 0:
                
                if len(origin_segs) > 1:
                    to_list = origin_streams.loc[origin_streams.ID.isin(origin_segs)].to.to_list()
                    upstream_seg = [x for x in origin_segs if x not in to_list]
                    
                    if len(upstream_seg) == 1:
                        junc_dict[stream[0]] = upstream_seg[0]
                    else:
                        print(f"Error: multiple upstream IDs returned for refactored stream {stream[0]}")
                
                else:
                    junc_dict[stream[0]] = origin_segs[0]
        
        # Subset refactored network
        refactored_network = refactored_network.loc[refactored_network.ID.isin(refactored_ids)]
        refactored_network = refactored_network.reset_index(drop=True)
        
        # Calculate origin stream to route tributaries/boundary streams
        refactored_network[junction_to] = refactored_network.ID.map(junc_dict)
        
        # Subset refactored stream attributes
        refactored_network = refactored_network.filter(items=['ID', crosswalk_col, junction_to, 'toID','order', 'gages', 'NHDWaterbodyComID','geometry'])
        
        # Rename fields and set data types
        refactored_network = refactored_network.rename(columns={'ID':'comid','toID':'to'})
        refactored_network.to = refactored_network.to.astype(int)
        refactored_network.order = refactored_network.order.astype(int)
        refactored_network.NHDWaterbodyComID = refactored_network.NHDWaterbodyComID.fillna(-9999)
        refactored_network.NHDWaterbodyComID = refactored_network.NHDWaterbodyComID.astype(int)
        refactored_network['Length'] = refactored_network.geometry.length

        return refactored_network


    def check_crosswalk_fractional_sum(self,refactored_network,crosswalk_col):
        
        '''Check that COMID crosswalk fractions all sum to one'''

        flag = False
        frac_sum = {}
        # Get list of NWM COMIDs in crosswalk and sum fractional crosswalks
        for _,xw in enumerate(refactored_network[crosswalk_col]):
            if ',' in xw:
                comid_list = xw.split(',')
                for xw_ in comid_list:
                    comid,frac = xw_.split('.')
                    if frac == '1':
                        frac = 1
                    else:
                        lfrac = len(frac)
                        frac = float(frac)/10**(lfrac-1)
                    
                    if comid in frac_sum.keys():
                        frac_sum[comid] += frac
                    else:
                        frac_sum[comid] = frac
            else:
                comid,frac = xw.split('.')
   
                if frac == '1':
                        frac = 1
                else:
                    lfrac = len(frac)
                    frac = float(frac)/10**(lfrac-1)
                if comid in frac_sum.keys():
                    frac_sum[comid] += frac
                else:
                    frac_sum[comid] = frac
    
        for k,v in frac_sum.items():
            if v >= 1.01 or v <= 0.99:
                flag = True
                print(f"{v*100}% of COMID {k} represented in refactored network")
        
        return flag
        
    
    def update_refact_ids(self,id,terminal_flag):

        '''Update refactored stream IDs with globally unique IDs'''

        terminal_flag = int(terminal_flag)
 
        if int(id) != terminal_flag:
            return int(id) + 1180001804 # max comid in NWM v2.1
        else:
            return int(id) 

    
    def replace_terminal_stream_id(self,stream_network,terminal_flag):

        '''Find terminal stream ID'''

        terminal_flag = int(terminal_flag)
        
        # find IDs only in 'to' column 
        current_terminal_ids = list(set(stream_network.to.unique()) - set(stream_network.comid.unique()))
        
        # Replace terminal IDs to terminal_flag
        stream_network.loc[stream_network.to.isin(current_terminal_ids),'to'] = terminal_flag

        return stream_network


    def sequence_missing_stream_segments(self,streams,relevant_streams,terminal_flag,extrapolate=False):
        
        missing_ids = list(set(streams.comid) - set(relevant_streams.comid))

        if len(missing_ids) > 0:

            # Find headwater node(s) (IDs only in 'ID' column)
            headwaters = list(set(streams.comid) - set(streams.to))
            
            # Create dict for faster lookups
            id_to = dict(zip(streams.comid,streams.to))

            # Find missing headwaters
            missing_hw = [m for m in headwaters if m in missing_ids]
            headwater_segs = {}
            # Extrapolate to streams outside of survey data extents
            if extrapolate == 'True':
                if len(missing_hw) > 0:
                    hw = 0
                    # For each headwater
                    for n in missing_hw:
                        hw += 1
                        # Add first downstream ID to dict
                        headwater_segs[hw] = [n]
                        # Get downstream ID
                        curr_downstream = id_to[n]
                        # Add first downstream ID to dict
                        headwater_segs[hw].append(curr_downstream)
                        # Remove headwater ID from missing
                        missing_ids.remove(n)
                        
                        # Collect all subsequent missing IDs
                        while curr_downstream in missing_ids:
                            # Set pointer to previous downstream ID
                            pointer = id_to[curr_downstream]
                            # Set downstream ID to pointer to ID
                            curr_downstream = id_to[pointer]
                            # Add curr_downstream ID to dict
                            headwater_segs[hw].append(curr_downstream)
                            # Remove curr_downstream ID from missing
                            missing_ids.remove(pointer)
                        
                        # Replace key with downstream ID with survey data
                        headwater_segs[curr_downstream] = headwater_segs.pop(hw)
            else:
                # Drop any segments outside of survey data extents
                missing_ids = [m for m in missing_ids if m not in missing_hw]
            
            # Iterate through remaining missing IDs   
            downstream_segs = {}   
            for h in headwaters:
                
                if len(missing_ids) == 0:
                    break
                
                curr_node = h
                # Walk downstream until the last segment or unitl all missing IDs are located
                while not curr_node == terminal_flag:   
                    if curr_node in missing_ids:
                        if prev_node in downstream_segs.keys():
                            # Add missing stream to upstream ID key list
                            downstream_segs[prev_node].append(curr_node)
            
                        else:
                            # Add missing stream to upstream ID key
                            downstream_segs[prev_node] = [curr_node]

                        # Update missing list
                        missing_ids.remove(curr_node)
                        
                        if len(missing_ids) == 0:
                            break
                    
                    else:
                        # Advance key to next ID
                        prev_node = curr_node
                    
                    # Set node to downstream ID
                    curr_node = id_to[curr_node]          
            
            return headwater_segs, downstream_segs


    def fill_in_missing_survey_data(self,headwater_segs, downstream_segs,streams,survey_pts,crosswalk_col, junction_to,refactor):

        # Create dict for faster lookups
        id_to = dict(zip(streams.comid,streams.to))
        if refactor == 'True':
            id_xs = dict(zip(streams.comid,streams[crosswalk_col]))
            id_junc = dict(zip(streams.comid,streams[junction_to]))

        for k,v in downstream_segs.items():
            upstream_survey = survey_pts.loc[survey_pts.comid==k]
            # Find to ID of the last missing segment in the list
            downstream_to = id_to[v[-1]]
            # Get the downstream survey using the downstream to ID
            downstream_survey = survey_pts.loc[survey_pts.comid==downstream_to]
            
            if len(downstream_segs[k]) > 1:
                # Split the list in half with the remainder going to upstream
                mid = len(downstream_segs[k]) - len(downstream_segs[k])//2
                
                # Add upstream attributes to relevant_streams for first half
                assign_upstream = downstream_segs[k][:mid]
                for missing_seg in assign_upstream:
                    upstream_survey.comid = missing_seg
                    upstream_survey.to = id_to[missing_seg]
                    if refactor == 'True':
                        upstream_survey[crosswalk_col] = id_xs[missing_seg]
                        upstream_survey[junction_to] = id_junc[missing_seg]
                    survey_pts = survey_pts.append(upstream_survey)
                
                # Add downstream attributes to relevant_streams for second half
                assign_downstream = downstream_segs[k][mid:]
                for missing_seg in assign_downstream:
                    downstream_survey.comid = missing_seg
                    downstream_survey.to = id_to[missing_seg]
                    if refactor == 'True':
                        downstream_survey[crosswalk_col] = id_xs[missing_seg]
                        downstream_survey[junction_to] = id_junc[missing_seg]
                    survey_pts = survey_pts.append(downstream_survey)
            else: 
                # Add upstream attributes to relevant_streams
                for missing_seg in v:
                    upstream_survey.comid = missing_seg
                    upstream_survey.to = id_to[missing_seg]
                    if refactor == 'True':
                        upstream_survey[crosswalk_col] = id_xs[missing_seg]
                        upstream_survey[junction_to] = id_junc[missing_seg]
                    survey_pts = survey_pts.append(upstream_survey)
        
        for k,v in headwater_segs.items():
            # Get downstream survey attributes
            survey = survey_pts.loc[survey_pts.comid==k]
            for missing_seg in v:
                survey.comid = missing_seg
                survey.to = id_to[missing_seg]
                if refactor == 'True':
                    survey[crosswalk_col] = id_xs[missing_seg]
                    survey[junction_to] = id_junc[missing_seg]
                survey_pts = survey_pts.append(survey)
        
        return survey_pts




