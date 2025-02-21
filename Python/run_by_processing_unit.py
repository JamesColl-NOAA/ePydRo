#!/usr/bin/env python3

import argparse
from os import getenv
from pathlib import Path

import geopandas as gpd
from database import Spatial
from dotenv import load_dotenv


def run_by_unit(args):

    survey_lines = args[0]
    set_proj     = args[1]
    id_col       = args[2]
    distance     = args[3]
    survey_table = args[4]
    outputs_dir  = args[5]
    group        = args[6]

    load_dotenv(dotenv_path="/repo/python/config/.env")
    output_streams_filename = getenv("output_streams_filename")
    aoi_filename = getenv("aoi_filename")
    save_output = getenv("save_output")
    testing = getenv("testing")
    refactor = getenv("refactor")
    crosswalk_col = getenv("crosswalk_col")
    junction_to = getenv("junction_to")
    extrapolate = getenv("extrapolate")
    terminal_flag = int(getenv("terminal_flag"))
    
    if testing == "True":  # run LRCA test bed
        if refactor == "True":  # run refactored LRCA test bed # TODO: Move this to it's own workflow
            
            print("Running refactored stream network for LCRA")
            # Read in streams subset
            refac_test_stream_path = Path(getenv("refac_test_stream_path"))
            streams = gpd.read_file(refac_test_stream_path)
            
            # Give refactored stream globablly unique IDs
            streams['old_to'] = streams.to
            streams['old_comid'] = streams.comid
            streams.to = streams.apply(lambda x: Spatial().update_refact_ids(x['to'],terminal_flag),axis=1)
            streams.comid = streams.apply(lambda x: Spatial().update_refact_ids(x['comid'],terminal_flag),axis=1)
           
        else:
            test_stream_path = Path(getenv("test_stream_path"))
            
            if test_stream_path.exists():
                # Read in streams subset
                streams = gpd.read_file(test_stream_path)
            else:
                # Read in full streams and subset
                nwm_v2_1_stream_path = Path(getenv("nwm_v2_1_stream_path"))
                
                streams_all = gpd.read_file(nwm_v2_1_stream_path)
                test_id_list = getenv("test_id_list")

                with open(test_id_list) as f:
                    comid_list = f.read().splitlines()

                comid_list = [int(i) for i in comid_list]
                streams = streams_all.loc[streams_all.ID.isin(comid_list)]
                Spatial().write_spatial(streams, test_stream_path)
                del streams_all

    else:

        # Subset surveys by group
        survey_lines = survey_lines.loc[survey_lines.rfc == group]

        if survey_lines.empty:
            return

        # Create convex_hull for group of survey lines
        proc_unit = Spatial().create_convex_hull(survey_lines)

        # Optional save aoi to spatial layer
        if save_output:
            aoi = gpd.GeoDataFrame(
                {"name": group, "geometry": proc_unit},
                geometry="geometry",
                crs=survey_lines.crs,
                index=[0],
            )
            print(f"Saving {group} aoi polygon to spatial layer")
            basename, ext = aoi_filename.split(".")
            aoi_filename = f"{basename}_{group}.{ext}"
            aoi_path = outputs_dir / aoi_filename
            Spatial().write_spatial(aoi, aoi_path)
            del aoi

        # Read in masked streams
        nwm_v2_1_stream_path = Path(getenv("nwm_v2_1_stream_path"))
        streams = gpd.read_file(nwm_v2_1_stream_path, mask=proc_unit)

    if refactor == "False":
        streams = streams[["ID", "to", "order_", "Length", "Lake", "gages", "geometry"]]
        streams = streams.rename(columns={"ID": "comid","order_": "order","Lake": "NHDWaterbodyComID"})
        streams = Spatial().replace_terminal_stream_id(streams,terminal_flag)
    
    streams.gages = streams.gages.str.strip()
    streams = streams.to_crs(survey_lines.crs)
    streams = streams.explode()
    streams = streams.reset_index(drop=True)

    # Subset streams with intersecting surveys
    print(f"Subset streams with surveys for {group}")
    relevant_streams = Spatial().subset_streams_with_survey_intersections(
        streams, survey_lines, refactor,
        crosswalk_col, junction_to
    )

    # Optional save streams to spatial layer
    if save_output:
        print("Saving survey points to spatial layer")
        relevant_streams_save = Spatial().check_data_types(relevant_streams.copy())
        basename, ext = output_streams_filename.split(".")
        output_streams_filename = f"{basename}_{group.replace(' ','_')}.{ext}"
        streams_out_path = outputs_dir / output_streams_filename
        Spatial().write_spatial(relevant_streams_save, streams_out_path)
        del relevant_streams_save
    
    # Select cross section points using stream segment distance
    print(f"Collecting surveys closest to {distance.astype(float)*100}% stream segment length for {group}")
    selected_surveys = Spatial().get_xs_closest_to_stream_point(
        relevant_streams, survey_lines, id_col, distance, refactor,
        crosswalk_col, junction_to
    )
    
    # Select cross section points using stream segment distance
    print(f"Select cross section points using stream segment distance for {group}")
    survey_pts = Spatial().calculate_relative_survey_pt_dist(
        selected_surveys, id_col, survey_table, set_proj
    )
    
    # Check for serially complete network
    print("Checking for serially complete network")
    headwater_segs, downstream_segs = Spatial().sequence_missing_stream_segments(streams,relevant_streams,terminal_flag,extrapolate)

    # Enforce serially complete network
    print("Enforce serially complete network")
    survey_pts = Spatial().fill_in_missing_survey_data(headwater_segs, downstream_segs,streams,survey_pts,crosswalk_col, junction_to,refactor)
    
    if testing == "True":  # run LRCA test bed
        # TODO: adjacent streams need to be removed from relevant streams before running this on everything; NAs produced resulting in errors
        # Enforce data types
        print(f"Enforce data types for {group}")
        survey_pts = Spatial().check_data_types(survey_pts)
    
    return survey_pts


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Clean up survey data and prep for routelink netcdf"
    )
    parser.add_argument(
        "-survey_lines", "--survey-lines", help="Survey lines layer", required=True
    )
    parser.add_argument(
        "-id_col", "--id-col", help="Id column name", required=True, type=str
    )
    parser.add_argument(
        "-distance", "--distance", help="Distance along line", required=True, type=float
    )
    parser.add_argument(
        "-set_proj", "--set-proj", help="Define data projecton", required=True, type=str
    )
    parser.add_argument(
        "-survey_table", "--survey-table", help="Survey table dataframe", required=True
    )

    args = vars(parser.parse_args())

    run_by_unit(**args)
