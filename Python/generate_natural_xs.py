#!/usr/bin/env python3

import argparse
from multiprocessing import Pool
from os import getenv
from os.path import splitext
import geopandas as gpd
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from database import Spatial
from dotenv import load_dotenv
from run_by_processing_unit import run_by_unit


def generate_bathy(inputs_dir, outputs_dir):

    # Load env variables
    load_dotenv(dotenv_path="/repo/python/config/.env")
    jobs = int(getenv("jobs"))
    data_type = getenv("data_type") # TODO: replace with survey attribute
    proj_for_distance = getenv("proj_for_distance")
    final_survey_pts_filename = getenv("final_survey_pts_filename")
    final_survey_lines_filename = getenv("final_survey_lines_filename")
    netcdf_filename = getenv("netcdf_filename")
    set_proj = getenv("set_proj")
    id_col = getenv("id_col")
    distance_on_line = np.array(getenv("distance_on_line"))
    save_output = getenv("save_output")
    testing = getenv("testing")
    refactor = getenv("refactor")
    terminal_flag = int(getenv("terminal_flag"))
    inputs_dir = Path(inputs_dir)
    outputs_dir = Path(outputs_dir)
    survey_lines_filename = getenv("survey_lines_filename")
    merged1d_filename = getenv("merged1d_filename")

    # File names
    final_survey_pts_path = outputs_dir / final_survey_pts_filename
    final_survey_lines_path = outputs_dir / final_survey_lines_filename
    netcdf_path = outputs_dir / netcdf_filename

    # Read in survey table
    survey_table_path = inputs_dir / merged1d_filename
    survey_table = pd.read_csv(survey_table_path, sep=",")

    # Reading survey lines
    print("Reading survey lines")
    survey_lines_path = inputs_dir / survey_lines_filename
    if splitext(survey_lines_filename)[1] == '.feather':
        survey_lines = gpd.read_feather(survey_lines_path)
    else:
        survey_lines = gpd.read_file(survey_lines_path)

    # Reproject surveys to get units in meters
    survey_lines = survey_lines.to_crs(proj_for_distance)
    
    # Read in AOI for masking stream data
    if testing == "True":  # run LRCA test bed

        survey_pts = run_by_unit(
            [
                survey_lines,
                set_proj,
                id_col,
                distance_on_line,
                survey_table,
                outputs_dir,
                "LRCA test",
            ]
        )

    else:
        rfc_group = survey_lines.rfc.unique()

        # Run groups (currently by RFC) in parallel
        batch_args = [
            survey_lines,
            set_proj,
            id_col,
            distance_on_line,
            survey_table,
            outputs_dir,
        ]
        with Pool(jobs) as p:
            group_list = p.map(run_by_unit, [batch_args + [group] for group in rfc_group])

        group_list = [group for group in group_list if not group.empty]
        survey_pts = pd.concat(group_list)

    # Optional Save survey points to spatial layer
    if save_output:
        print("Saving survey points to spatial layer")
        Spatial().write_spatial(survey_pts, final_survey_pts_path)

    # Optional Save selected survey lines to spatial layer
    if save_output:
        print("Saving survey lines to spatial layer")
        select_ids = survey_pts.xid.unique()
        survey_lines = survey_lines.loc[survey_lines.xid.isin(select_ids)]
        Spatial().write_spatial(survey_lines, final_survey_lines_path)

    if refactor == "True":
        print("Saving refactored ID yaml file")
        refactor_id_yaml_filename = getenv("refactor_id_yaml_filename")
        refactor_id_yaml_path = outputs_dir / refactor_id_yaml_filename
        
        # dict where key = terminal stream and values = all IDs in domain (including key)
        refac_ids = np.unique(survey_pts.comid.to_list())
        # TODO: update to handle multiple terminal streams
        yaml_key = survey_pts.loc[survey_pts.to==terminal_flag].comid.unique().item()
        yaml_dict = {yaml_key:refac_ids.tolist()}
        
        # Save refactored IDs dict to yaml file
        with open(refactor_id_yaml_path, 'w') as file:
            yaml.dump(yaml_dict, file,explicit_start=False, default_flow_style=False)

    print("Export surveys to routelink netcdf file")
    if refactor == "True":
        crosswalk_col = getenv("crosswalk_col")
        junction_to = getenv("junction_to")
        survey_pts = survey_pts.filter(items=['comid', 'to', 'order', 'Length', 'n', 'NHDWaterbodyComID', 'gages',
        crosswalk_col, junction_to, 'xid_d', 'z'])
        survey_pts = survey_pts.rename(columns={'comid': 'rlink','to': 'rto'})
    else:
        survey_pts = survey_pts.filter(items=['comid', 'to', 'order', 'Length','n', 'NHDWaterbodyComID', 'gages',
        'xid_d', 'z'])
        survey_pts = survey_pts.rename(columns={'comid': 'link'})
    
    if testing == "True":
        version = "LCRA"
    else:
        version = "Full Diffusive Domain"
    
    Spatial().export_surveys(survey_pts, data_type, set_proj, version, netcdf_path)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Process surveys and export as routelink netcdf"
    )
    parser.add_argument(
        "-inputs_dir",
        "--inputs-dir",
        help="Path to inputs",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-outputs_dir",
        "--outputs-dir",
        help="Path to save outputs",
        required=True,
        type=str,
    )

    args = vars(parser.parse_args())

    generate_bathy(**args)
