#!/usr/bin/env python3

import argparse
from concurrent.futures import ProcessPoolExecutor as Pool
from os import getenv, listdir
from pathlib import Path

import geopandas as gpd
USE_PYGEOS=1
import pandas as pd
from database import Spatial
from preprocess import Preprocess
from dotenv import load_dotenv


def preprocess_survey_data(survey_path):

    # Load env variables
    load_dotenv(dotenv_path="/repo/python/config/.env")

    # Input Paths
    hdf_data_dir = Path(getenv("hdf_data_dir"))
    survey_lines_path = Path(getenv("survey_lines_filename"))

    # Outputs
    merged1d_path = Path(getenv("merged1d_filename"))
    survey_db_filename = getenv("survey_db_filename")
    refac_test_stream_path = Path(getenv("refac_test_stream_path"))

    # Create directories
    outputs_dir.mkdir(parents=True, exist_ok=True)

    # Parameters
    set_proj = getenv("set_proj")
    jobs = int(getenv("jobs"))
    id_col = getenv("id_col")
    refactor = getenv("refactor")
    overwrite_preprocessing = getenv("overwrite_preprocessing")

    # Read model database
    survey_db = pd.read_csv(survey_path / survey_db_filename)
    # Check that all surveys are in hdf folder
    survey_hdf_mismatch = list(set(survey_db.hdf_path) - set(listdir(hdf_data_dir)))
    if len(survey_hdf_mismatch) > 0:
        print(
            f"mismatch between hdf directory and survey db: missing {len(survey_hdf_mismatch)} files"
        )

    # Check for missing survey data
    if not merged1d_path.exists() or overwrite_preprocessing == 'True':
        
        # Generate args for each model
        arg_list = []
        for _, model in survey_db.iterrows():

            hdf_file = hdf_data_dir / str(model.hdf_path)
            units = model.model_unit
            
            if model.proj_type == "string":
                from_proj = model.proj
            else:
                print(f"model {model.hdf_path} has no projection in survey db")
                continue
            
            arg_list.append([hdf_file,units,from_proj,set_proj])

        # Extract survey data for each model
        with Pool(jobs) as p:
            pool_list = p.map(Preprocess().generate_py_xyz_from_hdf, [args for args in arg_list])
        
        # Remove empty survey datasets
        pool_list = [df for df in pool_list if not df.empty]
        
        # Combine survey datasets
        survey_data = pd.concat(pool_list, ignore_index=True)

        # Validate survey data conversion
        missing_surveys = list(set(listdir(hdf_data_dir)) - set(survey_data.model_name.unique()))
        if len(missing_surveys) > 0:
            print(f"{len(missing_surveys)} surveys failed preprocessing")

        del arg_list, pool_list, missing_surveys

        # Reindex survey ids
        survey_data = Preprocess().reindex_survey_ids(survey_data)
        
        # Remove duplicate survey points
        survey_data = survey_data[
            ~survey_data.duplicated(
                subset=["xid_d", "x", "y", "z", "n"], keep="first"
            )
        ]

        # Add attributes to survey table
        survey_db_fields = survey_db[["hdf_path", "rfc"]]
        survey_data = pd.merge(
            survey_data,
            survey_db_fields,
            how="left",
            left_on="model_name",
            right_on="hdf_path",
        )

        # Find surveys with single point
        count_table = survey_data.groupby(["xid"]).size().reset_index(name="counts")
        single_pt_surveys = list(count_table.loc[count_table.counts < 2].xid.unique())

        # Drop surveys with single point
        survey_data = survey_data.loc[~survey_data.xid.isin(single_pt_surveys)]

        # Save merged/filtered survey data
        print("Writing out aggregate survey data and source referece table")
        survey_data.to_csv(merged1d_path, index=False)

    # Convert survey data to lines
    if not survey_lines_path.exists() or overwrite_preprocessing == 'True':

        print("Convert survey data to lines")
        survey_lines = Spatial().survey_table_to_lines(survey_data, set_proj, id_col)

        # Add survey attributes to survey lines
        survey_data_fields = survey_data[["xid", "source", "model_name", "rfc"]]
        survey_data_fields = survey_data_fields.drop_duplicates()
        survey_lines = pd.merge(survey_lines, survey_data_fields, how="left", on="xid")

        # Save survey lines to spatial layer
        print("Saving survey lines to spatial layer")
        Spatial().write_spatial(survey_lines, survey_lines_path)
    
    # Subset refactored network
    if refactor == "True":
        if not refac_test_stream_path.exists() or overwrite_preprocessing == 'True':
            
            crosswalk_col = getenv("crosswalk_col")
            rpu_refac_test_stream_path = getenv("rpu_refac_test_stream_path")
            testing = getenv("testing")
            terminal_flag = int(getenv("terminal_flag"))
            junction_to = getenv("junction_to")
            
            if testing == "True":
                test_stream_path = Path(getenv("test_stream_path"))
                lcra_streams = gpd.read_file(test_stream_path)
                stream_ids = lcra_streams.ID.unique()
            
            else: # TODO: This layer doesn't exist yet
                full_domain_stream_path = Path(getenv("full_domain_stream_path"))
                full_domain_streams = gpd.read_file(full_domain_stream_path)
                stream_ids = full_domain_streams.ID.unique()

            # Subset refactored network and to column for incoming tributaries
            refactored_network = Spatial().subset_refactored_network(stream_ids, rpu_refac_test_stream_path,crosswalk_col,lcra_streams,junction_to)

            # Flag terminal stream
            refactored_network = Spatial().replace_terminal_stream_id(refactored_network,terminal_flag)
            
            # Check network continuity
            flag = Spatial().check_crosswalk_fractional_sum(refactored_network,crosswalk_col)
            if flag == "True":
                print(f"Refactored network not complete. Check crosswalk to ensure continuity is preserved betwen networks.")
            else:
                Spatial().write_spatial(refactored_network,refac_test_stream_path)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Preprocess survey data, assign to relevant streams, and export as routelink netcdf"
    )
    parser.add_argument(
        "-outputs_dir",
        "--outputs-dir",
        help="Path to save outputs",
        required=True,
        type=str,
    )

    args = vars(parser.parse_args())

    outputs_dir = Path(args["outputs_dir"])

    # Load env variables
    load_dotenv(dotenv_path="/repo/python/config/.env")

    # Input Paths
    survey_path = Path(getenv("survey_path"))

    # Preprocess survey data
    print('Preprocessing survey data')
    preprocess_survey_data(survey_path)