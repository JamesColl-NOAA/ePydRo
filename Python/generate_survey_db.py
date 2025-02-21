#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
from preprocess import Preprocess

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Populate survey_db"
    )
    parser.add_argument(
        "-working_dir",
        "--working-dir",
        help="Path of working directory",
        required=True,
        type=str,
    )

    parser.add_argument(
        "-models_dir",
        "--models-dir",
        help="Path to HEC-RAS models",
        required=True,
        type=str,
    )

    args = vars(parser.parse_args())

    working_dir = Path(args["working_dir"])
    models_dir = Path(args["models_dir"])
    survey_db_pre_path = working_dir / "models_ready_for_proj_data.csv"

    # Ignore models with known missing projections
    models_missing_proj = ['San Jacinto - HCFCD - TSARP',
                    'San Jacinto - Buffalo Bayou Houston Ship Channel',
                    'Rio Grande - Unknown',
                    'Rio Grande - Columbia - AHPS FIM',
                    'Rio Grande - Halff Associates - AHPS FIM',
                    'Backup']

    # Verify that all HDF files exist
    hdf_files, status = Preprocess().check_for_missing_hdf_files(working_dir,models_dir,models_missing_proj)
    print (status)

    # Create initial survey_db
    Preprocess().create_initial_survey_db(hdf_files,models_dir,survey_db_pre_path)

    ''' Manual step of adding prj files to initial survey_db '''

    survey_db_clean_filename = working_dir / "survey_db.csv"
    survey_db_filename = working_dir / "models_with_proj_data_edited.csv"
    survey_db = pd.read_csv(survey_db_filename)
    hdf_folder = working_dir / 'hdf_data'
    hdf_folder.mkdir(parents=True,exist_ok=True)

    # Create survey_db and copy hdf files to single directory
    survey_db_clean = Preprocess().create_survey_db(survey_db,hdf_folder)
    survey_db_clean.to_csv(survey_db_clean_filename,index=False)