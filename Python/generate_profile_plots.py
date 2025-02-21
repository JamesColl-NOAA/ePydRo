#!/usr/bin/env python3

import argparse
from pathlib import Path
USE_PYGEOS=1
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor as Pool


def generate_profile_plots_single(args):

    output_survey_pts = args[0]
    plot_dir          = args[1]
    x_axis            = args[2]
    stream            = args[3]
    
    # Save cross sections
    plot_dir.mkdir(parents=True,exist_ok=True)

    # Filter survey points
    xs = output_survey_pts.loc[output_survey_pts.comid==int(stream)]
    
    # Get elevation units
    elev_units = xs.elev_units.unique()
    if len(elev_units) == 1:
        elev_units = elev_units.item()
    else:
        'multiple units'
    
    # Get profile distance units
    if x_axis == 'xid_d_post':
        x_units = 'meters'
    else:
        x_units = 'unitless'
    
    # Build plot
    plt.plot(xs[x_axis], xs.z, marker='o', color='black',label=stream)
    plt.ylabel(f"Elevation ({elev_units})", fontsize=14)
    plt.xlabel(f"{x_axis} ({x_units})", fontsize=14)
    plt.legend()
    plt.savefig(plot_dir / f"survey xs for comid {stream} {x_axis}.png", facecolor='w')
    plt.close() 


def generate_profile_plots_2x(args):

    output_survey_pts = args[0]
    plot_dir          = args[1]
    x_xid_d_post      = args[2]
    x_xid_d           = args[3]
    stream            = args[4]

    plot_dir.mkdir(parents=True,exist_ok=True)
    
    # Filter survey points
    xs = output_survey_pts.loc[output_survey_pts.comid==int(stream)]
    # Get elevation units
    elev_units = xs.elev_units.unique()
    if len(elev_units) == 1:
        elev_units = elev_units.item()
    else:
        'multiple units'
    
    # Build plot
    fig, ax1 = plt.subplots()
    plt.ylabel(f"Elevation ({elev_units})")

    x_xid_d_line = ax1.plot(xs[x_xid_d], xs.z, marker='o', color='blue',label='xid_d_post')
    
    ax2 = ax1.twiny()
    
    x_xid_d_post_line = ax2.plot(xs[x_xid_d_post], xs.z, marker='o', color='black',label='xid_d')

    # Generate legend
    data = x_xid_d_line + x_xid_d_post_line
    labels = [i.get_label() for i in data]
    plt.legend(data, labels)
    
    # set x axis labels and title
    ax1.set_xlabel(f"{x_xid_d} (unitless)")
    ax2.set_xlabel(f"{x_xid_d_post} (meters)")
    plt.title(f"xs for stream {stream}")

    fig.tight_layout()
    plt.savefig(plot_dir / f"survey xs for comid {stream} {x_xid_d_post} and {x_xid_d}.png", facecolor='w')
    plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate profiles for natural cross sections"
    )
    parser.add_argument(
        "-outputs_dir",
        "--outputs-dir",
        help="Path to save outputs",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-jobs",
        "--jobs",
        help="Number of cores to run",
        required=False,
        type=int,
        default=1,
    )
    parser.add_argument(
        "-plots",
        "--plots",
        help="Specify which plots to generate",
        required=False,
        type=str,
        choices=['xid_d_post', 'xid_d', 'both','all'],
        default='all',
    )

    args = vars(parser.parse_args())

    outputs_dir = Path(args["outputs_dir"])
    jobs        = args["jobs"]
    plots       = args["plots"]
    
    # Load env variables
    output_survey_pts_path = outputs_dir / 'final_survey_pts.gpkg'
    output_survey_pts = gpd.read_file(output_survey_pts_path)
    output_survey_pts = pd.DataFrame(output_survey_pts.drop(columns=['geometry']))
    stream_plots = output_survey_pts.comid.unique()
    
    x_xid_d_post = 'xid_d_post'
    x_xid_d = 'xid_d'
    
    if (plots == 'xid_d_post') or (plots == 'all'):
        # Save cross sections (xid_d_post)
        plot_dir = outputs_dir / 'plots_xid_d_post'
        batch_args = [output_survey_pts,plot_dir,x_xid_d_post]
        # Generate profile plots
        with Pool(jobs) as p:
            pool_list = p.map(generate_profile_plots_single, [batch_args + [stream] for stream in stream_plots])
        
    if (plots == 'xid_d') or (plots == 'all'):
        # Save cross sections (xid_d)
        plot_dir = outputs_dir / 'plots_xid_d'
        batch_args = [output_survey_pts,plot_dir,x_xid_d]
        # Generate profile plots
        with Pool(jobs) as p:
            pool_list = p.map(generate_profile_plots_single, [batch_args + [stream] for stream in stream_plots])
    
    if (plots == 'both') or  (plots == 'all'):
        # Save cross sections (both xid_d and xid_d_post)
        plot_dir = outputs_dir / 'plots_both'
        batch_args = [output_survey_pts,plot_dir,x_xid_d_post,x_xid_d]
        # Generate profile plots
        with Pool(jobs) as p:
            pool_list = p.map(generate_profile_plots_2x, [batch_args + [stream] for stream in stream_plots])
    