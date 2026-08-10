#!/usr/bin/env python
""" This script takes as input a text file with the software versions reported in format:
- name: software_name1
  version: version1 
  homepage: software_homepage1
  workflow task: workflow_task_name1
- name: software_name2
  version: version2
  homepage: software_homepage2
  workflow task: workflow_task_name2

and outputs a markdown table with the format:

| Program       | Version  | Relevant Links     |
|:--------------|:---------|:-------------------|
| sofware_name1 | version1 | software_homepage1 |
| sofware_name2 | version2 | software_homepage2 |

"""

from pathlib import Path

import yaml
import click
import pandas as pd


AGILENT_SOFTWARE_DPPD = [
    "R",
    "Bioconductor",
    "DT",
    "dplyr",
    "stringr",
    "limma",
    "glue",
    "ggplot2",
    "matrixStats",
    "statmod",
    "dp_tools",
    "singularity",
    "Quarto",
    "nextflow"
]

AGILENT_SOFTWARE_DPPD = [s.lower() for s in AGILENT_SOFTWARE_DPPD]

ASSUMED_SOFTWARE = [{
    "name": "singularity",
    "version": 3.9,
    "homepage": "https://sylabs.io"
}]


## Used when the R library metadata doesn't encode any URLS
HOMEPAGE_MAP = {
    "statmod":"https://cran.r-project.org/web/packages/statmod/index.html",
}


def yaml_to_markdown(input_yaml: Path, filename: str, skip_de: bool):
    """ Using a software versions """
    with open(input_yaml, "r") as f:
        data = yaml.safe_load(f)

    data.extend(ASSUMED_SOFTWARE)
    df = pd.DataFrame(data)

    if skip_de:
        AGILENT_SOFTWARE_DPPD.remove('matrixstats')
        AGILENT_SOFTWARE_DPPD.remove('statmod')

    # Filter to direct software used (i.e. exclude dependencies of the software)
    df = df.loc[df["name"].str.lower().isin(AGILENT_SOFTWARE_DPPD)]

    assert len(AGILENT_SOFTWARE_DPPD) == len(df), f"Not all software accounted for! Missing: {set(AGILENT_SOFTWARE_DPPD) - set(df['name'].str.lower())}"

    df['homepage'] = df.apply(lambda row: HOMEPAGE_MAP[row['name']]
                              if row['homepage'] == "NO URLS ENCODED"
                              else row['homepage'], axis="columns")

    df = df[["name", "version", "homepage"]]
    print(df)
    df = df.rename({"name":"Program","version":"Version","homepage":"Relevant Links"},
                   axis="columns")

    # Sort by program name for deterministic output
    df = df.sort_values("Program")
    df.to_markdown("software_versions_GLmicroarray.md", index=False)

if __name__ == "__main__":
    @click.command()
    @click.argument("input_yaml", type=click.Path(exists=True))
    @click.argument("filename")
    @click.argument("skip_de", type=click.BOOL)
    def cli(input_yaml, filename, skip_de):
        yaml_to_markdown(Path(input_yaml), filename, skip_de)

    cli()
