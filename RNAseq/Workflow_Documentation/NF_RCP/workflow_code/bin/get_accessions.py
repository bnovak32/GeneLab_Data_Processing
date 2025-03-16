#!/usr/bin/env python3

import requests
import argparse
import re
import sys
import json

def get_osd_and_glds(accession, api_url):
    # Fetch data from the API
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from API: {e}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print("Error decoding JSON response from API", file=sys.stderr)
        sys.exit(1)

    osd_accession = None
    glds_accessions = []

    # Check if the accession is OSD or GLDS
    if accession.startswith('OSD-'):
        osd_accession = accession
        # Find the GLDS identifiers associated with this OSD accession
        for hit in data.get("hits", {}).get("hits", []):
            source = hit.get("_source", {})
            if source.get("Accession") == osd_accession:
                identifiers = source.get("Identifiers", "")
                glds_accessions = re.findall(r'GLDS-\d+', identifiers)
                break
    elif accession.startswith('GLDS-'):
        glds_accessions = [accession]
        # Find the OSD accession associated with this GLDS accession
        for hit in data.get("hits", {}).get("hits", []):
            source = hit.get("_source", {})
            identifiers = source.get("Identifiers", "")
            if accession in identifiers:
                osd_accession = source.get("Accession")
                break
    else:
        print("Invalid accession format. Please use 'OSD-###' or 'GLDS-###'.", file=sys.stderr)
        sys.exit(1)

    if not osd_accession or not glds_accessions:
        print(f"No data found for {accession}", file=sys.stderr)
        sys.exit(1)

    return osd_accession, glds_accessions

def main():
    parser = argparse.ArgumentParser(description="Retrieve OSD and GLDS accessions.")
    parser.add_argument('--accession', required=True, help="Accession in the format 'OSD-###' or 'GLDS-###'")
    parser.add_argument('--api_url', required=True, help="OSDR API URL")
    args = parser.parse_args()

    osd_accession, glds_accessions = get_osd_and_glds(args.accession, args.api_url)

    # Output the results in a way that Nextflow can capture
    print(f"{osd_accession}")
    print(f"{','.join(glds_accessions)}")

if __name__ == "__main__":
    main()
