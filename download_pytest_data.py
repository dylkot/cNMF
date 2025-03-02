#!/usr/bin/env python3

import os
import sys
import urllib.request
import tarfile
import subprocess

def download_file(url, dest_path):
    """
    Download a file from `url` to `dest_path`.
    """
    print(f"Downloading {url} -> {dest_path}")
    urllib.request.urlretrieve(url, dest_path)
    print("Download complete.")

def extract_file(archive_path, extract_dir):
    """
    Extract .zip or .tar(.gz) archives to `extract_dir`.
    Add more conditions if you handle other formats.
    """
    print(f"Extracting {archive_path} into {extract_dir}")
    if archive_path.endswith(".zip"):
        with zipfile.ZipFile(archive_path, 'r') as z:
            z.extractall(extract_dir)
    elif archive_path.endswith(".tar.gz") or archive_path.endswith(".tgz"):
        with tarfile.open(archive_path, 'r:gz') as t:
            t.extractall(extract_dir)
    elif archive_path.endswith(".tar"):
        with tarfile.open(archive_path, 'r:') as t:
            t.extractall(extract_dir)
    else:
        print(f"Skipping extraction for {archive_path} (unrecognized archive format).")

def main():
    # List of example data to download
    # Each entry could have "url", "filename", and possibly "extract" if it's an archive
    data_files = [
        {
            "url": "https://storage.googleapis.com/sabeti-public/dkotliar/cNMF/pytest_testing_data/example_pbmc_results_20250301.tar.gz",
            "filename": "example_pbmc_results_20250301.tar.gz",
            "extract": True,
        },
        {
            "url": "https://storage.googleapis.com/sabeti-public/dkotliar/cNMF/pytest_testing_data/example_sim_results_20250301.tar.gz",
            "filename": "example_sim_results_20250301.tar.gz",
            "extract": True
        }
    ]

    # Where to place the downloaded files
    data_dir = os.path.join("tests", "test_data")
    os.makedirs(data_dir, exist_ok=True)

    # 1) Download each file and (optionally) extract it
    for info in data_files:
        file_url = info["url"]
        file_name = info["filename"]
        local_path = os.path.join(data_dir, file_name)

        # Download
        download_file(file_url, local_path)

        # Extract if indicated
        if info.get("extract", False):
            extract_file(local_path, data_dir)

        # Remove tar.gz
        os.remove(os.path.join(data_dir, info['filename']))

if __name__ == "__main__":
    main()