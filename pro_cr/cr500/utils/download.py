#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:18:05 2019

@author: rodell
"""

"""
@author: rodell
"""

import argparse
import requests
from pathlib import Path
import context
import shutil
import cr500
import pandas as pd
import pdb

class NoDataException(Exception):
    pass

def download(
    filename,
    root="http://breeze.colorado.edu/ftp/XPIA/WC16/",
    dest_folder=None,
):
    """
    copy file filename from http://breeze.colorado.edu/ftp/XPIA/WC16/ to 
    the local directory.  If local file exists, report file size and quit.
    Parameters
    ----------
    filename: string
      name of file to fetch from 
    root: optional string 
          to specifiy a different download url
    dest_folder: optional string or Path object
          to specifify a folder besides the current folder to put the files
          will be created it it doesn't exist
    Returns
    -------
    Side effect: Creates a copy of that file in the local directory
    """
    url = "{}/{}".format(root, filename)
    url = url.replace("\\", "/")
    print("trying {}".format(url))
    #
    # use current directory if dest_dir not specified
    #
    if dest_folder is None:
        dest_path = Path()
    else:
        dest_path = Path(dest_folder).resolve()
        dest_path.mkdir(parents=True, exist_ok=True)
    #
    # filename may contain subfolders
    #
    filepath = Path(filename)
    filename = filepath.name
    filepath = dest_path / Path(filename)
    print(f"writing to: {filepath}")
    if filepath.exists():
        the_size = filepath.stat().st_size
        print(
            ("\n{} already exists\n" "and is {} bytes\n" "will not overwrite\n").format(
                filename, the_size
            )
        )
        return None

    tempfile = str(filepath) + "_tmp"
    temppath = Path(tempfile)
    try:
        with open(temppath, "wb") as localfile:
            print(f"writing temporary file {temppath}")
            response = requests.get(url, stream=True)
            #
            # treat a 'Not Found' response differently, since you want to catch
            # this and possibly continue with a new file
            #
            if not response.ok:
                if response.reason == "Not Found":
                    the_msg = 'requests.get() returned "Not found" with filename {}'.format(
                        filename
                    )
                    raise NoDataException(the_msg)
                else:
                    #
                    # if we get some other response, raise a general exception
                    #
                    the_msg = "requests.get() returned {} with filename {}".format(
                        response.reason, filename
                    )
                    raise RuntimeError(the_msg)
                    #
                # clean up the temporary file
                #
            for block in response.iter_content(1024):
                if not block:
                    break
                localfile.write(block)
        the_size = temppath.stat().st_size
        print("downloaded {}\nsize = {}".format(filename, the_size))
        shutil.move(str(temppath), str(filepath))
        if the_size < 10.0e3:
            print(
                "Warning -- your file is tiny (smaller than 10 Kbyte)\nDid something go wrong?"
            )
    except NoDataException as e:
        print(e)
        print("clean up: removing {}".format(temppath))
        temppath.unlink()
    return None


