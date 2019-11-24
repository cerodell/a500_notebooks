#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 11:42:52 2019

@author: rodell
"""





from cr500.utils import read_cabauw
import context
import os
import time



files=context.pro_data_dir.glob('cesar*.nc')



final_dict, doppler_list, dop_dict, final_dict = read_cabauw.read(files)





