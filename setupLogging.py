'''
Created on Oct 30, 2018

@author: Andrew Davidson aedavids@ucsc.edu
'''

import os
import json
import logging.config

_alreadyInitialized = False

def setupLogging (
        default_path='logging.ini.json',  #'logging.json',
        default_level=logging.WARN,
        env_key='LOG_CFG'):
    """
    configure python logging facilities. Be default uses config file logging.ini.json
    """
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        logger = logging.getLogger('root')
        logger.warn("loading log configuration:{}".format(default_path))
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)

