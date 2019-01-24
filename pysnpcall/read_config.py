import numpy as np

class read_config(object):

    def __init__(self, config_file):

        # Set class vaiables
        self.config_file = config_file

    def get_config_contents(self):

        config_dict = {}
        config_file_handle = open(self.config_file)
        config_values = config_file_handle.readlines()

        for config_value in config_values:
            if not config_value.startswith('#'):
                if len(config_value) > 1:
                    components = config_value.split('=')
                    config_dict[ components[0].strip() ] =  components[1].strip()

        return config_dict
