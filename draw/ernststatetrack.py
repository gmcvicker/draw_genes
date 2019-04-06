
import sys

from .statetrack import StateTrack

import rpy2.robjects as robjects


class ErnstStateTrack(StateTrack):
    """This is a convenience class which inherits from the StateTrack class
    in order to provide a bunch of default configuration options for
    Ernst Chromatin State tracks. This removes the burden of having to
    specify a large number of labels and colors for every such track in
    the configuration file."""
    
    def __init__(self, region, options):

        # set a bunch of default configuration options if they are not already
        # specified
        defaults = {'n_state' : 15,
                    'state_label_0' : 'none',
                    'state_color_0' : "white",

                    'state_label_1' : '1 strong promoter',
                    'state_color_1' : "#e41a1c",
                    
                    'state_label_2' : '2 weak promoter',
                    'state_color_2' : "#fbb4ae",

                    "state_label_3" : "3 poised promoter",
                    "state_color_3" : "#377eb8",

                    "state_label_4" : "4 strong enhancer",
                    "state_color_4" : "#4daf4a",

                    "state_label_5" : "5 strong enhancer",
                    "state_color_5" : "#984ea3",

                    "state_label_6" : "6 weak enhancer",
                    "state_color_6" : "#ccebc5",

                    "state_label_7" : "7 weak enhancer",
                    "state_color_7" : "#decbe4",

                    "state_label_8" : "8 insulator",
                    "state_color_8" : "#ff7f00",

                    "state_label_9" : "9 txn transition",
                    "state_color_9" : "#ffff33",

                    "state_label_10" : "10 txn elongation",
                    "state_color_10" : "#f781bf",

                    "state_label_11" : "11 repressed",
                    "state_color_11" : "#a65628",

                    "state_label_12" : "12 weak txn",
                    "state_color_12" : "#fddaec",

                    "state_label_13" : "13 heterochrom/lo",
                    "state_color_13" : "grey40",

                    "state_label_14" : "14 repetitive/cnv",
                    "state_color_14" : "grey60",

                    "state_label_15" : "15 repetitive/cnv",
                    "state_color_15" : "grey80"}

        for key, val in list(defaults.items()):
            if key not in options:
                options[key] = val
                        
        super(ErnstStateTrack, self).__init__(region, options)
        

