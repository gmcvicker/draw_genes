
import sys

import numpy as np
import rpy2.robjects as robjects

import genome.db

from continuoustrack import ContinuousTrack
from basellrtrack import BaseLLRTrack


class LLRTrack(BaseLLRTrack):
    """LLRs are a type of continuous track that are symmetric about 0.
    Draws negative values below an axis at 0.0 (rather than setting the
    axis at the minimum value for the track and drawing all values above
    it)"""
    def __init__(self, region, options):
        source = "gdb"

        if "source" in options:
            source = options['source']

        if source == "gdb":
            track_name = options['track']
            gdb = options['gdb']
            track = gdb.open_track(track_name)
            values = track.get_nparray(region.chrom, start=region.start,
                                       end=region.end)
            track.close()

        if source == "wig":
            path = options['path']
            values = genome.wig.read_ints(path, region)


        if "scale" in options:
            scale = float(options['scale'])
            # sys.stderr.write("scaling values by %.2f\n" % scale)
            values = values * scale

        super_init = super(LLRTrack, self).__init__
        super_init(values, region, options)
        


    

                


