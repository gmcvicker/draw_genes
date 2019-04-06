import sys

import numpy as np

from .continuoustrack import ContinuousTrack


class GCContentTrack(ContinuousTrack):
     def __init__(self, region, options):
          gdb = options['gdb']
          
          if 'track' in options:
               track_name = options['track']
          else:
               track_name = "seq"
          
          track = gdb.open_track(track_name)

          # retrieve data from sequence track
          seq_ascii_vals = track.get_nparray(region.chrom, start=region.start,
                                             end=region.end)
          track.close()

          # convert to 0s and 1s, with Gs and Cs as 1s
          values = (seq_ascii_vals == ord("G")) | (seq_ascii_vals == ord("C"))

          # plot values as continuous track
          super_init = super(GCContentTrack, self).__init__
          super_init(values, region, options)

