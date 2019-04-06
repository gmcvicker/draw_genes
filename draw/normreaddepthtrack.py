import sys

import numpy as np

import genome.trackstat

from .continuoustrack import ContinuousTrack
from .basellrtrack import BaseLLRTrack

class NormReadDepthTrack(BaseLLRTrack):
     def __init__(self, region, options):
          gdb = options['gdb']

          track_name1 = options['track1']
          track_name2 = options['track2']

          pseudo_count = float(options['pseudocount'])
          
          track1 = gdb.open_track(track_name1)
          track2 = gdb.open_track(track_name2)
          
          values1 = track1.get_nparray(region.chrom, start=region.start,
                                      end=region.end)
          
          values2 = track2.get_nparray(region.chrom, start=region.start,
                                      end=region.end)

          if "scale_factor1" in options:
               scale_factor1 = float(options['scale_factor1'])
               try:
                    stat = gdb.get_track_stat(track1)
                    scale1 = scale_factor1 / float(stat.sum)
                    values1 = values1 * scale1
                    sys.stderr.write("  total reads %d, using "
                                     "scale %.3f\n" % (stat.sum, scale1))
               except ValueError as err:
                    sys.stderr.write("  WARNING: cannot scale values for "
                                     "track %s because track stats have "
                                     "not been calculated. run "
                                     "set_track_stats.py first.\n" %
                                     track.name)
          track1.close()          

          if "scale_factor2" in options:
               scale_factor2 = float(options['scale_factor2'])
               try:
                    stat = gdb.get_track_stat(track2)
                    scale2 = scale_factor2 / float(stat.sum)
                    values2 = values2 * scale2
                    sys.stderr.write("  total reads %d, using "
                                     "scale %.3f\n" % (stat.sum, scale2))
               except ValueError as err:
                    sys.stderr.write("  WARNING: cannot scale values for "
                                     "track %s because track stats have "
                                     "not been calculated. run "
                                     "set_track_stats.py first.\n" %
                                     track.name)
          
          track2.close()

          ratio = np.log2((values1 + pseudo_count) / (values2 + pseudo_count))
          
          sys.stderr.write("  %d > 0; %d < 0; %d == 0\n" %
                           (np.sum(ratio > 0.0), np.sum(ratio < 0.0),
                            np.sum(ratio == 0.0)))
          
          
          
          super_init = super(NormReadDepthTrack, self).__init__
          super_init(ratio, region, options)
          

