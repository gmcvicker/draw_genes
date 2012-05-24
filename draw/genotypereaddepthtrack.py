import sys

import numpy as np

import genome.db
import genome.wig

from continuoustrack import ContinuousTrack

class GenotypeReadDepthTrack(ContinuousTrack):     
     def __init__(self, region, options):
          track_attrib_name = options['track_attrib']
          total_attrib_name = options['total_attrib']

          track_names = region.attrib[track_attrib_name].split(",")
          total_count = int(region.attrib[total_attrib_name])

          if track_names[0] == ".":
               # '.' is just a place holder
               track_names = []

          if total_count > 0 and len(track_names) == 0:
               raise ValueError("no track names but total count > 0")

          if len(track_names) > 0 and total_count == 0:
               raise ValueError("track names but total count == 0")

          values = np.zeros(region.length())

          if total_count == 0:
               values[:] = np.nan
          else:
               # get GenomeDB instance
               gdb = options['gdb']

               # sum values from all tracks
               for track_name in track_names:
                    track = gdb.open_track(track_name)
                    
                    values[:] += track.get_nparray(region.chrom,
                                                   start=region.start,
                                                   end=region.end)
                    track.close()

               # rescale by total number of sequenced reads for these
               # individuals
               scale = 1e9 / total_count
               values *= scale
          
          log_scale = False
          if "log_scale" in options:
               log_scale = self.parse_bool_str(options['log_scale'])

          if log_scale:
               # add one to values, but avoid possible overflow of
               # 8 bit values
               f = values < 255
               values[f] += 1
               values = np.log2(values)

          super_init = super(GenotypeReadDepthTrack, self).__init__          
          super_init(values, region, options)
