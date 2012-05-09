import sys

import numpy as np

import genome.db
import genome.wig
import genome.trackstat

from continuoustrack import ContinuousTrack

class ReadDepthTrack(ContinuousTrack):     
     def __init__(self, region, options):

          source = "gdb"

          if "source" in options:
               source = options['source']

          if source == "gdb":
               track_name = options['track']
               gdb = genome.db.GenomeDB()
               track = gdb.open_track(track_name)
               values = track.get_nparray(region.chrom, start=region.start,
                                          end=region.end)

               if "scale_factor" in options:
                    scale_factor = float(options['scale_factor'])                    
                    stat = gdb.get_track_stat(track)
                    scale = scale_factor / float(stat.sum)
                    values = values * scale
                    sys.stderr.write("  total reads %d, using scale %.3f\n" %
                                     (stat.sum, scale))
               track.close()

          if source == "wig":
               path = options['path']
               values = genome.wig.read_ints(path, region)
          
          log_scale = False
          if "log_scale" in options:
               log_scale = self.parse_bool_str(options['log_scale'])

          if log_scale:
               # add one to values, but avoid possible overflow of
               # 8 bit values
               f = values < 255
               values[f] += 1
               values = np.log2(values)

          super_init = super(ReadDepthTrack, self).__init__          
          super_init(values, region, options)
          

