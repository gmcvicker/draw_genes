import numpy as np
import rpy2.robjects as robjects

from track import Track
from numerictrack import NumericTrack


class PointsTrack(NumericTrack):
    """This class is for drawing points defined at a subset of genomic
    positions. Sites in the track that are set to nan are ignored.
    This track can be used to display p-values for SNPs etc.
    """

    def __init__(self, region, options):

        # retrieve values from genome db. If desired could later
        # modify this to have 'source' option and allow values to be
        # read from other file types
        track_name = options['track']
        gdb = options['gdb']
        track = gdb.open_track(track_name)
        values = track.get_nparray(region.chrom, start=region.start,
                                   end=region.end)

        # get positions of defined values:
        defined_idx = np.where(~np.isnan(values)[0])
        self.pos = self.region.start + defined_idx

        self.values = values[defined_idx]

        if 'neg_log_transform' in options and \
          self.parse_bool_str(options['neg_log_transform']):
          self.values = -np.log10(self.values)

        super(PointsTrack, self).__init__(region, options)
        
        self.set_y_range(options)


        
    def draw_track(self, r):
        yscale = self.height / (self.max_val - self.min_val)

        vals = (self.values - self.min_val) * yscale + self.bottom

        if len(self.values) > 0:
            r.points(robjects(FloatVector(self.pos))
                     robjects(FloatVector(self.values)),
                     col=self.color,
                     bg=self.color
                     self.pch=21)

        self.draw_y_axis(r, self.n_ticks)
        
                              
        

                    

        
    
