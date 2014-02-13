import numpy as np
import rpy2.robjects as robjects
import sys

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
        defined_idx = np.where(~np.isnan(values))[0]
        self.pos = region.start + defined_idx

        values = values[defined_idx]
        values[values < 1e-30] = 1e-30

        self.set_colors(options)
        
        if len(values) > 0:
            sys.stderr.write("range: %g %g\n" % (np.min(values), np.max(values)))
        else:
            sys.stderr.write("WARNING: no defined values in region\n")
        

        if "threshold" in options:
            self.threshold = float(options['threshold'])

            if "below_thresh_color" in options:
                self.below_thresh_color = options['below_thresh_color'].replace('"', '')
            else:
                self.below_thresh_color = self.color

            if "above_thresh_color" in options:
                self.above_thresh_color = options['above_thresh_color'].replace('"', '')
            else:
                self.above_thresh_color = self.color

            if "draw_thresh_line" in options:
                self.draw_thresh_line = self.parse_bool_str(options['draw_thresh_line'])
            else:
                self.draw_thresh_line = False
        else:
            self.threshold = None
            self.below_thresh_color = self.color
            self.above_thresh_color = self.color
            self.draw_thresh_line = False
            
        if 'neg_log_transform' in options:
            self.neg_log_transform = self.parse_bool_str(options['neg_log_transform'])
        else:
            self.neg_log_transform = False

        if self.neg_log_transform:
          values = -np.log10(values)

          if self.threshold is not None:
              self.threshold = -np.log10(self.threshold)

        super(PointsTrack, self).__init__(values, region, options)

        self.set_y_range(options)


        
    def draw_track(self, r):

        sys.stderr.write("DRAWING POINTS\n")
        
        yscale = self.height / (self.max_val - self.min_val)

        vals = (self.values - self.min_val) * yscale + self.bottom

        if len(self.values) > 0:
            r.points(robjects.FloatVector(self.pos),
                     robjects.FloatVector(vals),
                     col=self.color,
                     bg=self.color, cex=0.5,
                     pch=21)

            # draw line at 0
            zero_val = -self.min_val * yscale + self.bottom
            r.lines(robjects.FloatVector([self.region.start, self.region.end]),
                    robjects.FloatVector([zero_val, zero_val]),
                    col="grey50")

            if self.threshold is not None:
                if self.draw_thresh_line:
                    # draw line showing threshold
                    thresh_val = (self.threshold - self.min_val) * yscale + self.bottom
                    r.lines(robjects.FloatVector([self.region.start, self.region.end]),
                            robjects.FloatVector([thresh_val, thresh_val]),
                            lty=2, col="red")

                # redraw points below threshold a different color
                if self.neg_log_transform:
                    below_thresh = self.values >= self.threshold
                    above_thresh = self.values <= self.threshold
                else:
                    below_thresh = self.values <= self.threshold
                    above_thresh = self.values >= self.threshold

                if np.any(below_thresh):
                    r.points(robjects.FloatVector(self.pos[below_thresh]),
                             robjects.FloatVector(vals[below_thresh]),
                             col=self.below_thresh_color,
                             bg=self.below_thresh_color, cex=0.5,
                             pch=21)
                    
                if np.any(above_thresh):
                    sys.stderr.write("drawing points above threshold:%s\n" %
                                     ",".join(["%d" % x for x in self.pos[above_thresh]]))
                    r.points(robjects.FloatVector(self.pos[above_thresh]),
                             robjects.FloatVector(vals[above_thresh]),
                             col=self.above_thresh_color,
                             bg=self.above_thresh_color, cex=0.5,
                             pch=21)
                    

        self.draw_y_axis(r, self.n_ticks)
        
                              
        

                    

        
    
