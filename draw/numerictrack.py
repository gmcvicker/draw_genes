
import sys
import numpy as np

import rpy2.robjects as robjects

from .track import Track


class NumericTrack(Track):
    """This is an abstract BaseClass containing functions that are in common
    to tracks with numeric data such as drawing and labeling the y axis."""

    
    def __init__(self, values, region, options):
        super(NumericTrack, self).__init__(region, options)

        self.values = values
        
        if 'n_ticks' in options:
            self.n_ticks = int(options['n_ticks'])
        else:
            self.n_ticks = 3


    
    def set_y_range(self, options):
        """Sets the maximum and minimum values of the y-axis"""
        # first set maximum and minimum values to range of data
        is_nan = np.isnan(self.values)

        if np.any(~is_nan):
            self.max_val = np.max(self.values[~is_nan])
            self.min_val = np.min(self.values[~is_nan])
        else:
            self.max_val = 0.0
            self.min_val = 0.0

        # Make sure that y-axis maximum and minimum are at least
        # those specified by soft_min_val and soft_max_val.
        # These are 'soft' limits in that the axis is allowed to
        # grow beyond them to accomodate larger/smaller values.
        if "soft_max_val" in options:
            # this is the minimum maximum y-axis value
            if float(options['soft_max_val']) > self.max_val:
                self.max_val = float(options['soft_max_val'])

        if "soft_min_val" in options:
            #  maximum allowed minimum value for y-axis
            if float(options['soft_min_val']) < self.min_val:
                self.min_val = float(options['soft_min_val'])

        # If max and min values are specified these are 'hard'
        # limits. Just set the range of the x-axis to these values
        if 'max_val' in options:
            self.max_val = float(options['max_val'])

        if 'min_val' in options:
            self.min_val = float(options['min_val'])


            
    def draw_y_axis(self, r, n_ticks=3):
        region_len = self.region.end - self.region.start + 1
        tick_width = region_len * 0.005
        label_scale = region_len * 0.0075

        if n_ticks < 2:
            return 

        # determine where tick marks should be
        x0 = [self.region.start] * n_ticks
        x1 = [self.region.start - tick_width] * n_ticks
        span = (self.max_val - self.min_val) / (n_ticks-1)

        if span >= 10.0:
            # use integers
            span = round(span)
        
        y = np.arange(n_ticks) * span + self.min_val

        if self.max_val == self.min_val:
            yscale = self.height
        else:
            yscale = self.height / (self.max_val - self.min_val)

        y_transform = (y - self.min_val) * yscale + self.bottom

        # draw tick marks
        r.segments(x0=robjects.FloatVector(x0),
                   x1=robjects.FloatVector(x1),
                   y0=robjects.FloatVector(y_transform),
                   y1=robjects.FloatVector(y_transform))

        # draw y-axis line
        r.segments(x0=robjects.FloatVector([x0[0]]),
                    x1=robjects.FloatVector([x0[0]]),
                    y0=robjects.FloatVector([y_transform[0]]),
                    y1=robjects.FloatVector([y_transform[-1]]))

        # draw labels beside tick marks
        if span >= 10.0:
            labels = [str(int(val)) for val in y]
        elif span >= 0.5:
            labels = ["%.1f" % val for val in y]
        elif span > 0.01:
            labels = ["%.2f" % val for val in y]
        else:
            labels = [str("%.2e" % val) for val in y]
        
        # x_label = self.region.start - tick_width * 0.5
        x_label = self.region.start

        r.text(x=x_label, y=robjects.FloatVector(y_transform),
               pos=2, cex=self.cex * 0.75,
               labels=robjects.StrVector(labels))


    
    
    
