
import sys

import numpy as np
import rpy2.robjects as robjects

from track import Track


def smooth_values(vals, win_sz):
    win = np.ones(win_sz, dtype=np.float32)
    # take moving average using win_sz sliding window
    return np.convolve(win / win.sum(), vals, mode='same')
    
    

class ContinuousTrack(Track):

    def __init__(self, values, region, options):
        super(ContinuousTrack, self).__init__(region, options)

        if 'color' in options:
            self.color = options['color'].replace('"', '')
        else:
            # use a default color
            self.color = "#E41A1C"

        if 'smooth' in options:
            smooth_window_sz = int(options['smooth'])
        else:
            smooth_window_sz = 1
        
        if smooth_window_sz > 1:
            self.values = smooth_values(values, smooth_window_sz)
        else:
            self.values = values

        self.set_y_range(options)


        if 'n_ticks' in options:
            self.n_ticks = int(options['n_ticks'])
        else:
            self.n_ticks = 3

        if 'draw_border' in options:
            self.draw_border = self.parse_bool_str(options['draw_border'])
        else:
            self.draw_border = True


    def set_y_range(self, options):
        """Sets the maximum and minimum values of the y-axis"""
        # first set maximum and minimum values to range of data
        is_nan = np.isnan(self.values)

        if np.any(~is_nan):
            self.max_val = np.max(self.values[~is_nan])
        else:
            self.max_val = 0.0

        if np.any(~is_nan):
            self.min_val = np.min(self.values[~is_nan])
        else:
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
            


    def get_segments(self, vals):
        """Returns set of lists representing contiguous
        segments with the same values."""
        x1 = []
        x2 = []
        y = []

        region_len = self.region.end - self.region.start + 1
        
        if vals.size != region_len:
            raise ValueError("expected number of values (%d) to match "
                             "region len (%d)" % (vals.size, region_len))
        
        cur_start = cur_end = cur_val = None
        pos = self.region.start-1

        n_undef = 0
        n_def = 0
        
        for val in vals:
            pos += 1
            
            if np.isnan(val):
                # undefined region
                n_undef += 1
                if cur_val is not None:
                    # end current segment
                    cur_end = pos-1
                    y.append(cur_val)
                    x1.append(cur_start)
                    x2.append(cur_end+1)
                cur_start = cur_start = cur_val = None
            else:
                n_def += 1
                if cur_val is None or cur_val != val:
                    # start a new segment...
                    if cur_start is not None:
                        # ..but end current segment first
                        cur_end = pos-1
                        y.append(cur_val)
                        x1.append(cur_start)
                        # add 1 because drawn coordinates are "between"
                        # start/end
                        x2.append(cur_end+1)

                    cur_start = pos
                    cur_val = val
        if cur_val is not None:
            # end final segment
            cur_end = pos - 1
            y.append(cur_val)
            x1.append(cur_start)
            x2.append(cur_end+1)

            #sys.stderr.write("%d segments for regions of size %d "
            #             "(undef=%d, def=%d)\n" %
            #             (len(x1), region_len, n_undef, n_def))
        return (x1, x2, y)


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



    def get_polygon_coords(self, x1, x2, y, axis=None):
        if len(x1) < 1:
            # no segments in this region
            return ([], [])

        if not axis:
            axis = self.bottom

        # start with the first segment
        p_x = [x1[0], x1[0], x2[0]]
        p_y = [axis, y[0], y[0]]
        prev_x = x2[0]
        
        for i in range(1, len(x1)):
            if x1[i] > (prev_x + 1):
                # discontiguous segments
                # finish prev segment
                p_x.append(prev_x)
                p_y.append(axis)

                # start a new segment
                p_x.append(x1[i])
                p_y.append(axis)

            p_x.append(x1[i])
            p_y.append(y[i])

            p_x.append(x2[i])
            p_y.append(y[i])

            prev_x = x2[i]

        # finish final segment
        p_x.append(prev_x)
        p_y.append(axis)

        return (p_x, p_y)
            
        
    def draw_track(self, r):
        yscale = self.height / (self.max_val - self.min_val)

        vals = (self.values - self.min_val) * yscale + self.bottom

        # identify contiguous segments with same values
        (x1, x2, y) = self.get_segments(vals)

        # new way of drawing: convert contiguous segments
        # to polygon coordinates
        (x, y) = self.get_polygon_coords(x1, x2, y)

        if self.draw_border:
            border = "black"
        else:
            border = robjects.r("NA")
            
        if len(x) > 0:
            r.polygon(robjects.FloatVector(x),
                      robjects.FloatVector(y),
                      col=self.color,
                      border=border)
                      
        # old way, 
        # n_seg = len(x1)
        # if n_seg > 0:
        #     r.rect(robjects.FloatVector(x1),
        #            robjects.FloatVector([self.bottom]),
        #            robjects.FloatVector(x2),
        #            robjects.FloatVector(y), col=self.color,
        #            border=robjects.NA_Logical)

        self.draw_y_axis(r, self.n_ticks)


