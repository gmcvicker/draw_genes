
import sys

import numpy as np
import rpy2.robjects as robjects

from numerictrack import NumericTrack


class ContinuousTrack(NumericTrack):
    """Class for drawing numeric values that are to be displayed as 
    a continuous curve or profile"""

    def __init__(self, values, region, options):
        super(ContinuousTrack, self).__init__(values, region, options)

        if 'smooth' in options:
            smooth_window_sz = int(options['smooth'])
        else:
            smooth_window_sz = 1
        
        if smooth_window_sz > 1:
            if 'smoother' in options:
                smoother = options['smoother']
            else:
                smoother = 'average'
            self.values = self.smooth_values(values, smooth_window_sz, smoother)
        else:
            self.values = values

        self.set_y_range(options)

        if 'n_ticks' in options:
            self.n_ticks = int(options['n_ticks'])
        else:
            self.n_ticks = 3

    
    def smooth_values(self, vals, win_sz, method):
        if method == 'average':
            return self.moving_average(vals, win_sz)
        elif method == 'savitsky-golay':
            return self.savitsky_golay(vals, win_sz)
        else:
            return vals

        
    # adapted from http://www.scipy.org/Cookbook/SavitskyGolay
    def savitsky_golay(self, vals, win_sz, order=2):
        try:
            win_sz = np.abs(np.int(win_sz))
            order = np.abs(np.int(order))
        except ValueError, msg:
            raise ValueError("smoothing window size and order "
                             "must be of type int")

        if win_sz < 1:
            raise TypeError("smoothing window size size must be a "
                            "positive odd number")

        if win_sz < order + 2:
            if order % 2 == 0:
                win_sz = order + 3
            else:
                win_sz = order + 2
            sys.stderr.write("  WARNING: smoothing window size is too small "
                             "for the polynomial's order; setting to minimum "
                             "value of %s.\n" % win_sz)

        if win_sz % 2 != 1:
            win_sz += 1
            sys.stderr.write("  WARNING: smoothing window size must be odd; "
                             "incrementing by one.\n")

        order_range = range(order+1)
        half_window = (win_sz -1) // 2

        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[0]

        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = vals[0] - np.abs(vals[1:half_window+1][::-1] - vals[0])
        lastvals = vals[-1] + np.abs(vals[-half_window-1:-1][::-1] - vals[-1])

        vals = np.concatenate((firstvals, vals, lastvals))
        vals = np.convolve(m, vals, mode='valid')
        vals[vals < 0.0] = 0.0
        return vals

    
    def moving_average(self, vals, win_sz):
        win = np.ones(win_sz, dtype=np.float32)
        # take moving average using win_sz sliding window
        return np.convolve(win / win.sum(), vals, mode='same')

            

    def get_segments(self, vals):
        """Returns set of lists representing contiguous
        segments with the same values."""
        x1 = []
        x2 = []
        y = []

        # region_len = self.region.end - self.region.start + 1
        
        # if vals.size != region_len:
        #     raise ValueError("expected number of values (%d) to match "
        #                      "region len (%d)" % (vals.size, region_len))
        
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
        return (np.array(x1), np.array(x2), np.array(y))




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

        return (np.array(p_x), np.array(p_y))
            
        
    def draw_track(self, r):

        if self.max_val == self.min_val:
            yscale = self.height
        else:
            yscale = self.height / (self.max_val - self.min_val)

        vals = (self.values - self.min_val) * yscale + self.bottom

        # break region into smaller blocks because
        # some programs (e.g. illustrator) don't like it when
        # polygons have too many points
        block_sz = 10000
        for block_start in range(0, vals.size, block_sz):
            block_end = min(block_start + block_sz, vals.size)

            block_vals = vals[block_start:block_end]
            
            # identify contiguous segments with same values
            (x1, x2, y) = self.get_segments(block_vals)

            # new way of drawing: convert contiguous segments
            # to polygon coordinates
            (x, y) = self.get_polygon_coords(x1, x2, y)

            if len(x) > 0:
                r.polygon(robjects.FloatVector(x + block_start),
                          robjects.FloatVector(y),
                          col=self.color,
                          border=self.border_color)
                      
        # old way, 
        # n_seg = len(x1)
        # if n_seg > 0:
        #     r.rect(robjects.FloatVector(x1),
        #            robjects.FloatVector([self.bottom]),
        #            robjects.FloatVector(x2),
        #            robjects.FloatVector(y), col=self.color,
        #            border=robjects.NA_Logical)

        self.draw_y_axis(r, self.n_ticks)


