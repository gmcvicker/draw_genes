
import sys

import numpy as np
import rpy2.robjects as robjects

import genome.db
from continuoustrack import ContinuousTrack


class BaseLLRTrack(ContinuousTrack):
    """LLRs are a type of continuous track that are symmetric about 0.
    Draws negative values below an axis at 0.0 (rather than setting the
    axis at the minimum value for the track and drawing all values above
    it)"""
    def __init__(self, values, region, options):
        super_init = super(BaseLLRTrack, self).__init__
        super_init(values, region, options)

        if "pos_color" in options:
            self.pos_color = options['pos_color'].replace('"', '')
        else:
            self.pos_color = self.color
        
        if "neg_color" in options:
            self.neg_color = options['neg_color'].replace('"', '')
        else:
            self.neg_color = self.color
        


    
    def draw_track(self, r):
        if self.border:
            self.draw_border(r)

        if self.min_val > 0.0:
            min_val = 0.0
        else:
            min_val = self.min_val

        if self.max_val < 0.0:
            max_val = 0.0
        else:
            max_val = self.max_val

        if max_val == min_val:
            yscale = 1.0
        else:
            yscale = self.height / (self.max_val - min_val)

        # put the axis at 0
        axis = (-min_val * yscale) + self.bottom
        vals = (self.values * yscale) + axis

        (x1, x2, y) = self.get_segments(vals)

        n_seg = len(x1)
        if n_seg > 0:
            
            x1 = np.array(x1, dtype=np.float32)
            # add one because drawn coordinates are "between"
            # start and end
            # x2 = np.array(x2, dtype=np.float64) + 1.0
            x2 = np.array(x2, dtype=np.float32)
            y = np.array(y, dtype=np.float32)

            pos_vals = ((y >= axis) & (~np.isnan(y)))
            neg_vals = ((y < axis) & (~np.isnan(y)))

            # color positive and negative values separately,
            # don't draw 0 values
            if np.any(pos_vals):
                # draw values above 0
                (p_x, p_y) = self.get_polygon_coords(x1[pos_vals],
                                                     x2[pos_vals],
                                                     y[pos_vals],
                                                     axis=axis)
                r.polygon(robjects.FloatVector(p_x),
                          robjects.FloatVector(p_y),
                          col=self.pos_color, border=self.pos_color)

                
            if np.any(neg_vals):
                # draw values below zero
                (p_x, p_y) = self.get_polygon_coords(x1[neg_vals],
                                                     x2[neg_vals],
                                                     y[neg_vals],
                                                     axis=axis)
                r.polygon(robjects.FloatVector(p_x),
                          robjects.FloatVector(p_y),
                          col=self.neg_color, border=self.neg_color)

            self.draw_y_axis(r, self.n_ticks)

                


