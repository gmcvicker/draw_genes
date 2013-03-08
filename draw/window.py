
import numpy as np
import re
import rpy2.robjects as robjects



def n_digits(x):
    s = 1
    n_digit = 0

    while s < x:
        s *= 10
        n_digit += 1
    
    return n_digit


def add_commas(x):
    return re.sub(r'(\d{3})(?=\d)', r'\1,', str(x)[::-1])[::-1]




class Window(object):
    """Instances of this class represent genomic windows. Tracks are added
    to the window and track positioning and drawing is handled by this class.
    The details of how the track is drawn is left up to the implementation
    of the Track.draw method."""
    
    def __init__(self, region, margin=0.10, draw_grid=True,
                 vert_lines=[], draw_midline=False, cex=1.0):
        self.region = region
        self.margin = margin
        self.draw_grid = draw_grid
        self.draw_midline = draw_midline
        self.vert_lines = vert_lines
        self.cex = cex
        self.tracks = []

    def add_track(self, track):
        """Add a track to this genome window"""
        self.tracks.append(track)

        

    def draw_gridlines(self, r, top, bottom):
        region_len = self.region.end-self.region.start+1

        # we want to draw lines on nice round numbers
        n_dig = n_digits(region_len)
        if n_dig < 2:
            return

        n_dec = n_dig - 2
        if n_dec < 1:
            n_dec = 1

        grid_sz = 10.0 ** n_dec
        start = np.around(self.region.start, decimals=-n_dec)
        end = self.region.end

        x_left = np.arange(start, end, 2*grid_sz)
        x_right = x_left + grid_sz

        if x_left.size > 0:
            r.rect(robjects.FloatVector(x_left),
                   top,
                   robjects.FloatVector(x_right),
                   bottom,
                   col="grey90",
                   border=robjects.r("NA"))
    
        
    def draw_axis(self, r):
        region_len = self.region.end-self.region.start+1

        # we want to draw lines on nice round numbers
        n_dig = n_digits(region_len)
        if n_dig < 2:
            return

        n_dec = n_dig - 2
        if n_dec < 1:
            n_dec = 1

        step_size = (10.0 ** n_dec) * 2
        
        start = np.around(self.region.start, decimals=-n_dec)
        end = self.region.end + step_size
        ticks = np.arange(start, end, step_size, dtype=np.uint32)

        label_dec = 6 - n_dig
        if label_dec < 0:
            label_dec = 0

        labels = [add_commas(x) for x in ticks]
            
        r.axis(r.c(1), at=robjects.IntVector(ticks),
               labels=labels, **{'cex.axis' : self.cex} )
            
            


    def get_height(self):
        height = 0
        for track in self.tracks:
            height += track.height + self.margin
        return height
    

    def draw(self, r):
        """Plots this genome window"""
        if len(self.tracks) == 0:
            return

        cur_y = -self.margin

        # position each of the tracks and calculate
        # total y dimension
        for track in self.tracks:
            top = cur_y
            bottom = top - track.height
            track.set_position(self.region.start, self.region.end,
                               top, bottom)
            cur_y = bottom - self.margin

        # create an empty plot
        xlim = r.c(self.region.start, self.region.end)

        top = 0
        bottom = cur_y
        ylim = r.c(bottom, top)

        xlab = self.region.chrom.name + " position"

        r.plot(r.c(0), r.c(0), type="n", xlim=xlim, ylim=ylim,
               yaxt="n", xaxt="n", xlab=xlab, ylab="", bty="n",
               **{"mar" : r.c(5.1, 0.1, 0.1, 0.1)})

        self.draw_axis(r)
        
        if self.draw_grid:
            self.draw_gridlines(r, top, bottom)

        # draw a vertical line at the midpoint
        if self.draw_midline:
            region_mid = (self.region.start + self.region.end) / 2
            r.lines(robjects.FloatVector([region_mid, region_mid]),
                    robjects.FloatVector([top, bottom]), col="grey70")

        # draw vertical lines where specified
        for x in self.vert_lines:
            r.lines(robjects.FloatVector([x, x]), robjects.FloatVector([top, bottom]),
                    col="grey")

        # now draw each of the tracks
        for track in self.tracks:
            track.draw(r)

