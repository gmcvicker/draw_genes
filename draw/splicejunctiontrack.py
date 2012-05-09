import sys
import math

import numpy as np
import rpy2.robjects as robjects

from track import Track



class SpliceJunctionTrack(Track):
    def __init__(self, gdb, region, height=None,
                 border=False, color="black"):
        sys.stderr.write("creating splice junction track\n")

        super(SpliceJunctionTrack, self).__init__(height=height, border=border)
        self.color = color
        self.region = region
        self.gdb = gdb

        spliced_read_adp = self.gdb.get_adaptor("spliced_read")
        reads = spliced_read_adp.fetch_by_region(self.region)

        self.junctions = []
        for read in reads:
            best = read.best_junctions()
            for junction in best:
                junction.read_count = read.read_count
            self.junctions.extend(best)


        padding = region.length() * 0.01
        self.assign_feature_rows(self.junctions, padding=padding,
                                 use_strands=True)

        if self.height is None:
            self.height = float(self.n_row) * 0.2

        

    def draw_track(self, r):
        if len(self.junctions) < 1:
            return

        row_width = float(self.height) / float(self.n_row)
        row_spacing = 0.1 * row_width
        region_len = self.right - self.left + 1
        text_offset = 0.015 * region_len

        for junction in self.junctions:
            if junction.read_count >= 5:
                lwd = 1.5
            else:
                lwd = 1
                
            # color of junction depends on whether it is canonical or not
            if junction.is_canonical:
                color = self.color
            elif junction.is_non_canonical:
                color = "green"
            else:
                color = "grey50"

            row = self.row_assignment[junction]

            top = self.top - (float(row)*row_width)
            bottom = top - row_width + row_spacing
            midpoint = (top + bottom)*0.5            

            # orient the junction if possible
            if junction.strand == 1:
                r.segments(x0=junction.start, y0=top,
                           x1=junction.end, y1=top, col=color,
                           lwd=lwd)
            elif junction.strand == -1:
                r.segments(x0=junction.start, y0=bottom,
                           x1=junction.end, y1=bottom, col=color,
                           lwd=lwd)
            else:
                r.segments(x0=junction.start, y0=midpoint,
                           x1=junction.end, y1=midpoint, col=color,
                           lwd=lwd)

            x0 = [junction.start, junction.end]
            x1 = [junction.start, junction.end]
            y0 = [bottom, bottom]
            y1 = [top, top]
            
            r.segments(x0=robjects.FloatVector(x0),
                       y0=robjects.FloatVector(y0),
                       x1=robjects.FloatVector(x1),
                       y1=robjects.FloatVector(y1), col=color,
                       lwd=lwd)

            # add text giving read count 

            r.text(x=(junction.end + text_offset), y=midpoint,
                   labels=str(junction.read_count), col="black",
                   cex=0.5)
            
        
