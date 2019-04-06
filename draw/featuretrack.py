
import sys

from .track import Track

import genome.coord
import rpy2.robjects as robjects

MIN_FEAT_LEN = 5

class FeatureTrack(Track):
    """Class for drawing a single transcript"""

    def __init__(self, region, options):        
        super(FeatureTrack, self).__init__(region, options)

        if 'draw_labels' in options:
            self.draw_labels = self.parse_bool_str(options['draw_labels'])
        else:
            self.draw_labels = False
        


        if 'cex' in options:
            self.cex = float(options['cex'])
        else:
            self.cex = 1.0
        
        
        self.track_name = options['track']
        self.features = []

        # get overlapping features from database
        gdb = options['gdb']
        track = gdb.open_track(self.track_name)
        table = track.h5f.getNode("/" + self.region.chrom.name)

        qry = "(start <= %d) & (end >= %d)" % (self.region.end,
                                               self.region.start)
        for row in table.where(qry):
            feat = genome.coord.Coord(region.chrom, row['start'],
                                      row['end'], strand=row['strand'],
                                      score=row['score'], name=row['name'])
            self.features.append(feat)


        # set padding here if labels are used...
        if self.draw_labels:
            padding = region.length() * 0.1
        else:
            padding = region.length() * 0.025
            
        self.assign_feature_rows(self.features, use_strands=False,
                                 padding=padding)
        
        if self.height <= 0.0:
            # assign height based on how many rows of features
            # there are
            self.height = float(self.n_row) * 0.5

        track.close()

    
    def draw_track(self, r):
        if self.n_row > 0:
            y_scale = self.height / float(self.n_row)
        else:
            y_scale = 1.0
        feat_height = y_scale * 0.8
        margin_height = y_scale - feat_height
        region_len = self.region.end - self.region.start + 1.0

        for feat in self.features:
            row = self.row_assignment[feat]
            top = self.top - ((row * feat_height) + (margin_height)*row)
            bottom = top - feat_height

            if feat.strand == 1:
                color = self.fwd_color
                border_color = self.fwd_border_color
                label = "> %s" % feat.name
            elif feat.strand == -1:
                color = self.rev_color
                border_color = self.rev_border_color
                label = "< %s" % feat.name
            else:
                color = self.color
                border_color = self.border_color
                label = feat.name
                
            r.rect(feat.start - 0.5, bottom, feat.end + 0.5, top,
                   col=color, border=border_color)

            if self.draw_labels:     
                # draw feature label?
                mid = (top + bottom) * 0.5
                label_len = float(len(label))
                offset = (0.0075 * label_len) * region_len
                r.text(x=(feat.end + offset), y=mid,
                       labels=r.c(label), col=color, cex=self.cex)


