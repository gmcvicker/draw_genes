
import sys

from track import Track

import genome.db
import rpy2.robjects as robjects


class SegmentTrack(Track):
    """Class for drawing a single transcript"""

    def __init__(self, region, options):        
        super(SegmentTrack, self).__init__(region, options)

        if 'color' in options:
            self.color = options['color'].replace('"', '')
        else:
            self.color = "black"

        self.track_name = options['track']
        self.features = []

        # get overlapping features from database
        gdb = options['gdb']
        track = gdb.open_track(self.track_name)
        table = track.h5f.getNode("/" + self.region.chrom.name)

        qry = "(start <= %d) & (end >= %d)" % (self.region.end,
                                               self.region.start)
        for row in table.where(qry):
            feat = genome.coord.Coord(region.chrom, row['start'], row['end'])
            self.features.append(feat)


        if self.height <= 0.0:
            self.height = 1.0
        
    
    def draw_track(self, r):
        feat_height = 0.5 * self.height
        margin_height = self.height - feat_height

        if len(self.features) == 0:
            # no features
            return

        feat_left = []
        feat_right = []
        for feat in self.features:
            feat_left.append(feat.start)
            feat_right.append(feat.end)

        feat_top = self.top - margin_height/2
        feat_bottom = feat_top - feat_height

        # draw rectangle for each feature
        r.rect(robjects.FloatVector(feat_left),
               feat_bottom,
               robjects.FloatVector(feat_right),
               feat_top,
               col=self.color, border=robjects.IntVector("NA"))


