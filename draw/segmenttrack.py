
import sys

from track import Track

import genome.db
import rpy2.robjects as robjects

import numpy as np

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

        if 'track_type' in options:
            track_type = options['track_type']
        else:
            track_type = 'table'

        if track_type == 'table':
            self.add_table_features(track)
        elif track_type == 'flags':
            self.add_flag_features(track)
        else:
            raise ValueError("unknown track type '%s' expected 'table' or 'flags'")


        if self.height <= 0.0:
            self.height = 1.0


    def add_table_features(self, track):        
        table = track.h5f.getNode("/" + self.region.chrom.name)

        qry = "(start <= %d) & (end >= %d)" % (self.region.end,
                                               self.region.start)
        for row in table.where(qry):
            feat = genome.coord.Coord(self.region.chrom, row['start'], row['end'])
            self.features.append(feat)


    
    def add_flag_features(self, track):
        vals = track.get_nparray(self.region.chrom, self.region.start, self.region.end)

        if vals[0] == 1:
            # record start of region as start of first segment
            starts = [self.region.start]
        else:
            starts = []

        d = np.diff(vals)
        # starts are at transitions from 0 => 1 
        starts.extend(np.where(d == 1.0)[0] + self.region.start + 1)
        # and ends are at transitions from 1 => 0
        ends = list(np.where(d == -1.0)[0] + self.region.start + 1)

        if vals[-1] == 1:
            ends.append(self.region.end)

        if len(starts) != len(ends):
            raise ValueError("Expected same number of starts and ends")

        for s, e in zip(starts, ends):
            feat = genome.coord.Coord(self.region.chrom, s, e)
            self.features.append(feat)
        
            
        
        
            
    
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
               col=self.color, border=self.border_color)


