
import sys

from track import Track


import genome.db
import rpy2.robjects as robjects



class StateTrack(Track):
    """Class for drawing set of discrete non-overlapping states"""

    def __init__(self, region, options):
        super(StateTrack, self).__init__(region, options)

        self.n_state = int(options['n_state'])
        
        if self.n_state < 1:
            raise ValueError("expected n_state to be >= 1")

        # assign a color to each state
        self.state_colors = {}
        for i in range(self.n_state + 1):
            color_label = "state_color_%d" % i
            if color_label in options:
                self.state_colors[i] = options[color_label].replace('"', '')
            else:
                self.state_colors[i] = "white"

        # assign a label to each state
        self.state_labels = {}
        for i in range(self.n_state + 1):
            label = "state_label_%d" % i
            if label in options:
                self.state_labels[i] = options[label]
            else:
                self.state_labels[i] = ""
        
        self.track_name = options['track']
        self.features = []

        # get GenomeDB instance
        gdb = options['gdb']
        track = gdb.open_track(self.track_name)

        self.features = self.__create_features(region, track)

        track.close()


    def __create_features(self, region, track):
        vals = track.get_nparray(region.chrom, region.start, region.end)

        region_len = region.end - region.start + 1
        
        features = []

        cur_feat = None
        for i in range(region_len):
            if (cur_feat is None) or (vals[i] != cur_feat.state_id):
                if cur_feat:
                    # update end of previous feature
                    cur_feat.end = region.start + i - 1
                
                cur_feat = genome.coord.Coord(region.chrom,
                                              region.start + i,
                                              region.start + i)
                cur_feat.state_id = vals[i]

                features.append(cur_feat)

        if cur_feat:
            # update end of last feature
            cur_feat.end = region.end

        return features
        
            



    def draw_track(self, r):
        feat_height = 0.5 * self.height
        margin_height = self.height - feat_height

        region_len = self.region.end - self.region.start + 1.0

        for feat in self.features:
            top = self.top - margin_height/2
            bottom = top - feat_height

            # color based on state number of feature
            if feat.state_id in self.state_colors:
                color = self.state_colors[feat.state_id]
            else:
                color = "grey50"
            
            r.rect(feat.start, bottom, feat.end, top, col=color,
                   border=robjects.r("NA"))

        # draw a label for the entire track
        self.draw_track_label(r)
        
        # now draw labels on top of features            
        label_offsets = (0.5, 0.0, -0.5)
        offset_idx = 0

        for feat in self.features:
            label_offset = label_offsets[offset_idx]
            offset_idx += 1
            if offset_idx >= len(label_offsets):
                offset_idx = 0
        
            if feat.state_id in self.state_labels:
                # add label to state
                label = self.state_labels[feat.state_id]
                mid_y = (top + bottom) * 0.5 + label_offset
                mid_x = (feat.start + feat.end) * 0.5
                r.text(x=mid_x, y=mid_y, labels=r.c(label),
                       cex=self.cex)
            else:
                sys.stderr.write("no label for state %d\n" % feat.state_id)
            
