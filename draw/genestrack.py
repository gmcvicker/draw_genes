
import sys

from genome import coord
from track import Track
from transcripttrack import TranscriptTrack

        
class GenesTrack(Track):
    """Class for drawing all of the genes in a region"""
    
    def __init__(self, all_transcripts, region, options):
        # call superclass constructor
        super(GenesTrack, self).__init__(region, options)

        self.transcripts = all_transcripts

        self.n_fwd_rows = None
        self.n_rev_rows = None
        self.row_assignment = None
        self.color = options['color']
        self.utr_color = options['utr_color']
        
        # get all transcripts that overlap region
        self.overlap_trs = coord.get_coord_overlaps(self.region,
                                                    self.transcripts,
                                                    use_strand=False)

        # assign rows to the transcripts
        padding = region.length() * 0.2
        self.assign_feature_rows(self.overlap_trs, padding=padding)

        if self.height <= 0.0:
            # assign height based on how many rows of transcripts
            # there are
            self.height = float(self.n_row) * 0.5


    def draw_track(self, r):
        # make margin 20% the width of a transcript
        y_scale  = self.height / float(self.n_row)
        tr_height = y_scale * 0.80
        margin_height = y_scale - tr_height

        for tr in self.overlap_trs:
            # give each transcript its own track
            tr_options = {'color' : self.color,
                          'utr_color' : self.utr_color,
                          'border' : 'false',
                          'height' : str(tr_height)}
            
            tr_track = TranscriptTrack(tr, self.region, tr_options)

            # position transcript track
            row = self.row_assignment[tr]
            top = self.top - ((row * tr_height) + (margin_height*row))
            bottom = top - tr_height
            tr_track.set_position(self.region.start, self.region.end,
                                  top, bottom)

            # draw transcript
            tr_track.draw_track(r)
