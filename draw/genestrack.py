
import sys

from genome import coord
import genome.gene
from .track import Track
from .transcripttrack import TranscriptTrack

        
class GenesTrack(Track):
    """Class for drawing all of the genes in a region"""
    
    def __init__(self, all_genes, region, options):
        # call superclass constructor
        super(GenesTrack, self).__init__(region, options)

        self.genes = all_genes
        self.n_fwd_rows = None
        self.n_rev_rows = None
        self.row_assignment = None
        self.color = options['color']
        self.utr_color = options['utr_color']
        self.longest_isoform_only = self.parse_bool_str(options['longest_isoform_only'])

        if 'draw_label' in options:
            self.draw_label = self.parse_bool_str(options['draw_label'])
        else:
            self.draw_label = True

        sys.stderr.write("%d genes total\n" % len(all_genes))
            
        # get all genes that overlap region
        self.overlap_genes = coord.get_coord_overlaps(self.region,
                                                      self.genes,
                                                      use_strand=False)

        sys.stderr.write("%d genes overlap region\n" % len(self.overlap_genes))

        self.overlap_trs = []
        for g in self.overlap_genes:
            if self.longest_isoform_only:
                # only longest transcript for each gene
                self.overlap_trs.append(g.get_longest_transcript())
            else:
                # use all transcripts
                self.overlap_trs.extend(gene.transcripts)
        sys.stderr.write("%d transcripts overlap region" % len(self.overlap_trs))
        
        # assign rows to the transcripts
        if self.draw_label:
            padding = region.length() * 0.2
        else:
            padding = region.length() * 0.01
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

            if self.draw_label:
                draw_label_str = "true"
            else:
                draw_label_str = "false"
            
            tr_options = {'color' : self.color,
                          'utr_color' : self.utr_color,
                          'border' : 'false',
                          'height' : str(tr_height),
                          'draw_label' : draw_label_str}
            
            tr_track = TranscriptTrack(tr, self.region, tr_options)

            # position transcript track
            row = self.row_assignment[tr]
            top = self.top - ((row * tr_height) + (margin_height*row))
            bottom = top - tr_height
            tr_track.set_position(self.region.start, self.region.end,
                                  top, bottom)

            # draw transcript
            tr_track.draw_track(r)
