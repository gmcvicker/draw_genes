
from track import Track


import numpy as np
import rpy2.robjects as robjects


class TranscriptTrack(Track):
    """Class for drawing a single transcript"""

    def __init__(self, transcript, region, options):
        super(TranscriptTrack, self).__init__(region, options)
        self.transcript = transcript
        self.color = options['color'].replace('"', '')
        self.utr_color = options['utr_color'].replace('"', '')

        if 'draw_label' in options:
            self.do_label = self.parse_bool_str(options['draw_label'])
        else:
            self.do_label = True

    
    def draw_coding_region(self, r, coord):
        # add 1 b/c drawn coordinates are "between" start and end
        x_coord = r.c(coord.start, coord.end+1,
                      coord.end+1, coord.start)
        y_coord = r.c(self.top, self.top,
                      self.bottom, self.bottom)
        
        r.polygon(x_coord, y_coord, border=self.color,
                  col=self.color)

    def draw_non_coding_region(self, r, coord):

        # make non-coding portions of exons
        # 2/3 height of coding portions
        h = self.top - self.bottom
        mid = (self.top + self.bottom) * 0.5
        nc_top = mid + h*0.33
        nc_bottom = mid - h*0.33

        # add one because drawn coordinates are "between" start and end
        x_coord = r.c(coord.start, coord.end+1,
                      coord.end+1, coord.start)
        y_coord = r.c(nc_top, nc_top,
                      nc_bottom, nc_bottom)

        r.polygon(x_coord, y_coord, border=self.color,
                  col=self.utr_color)
        

    def draw_exon(self, r, exon):
        cds_start = self.transcript.cds_start
        cds_end = self.transcript.cds_end

        if(self.transcript.is_coding() and
           cds_start <= exon.end and
           cds_end >= exon.start):
            # at least some portion of this exon is coding

            coord = exon.copy()

            if cds_start > exon.start:
                cod_start = cds_start
                # draw non-coding region before cds start
                cod_end = cds_end
                coord.start = exon.start
                coord.end = cds_start-1
                self.draw_non_coding_region(r, coord)
            else:
                cod_start = exon.start
            
            if cds_end < exon.end:
                # draw non-coding region after cds end
                cod_end = cds_end
                # add 1 because drawn region is between start and end
                coord.start = cds_end+1
                coord.end = exon.end
                self.draw_non_coding_region(r, coord)
            else:
                cod_end = exon.end
                            
            # draw coding portion:
            coord.start = cod_start
            coord.end = cod_end

            self.draw_coding_region(r, coord)
        else:
            # this is a non-coding exon
            self.draw_non_coding_region(r, exon)
            

    def draw_intron(self, r, intron):
        # add 1 because drawn region is "between" start and end
        x = r.c(intron.start, intron.end + 1)
        midpoint = (self.top + self.bottom) * 0.5
        y = r.c(midpoint, midpoint)
        r.lines(x, y, col=self.color)


        # draw arrows indicating direction of transcription
        region_len    = self.region.end - self.region.start + 1.0
        arrow_width   = region_len * 0.005
        arrow_height  = self.height * 0.33
        arrow_spacing = arrow_width
        arrow_top     = midpoint + arrow_height * 0.5
        arrow_bottom  = midpoint - arrow_height * 0.5

        if self.transcript.strand == 1:
            arrow_start = np.arange(intron.start + arrow_spacing,
                                    intron.end - arrow_width,
                                    arrow_spacing + arrow_width)
            arrow_end = arrow_start + arrow_width

            if arrow_start.size > 0:
                r.segments(x0=robjects.FloatVector(arrow_start),
                           y0=r.c(arrow_top),
                           x1=robjects.FloatVector(arrow_end),
                           y1=r.c(midpoint), col=self.color)
                r.segments(x0=robjects.FloatVector(arrow_start),
                           y0=r.c(arrow_bottom),
                           x1=robjects.FloatVector(arrow_end),
                           y1=r.c(midpoint), col=self.color)
            
        else:
            arrow_start = np.arange(intron.end - arrow_spacing - arrow_width,
                                    intron.start,
                                    -(arrow_spacing + arrow_width))
            arrow_end = arrow_start + arrow_width

            if arrow_start.size > 0:
                r.segments(x0=robjects.FloatVector(arrow_start),
                           y0=r.c(midpoint),
                           x1=robjects.FloatVector(arrow_end),
                           y1=r.c(arrow_top), col=self.color)
                r.segments(x0=robjects.FloatVector(arrow_start),
                           y0=r.c(midpoint),
                           x1=robjects.FloatVector(arrow_end),
                           y1=r.c(arrow_bottom), col=self.color)


    def draw_label(self, r):
        # draw label
        tr = self.transcript
        if tr.name:
            region_len = self.region.end - self.region.start + 1.0
            left_space = tr.start - self.left
            right_space = self.right - tr.end
            midpoint = (self.top + self.bottom)*0.5
            offset = region_len * 0.005
            
            if left_space > right_space:
                if left_space > 0:
                    # draw label to left of transcript start
                    r.text(x=tr.start-offset, y=midpoint, pos=2,
                           labels=r.c(tr.name), col=self.color,
                           cex=self.cex * 0.8)
            else:
                if right_space > 0:
                    # draw label to right of transcript end
                    r.text(x=tr.end + offset, y=midpoint, pos=4,
                           labels=r.c(tr.name), col=self.color,
                           cex=self.cex * 0.8)


    def draw_track(self, r):        
        tr = self.transcript
        
        # draw exons
        for ex in tr.exons:
            self.draw_exon(r, ex)

        # draw introns
        for intron in tr.get_introns():
            self.draw_intron(r, intron)

        if self.do_label:
            self.draw_label(r)

