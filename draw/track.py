import sys

import rpy2.robjects as robjects


class RowElement(object):
    def __init__(self, feature, prev=None, next=None, padding=0.0):
        self.feature = feature
        self.start = feature.start - padding
        self.end = feature.end + padding
        self.prev = prev
        self.next = next


class Row(object):
    def __init__(self, id):
        self.id = id
        self.elements = []
        self.first = None
        self.last = None


    def add_element(self, new_elem):
        """adds a feature to this row, returns True on success
        or False if there was no room for the feature"""

        if self.first is None:
            # create a new element
            #sys.stderr.write("added feat %d-%d as only element of row %d\n" %
            #                 (new_elem.start, new_elem.end, self.id))
            new_elem.prev = None
            new_elem.next = None
            self.first = self.last = new_elem
            return True

        cur = self.first
        while(cur and cur.start < new_elem.start):
            cur = cur.next

        if cur is None:
            if self.last.end < new_elem.start:
                # add feature to end of list
                # sys.stderr.write("added feat %d-%d to end of row %d\n" %
                #                 (new_elem.start, new_elem.end, self.id))
                new_elem.prev = self.last
                new_elem.next = None
                self.last.next = new_elem
                self.last = new_elem
                return True
            else:
                # last element overlaps with this one
                return False

        prev = cur.prev
        if prev is None:
            if new_elem.end < cur.start:
                # can insert feature at start of list
                # sys.stderr.write("added feat %d-%d to start of row %d "
                #                 "(before %d-%d)\n" %
                #                 (new_elem.start, new_elem.end, self.id,
                #                 cur.start, cur.end))
                new_elem.prev = None
                new_elem.next = cur
                cur.prev = new_elem
                self.first = new_elem
                return True
            else:
                # first feature overlaps with this one
                return False
        
        if prev.end < new_elem.start and cur.start > new_elem.end:
            # can insert feature between these two elements
            #sys.stderr.write("added feat %d-%d between %d-%d and %d-%d "
            #                 "in row %d\n" %
            #                 (new_elem.start, new_elem.end,
            #                  prev.start, prev.end,
            #                  cur.start, cur.end, self.id))
            
            new_elem.prev = prev
            new_elem.next = cur
            prev.next = new_elem
            cur.prev = new_elem
            return True

        # could not find place to insert new element in this row
        return False
            
        
        
        
        
        
class Track(object):
    """An abstract base class for all genome tracks that can be added
    to a Window and drawn. Subclasses should provide an implementation
    of the draw method."""
    
    def __init__(self, region, options):
        self.init_attrib(region, options)



    def init_attrib(self, region, options):
        if 'height' in options:
            self.height = float(options['height'])
        else:
            self.height = 1.0

        if 'cex' in options:
            self.cex = float(options['cex'])
        else:
            self.cex = 1.0

        if 'track_label' in options:
            self.track_label = options['track_label'].rstrip()
        else:
            self.track_label = ""

        self.set_colors(options)

        self.region = region
        self.left = None
        self.right = None
        self.top = None
        self.bottom = None

        self.n_fwd_row = 0
        self.n_rev_row = 0
        self.n_row = 0
        self.row_assignment = {}
        
        

    def set_colors(self, options):
        """sets the color and border attributes of this track,
        using options in provided dictionary"""
        
        # set the color and border color for the track
        if 'color' in options:
            self.color = options['color'].replace('"', '')
        else:
            # use a default color
            self.color = "#E41A1C"

        # some tracks allow different fill for fwd / rev strand
        if 'fwd_color' in options:
            self.fwd_color = options['fwd_color'].replace('"', '')
        else:
            # use default color for fwd strand
            self.fwd_color = self.color

        if 'rev_color' in options:
            self.rev_color = options['rev_color'].replace('"', '')
        else:
            self.rev_color = self.color
        

        # set colors for border
        if 'draw_border' in options:
            self.draw_border = self.parse_bool_str(options['draw_border'])
        else:
            self.draw_border = True
            
        if self.draw_border:            
            if 'border_color' in options:
                self.border_color = options['border_color'].replace('"', '')
            else:
                # use same color as fill
                self.border_color = self.color

            if 'fwd_border_color' in options:
                self.fwd_border_color = \
                    options['fwd_border_color'].replace('"', '')
            else:
                self.fwd_border_color = self.fwd_color

            if 'rev_border_color' in options:
                self.rev_border_color = \
                    options['rev_border_color'].replace('"', '')
            else:
                self.rev_border_color = self.rev_color
        else:
            # no not draw any border
            self.border_color = robjects.r("NA")
            self.rev_border_color = robjects.r("NA")
            self.fwd_border_color = robjects.r("NA")
                

        

    def parse_bool_str(self, bool_str):
        """Utility method for parsing boolean option strings"""
        if bool_str.lower() in ('on', 'yes', 'true', '1'):
            return True
        elif bool_str.lower() in ('off', 'no', 'false', '0'):
            return False
        else:
            raise ValueError("expected boolean string to be one of: %s",
                             ", ".join(['on', 'yes', 'true', '1',
                                        'off', 'no', 'false', '0']))
        return False

    def set_position(self, left, right, top, bottom):
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom


    def draw_track_label(self, r, color="black"):
        """draws a label on the left side of the track"""
        
        y_mid = self.top - (self.top - self.bottom) / 2

        # start text coordinate 1% of the way into the track
        x = self.left + (self.right - self.left + 1) * 0.01
        
        r.text(x=x, y=y_mid, pos=4, labels=self.track_label,
               cex=self.cex, col=color)



    def draw_track(self, r):
        # this is where the main work in drawing should be implemented
        # by child classes
        sys.stderr.write("WARNING: draw_track function not implemented\n")
        pass


    def draw(self, r):
        self.draw_track(r)

        if self.track_label:
            self.draw_track_label(r)


    def assign_feature_rows(self, features, use_strands=True,
                            padding=0.0):
        fwd_rows = []
        rev_rows = []

        self.n_row = 0
        self.n_rev_row = 0
        self.n_fwd_row = 0

        row_ids = {}

        for feat in features:
            row_elem = RowElement(feat, padding=padding)
            
            if use_strands and feat.strand == -1:
                placed_elem = False
                for row in rev_rows:
                    if(row.add_element(row_elem)):
                        # managed to fit feature into this row
                        placed_elem = True
                        row_ids[feat] = row.id
                        break

                if not placed_elem:
                    # create a new row
                    self.n_rev_row += 1
                    row = Row(id = -self.n_rev_row)
                    rev_rows.append(row)
                    row_ids[feat] = row.id
                    if not row.add_element(row_elem):
                        raise Exception("could not place feature in new row")

            else:
                placed_elem = False
                for row in fwd_rows:
                    if(row.add_element(row_elem)):
                        # managed to fit feature into this row
                        placed_elem = True
                        row_ids[feat] = row.id
                        break

                if not placed_elem:
                    # create a new row
                    self.n_fwd_row += 1
                    row = Row(id = self.n_fwd_row)
                    fwd_rows.append(row)
                    row_ids[feat] = row.id
                    if not row.add_element(row_elem):
                        raise Exception("could not place feature in new row")

        self.n_row = self.n_fwd_row + self.n_rev_row

        if use_strands:
            # add extra row for separation between fwd/rev strands
            self.n_row += 1

        self.row_assignment = {}
        for feat, row_id in row_ids.items():
            # update assignment of reverse-strand features so that
            # all numbers are positive (they just come after
            # fwd strand rows)
            if row_id < 0:
                self.row_assignment[feat] = -row_id + self.n_fwd_row
            else:
                self.row_assignment[feat] = row_id - 1

            
                    
        
