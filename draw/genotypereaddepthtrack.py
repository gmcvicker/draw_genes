import sys

import numpy as np

import genome.db

import genome.trackstat

from continuoustrack import ContinuousTrack

import rpy2.robjects as robjects


SNP_UNDEF = -1



# set default colors for each genotype class
# some nice genotype colors:
# dark blue, purple, orange
#   "203569" - dark blue
#   "98509F" - purple
#   "F47927" - orange
DEFAULT_REF_COLOR = "#203569"
DEFAULT_HET_COLOR = "#98509F"
DEFAULT_ALT_COLOR = "#F47927"

# or shades of blue:
#  "#2A3A42"
#  "#6C8796"
#  "#B9E5FB"


GENO_PROB_THRESH = 0.90


class GenotypeReadDepthTrack(ContinuousTrack):     
        
    def __init__(self, region, options):
        if "individual_file" not in options:
            raise ValueError("Config for track should specify "
                             "INDIVIDUAL_FILE option")

        self.individual_file = options['individual_file']

        # read individuals from file, group by genotype
        inds_by_geno = self.get_individuals_by_geno(region, options)
        
        ref_tracks = self.get_tracks(inds_by_geno['ref'], options)
        het_tracks = self.get_tracks(inds_by_geno['het'], options)
        alt_tracks = self.get_tracks(inds_by_geno['alt'], options)

        gdb = options['gdb']

        # TODO: could allow colors to be specified in track options
        self.ref_color = DEFAULT_REF_COLOR
        self.het_color = DEFAULT_HET_COLOR
        self.alt_color = DEFAULT_ALT_COLOR

        ref_vals, ref_total_mapped = \
          self.get_vals_and_total_mapped(ref_tracks, region, options)
        het_vals, het_total_mapped = \
          self.get_vals_and_total_mapped(het_tracks, region, options)
        alt_vals, alt_total_mapped = \
          self.get_vals_and_total_mapped(alt_tracks, region, options)

        sys.stderr.write("individuals, tracks, total_mapped_reads by genotype:\n"
                         "  %d/%d/%d, %d/%d/%d, %d/%d/%d\n" %
            (len(inds_by_geno['ref']), len(inds_by_geno['het']), 
             len(inds_by_geno['alt']), len(ref_tracks), len(het_tracks), 
             len(alt_tracks), ref_total_mapped, het_total_mapped,
             alt_total_mapped))

        ref_vals = self.rescale_values(ref_vals, ref_total_mapped, options)
        het_vals = self.rescale_values(het_vals, het_total_mapped, options)
        alt_vals = self.rescale_values(alt_vals, alt_total_mapped, options)
        
        if 'smooth' in options:
            smooth_win_sz = int(options['smooth'])
        else:
            smooth_win_sz = 1
        
        if smooth_win_sz > 1:
            if 'smoother' in options:
                smoother = options['smoother']
            else:
                smoother = 'average'

            self.ref_vals = self.smooth_values(ref_vals, smooth_win_sz, 
                                               smoother)
            self.het_vals = self.smooth_values(het_vals, smooth_win_sz, 
                                               smoother)
            self.alt_vals = self.smooth_values(alt_vals, smooth_win_sz, 
                                               smoother)
        else:
            self.ref_vals = ref_vals
            self.het_vals = het_vals
            self.alt_vals = alt_vals

        # draw the reference in back, alt in front?
        if 'ref_in_back' in options:
            self.ref_in_back = self.parse_bool_str(options['ref_in_back'])
        else:
            self.ref_in_back = True
            
        self.init_attrib(region, options)
        self.set_y_range(options)

        if 'n_ticks' in options:
            self.n_ticks = int(options['n_ticks'])
        else:
            self.n_ticks = 3



        

        
                  
    def read_all_individuals(self):
        ind_list = []
        f = open(self.individual_file)
        for l in f:
            ind = l.split()[0].replace("NA", "")
            ind_list.append(ind)
        f.close()

        return ind_list


    
    def get_snp_genotypes(self, region, individuals, options):
        """Retrieves genotypes for all individuals for this 
        region's SNP"""
        gdb = options['gdb']

        if 'snp_index_track' in options:
            snp_index_trackname = options['snp_index_track']
        else:
            sys.stderr.write("no SNP_INDEX_TRACK specified in config, "
                             "assuming /impute2/snp_index\n")
            snp_index_trackname = 'impute2/snp_index'

        if 'geno_prob_track' in options:
            geno_prob_trackname = options['geno_prob_track']
        else:
            sys.stderr.write("no GENO_PROB_TRACK specified in config, "
                             "assuming /impute2/yri_geno_probs\n")
            geno_prob_trackname = 'impute2/yri_geno_probs'
            
        snp_index_track = gdb.open_track(snp_index_trackname)
        geno_track = gdb.open_track(geno_prob_trackname)
        geno_tab = geno_track.h5f.getNode("/%s" % region.chrom.name)

        snp_positions = [int(x) for x in region.snp_pos.split(",")]

        # just use first snp specified
        if len(snp_positions) < 1:
            raise ValueError("regions must specify at least one snp_pos")
        
        snp_pos = snp_positions[0]
        
        i = snp_index_track.get_val(region.chrom, snp_pos)

        if i == SNP_UNDEF:
            raise ValueError("there is no SNP at position %s:%d\n" %
                             (region.chrom.name, snp_pos))

        geno_probs = geno_tab[i,]

        return geno_probs



    def get_individuals_by_geno(self, region, options):
        """Returns a dictionary with four keys: 'ref', 'het', 'alt',
        'unk'.  The values for each key are lists of individual
        identifiers with the genotypes: homozygous reference,
        heterozygous, homozygous alternative, unknown."""
        individuals = self.read_all_individuals()

        genotypes = self.get_snp_genotypes(region, individuals, options)

        ind_dict = {'ref' : [], 'het' : [], 'alt' : [], 'unk' : []}
        
        for i in range(len(individuals)):
            # look at genotype probs for individual
            ref_prob = genotypes[i*3]
            het_prob = genotypes[i*3+1]
            alt_prob = genotypes[i*3+2]

            tot_prob = ref_prob + het_prob + alt_prob

            if tot_prob > 1.01 or tot_prob < 0.99:
                sys.stderr.write("WARNING: genotype probabilities for "
                                 "individual %s do not add to 1.0: "
                                 "%.2f %.2f %.2f\n" % (individual[i], ref_prob, 
                                                       het_prob, alt_prob))
                ind_dict['unk'].append(individuals[i])
            elif ref_prob > GENO_PROB_THRESH:
                ind_dict['ref'].append(individuals[i])
            elif het_prob > GENO_PROB_THRESH:
                ind_dict['het'].append(individuals[i])
            elif alt_prob > GENO_PROB_THRESH:
                ind_dict['alt'].append(individuals[i])
            else:
                ind_dict['unk'].append(individuals[i])
                
        return ind_dict



    def get_vals_and_total_mapped(self, tracks, region, options):
        gdb = options['gdb']

        # sum values and total mapped reads from all tracks
        values = np.zeros(region.length())        
        total_mapped_reads = 0
        
        for track in tracks:
            stat = genome.trackstat.get_stats(gdb, track)
            total_mapped_reads += stat.sum
                            
            values[:] += track.get_nparray(region.chrom,
                                           start=region.start,
                                           end=region.end)
            track.close()

        return values, total_mapped_reads


    def rescale_values(self, values, total_mapped, options):
        # rescale by total number of sequenced reads for these
        # individuals
        if total_mapped > 0:
            # convert to FPKM
            scale = 1e9 / float(total_mapped)
            values *= scale
          
            log_scale = False
            if "log_scale" in options:
                log_scale = self.parse_bool_str(options['log_scale'])

            if log_scale:
                # add one to values, but avoid possible overflow of
                # 8 bit values
                f = values < 255
                values[f] += 1
                values = np.log2(values)
        else:
            values[:] = np.nan

        return values


    def set_y_range(self, options):
        """Sets the maximum and minimum values of the
        y-axis. Overrides method from parent class to consider values
        from all three genotypes."""
        # first set maximum and minimum values to range of data

        max_list = []
        min_list = []

        # find min and max values across all three genotypes
        for v in [self.ref_vals, self.het_vals, self.alt_vals]:
            is_nan = np.isnan(v)
            if np.any(~is_nan):
                max_list.append(np.max(v[~is_nan]))
                min_list.append(np.min(v[~is_nan]))
            else:
                max_list.append(0.0)
                min_list.append(0.0)

        self.max_val = max(max_list)
        self.min_val = min(min_list)
                
        # Make sure that y-axis maximum and minimum are at least
        # those specified by soft_min_val and soft_max_val.
        # These are 'soft' limits in that the axis is allowed to
        # grow beyond them to accomodate larger/smaller values.
        if "soft_max_val" in options:
            # this is the minimum maximum y-axis value
            if float(options['soft_max_val']) > self.max_val:
                self.max_val = float(options['soft_max_val'])

        if "soft_min_val" in options:
            #  maximum allowed minimum value for y-axis
            if float(options['soft_min_val']) < self.min_val:
                self.min_val = float(options['soft_min_val'])

        # If max and min values are specified these are 'hard'
        # limits. Just set the range of the x-axis to these values
        if 'max_val' in options:
            self.max_val = float(options['max_val'])

        if 'min_val' in options:
            self.min_val = float(options['min_val'])


            

    

    def get_tracks(self, individuals, options):
        gdb = options['gdb']

        track_list = []
        for ind in individuals:
            track_name = options['track'].replace("@INDIVIDUAL@", ind)

            if gdb.has_track(track_name):
                track_list.append(gdb.open_track(track_name))
            else:
                # skip tracks where there are no matching individuals
                pass
            

        return track_list




    def draw_track(self, r):
        
        if self.max_val == self.min_val:
            yscale = self.height
        else:
            yscale = self.height / (self.max_val - self.min_val)

        # reference drawn in back or front?
        if self.ref_in_back:
            geno_colors = [self.ref_color, self.het_color, self.alt_color]
            geno_vals = [self.ref_vals, self.het_vals, self.alt_vals]
        else:
            geno_colors = [self.alt_color, self.het_color, self.ref_color]
            geno_vals = [self.alt_vals, self.het_vals, self.ref_vals]
            
        for gcol, gvals in zip(geno_colors, geno_vals):
            vals = (gvals - self.min_val) * yscale + self.bottom

            # break region into smaller blocks because
            # some programs (e.g. illustrator) don't like it when
            # polygons have too many points
            block_sz = 10000
            for block_start in range(0, vals.size, block_sz):
                block_end = min(block_start + block_sz, vals.size)

                block_vals = vals[block_start:block_end]

                # identify contiguous segments with same values
                (x1, x2, y) = self.get_segments(block_vals)

                # convert contiguous segments to polygon coordinates
                (x, y) = self.get_polygon_coords(x1, x2, y)

                if len(x) > 0:
                    r.polygon(robjects.FloatVector(x + block_start),
                              robjects.FloatVector(y),
                              col=gcol, border=gcol)

        self.draw_y_axis(r, self.n_ticks)




