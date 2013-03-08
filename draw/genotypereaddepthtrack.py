import sys

import numpy as np

import genome.db

import genome.trackstat

from continuoustrack import ContinuousTrack

SNP_UNDEF = -1

INDIVIDUAL_FILE = "/data/share/10_IND/IMPUTE/samples.txt"


class GenotypeReadDepthTrack(ContinuousTrack):     
        
    def __init__(self, region, options):
        inds = self.get_individuals(region, options)
        
        tracks = self.get_tracks(inds, options)

        # sum values and total mapped reads from all tracks
        values = np.zeros(region.length())
        total_mapped_reads = 0

        gdb = options['gdb']

        if 'color' in options:
            pass
        else:
            geno_class = options['genotype']

            # set default colors for each genotype class
            if geno_class == "ref":
                options['color'] = "#2A3A42"
            elif geno_class == "het":
                options['color'] = "#6C8796"
            elif geno_class == "alt":
                options['color'] = "#B9E5FB"

        for track in tracks:
            stat = genome.trackstat.get_stats(gdb, track)
            total_mapped_reads += stat.sum
                            
            values[:] += track.get_nparray(region.chrom,
                                           start=region.start,
                                           end=region.end)
            track.close()
                
        # rescale by total number of sequenced reads for these
        # individuals
        if total_mapped_reads > 0:
            scale = 1e9 / float(total_mapped_reads)
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
        
        # do the actual drawing...
        super_init = super(GenotypeReadDepthTrack, self).__init__          
        super_init(values, region, options)

                  
    def read_all_individuals(self):
        ind_list = []
        f = open(INDIVIDUAL_FILE)
        for l in f:
            ind = l.split()[0].replace("NA", "")
            ind_list.append(ind)
        f.close()

        return ind_list


    
    def get_snp_genotypes(self, region, individuals, options):
        """Retrieves genotypes for all individuals for this 
        region's SNP"""
        gdb = options['gdb']

        snp_index_track = gdb.open_track('impute2/snp_index')
        geno_track = gdb.open_track("impute2/all_geno_probs")
        geno_tab = geno_track.h5f.getNode("/%s" % region.chrom.name)

        snp_positions = [int(x) for x in region.snp_pos.split(",")]
        chisq_stat = [float(x) for x in region.chisq_stat.split(",")]

        # use SNP with highest chisq stat to determine genotype classes
        snp_pos = None
        best_stat = None
        for pos, stat in zip(snp_positions, chisq_stat):
            if best_stat is None or stat > best_stat:
                best_stat = stat
                snp_pos = pos
        
        i = snp_index_track.get_val(region.chrom, snp_pos)

        if i == SNP_UNDEF:
            raise ValueError("there is no SNP at position %s:%d\n" %
                             region.chrom.name, snp_pos)

        geno_probs = geno_tab[i,]

        return geno_probs



    def get_individuals(self, region, options):
        """Returns a list of individuals who match the genotype
        class specified in the options"""
        individuals = self.read_all_individuals()

        genotypes = self.get_snp_genotypes(region, individuals,
                                           options)

        geno_class = options['genotype']

        if geno_class == "ref":
            geno_idx = 0
        elif geno_class == "het":
            geno_idx = 1
        elif geno_class == "alt":
            geno_idx = 2
        else:
            raise ValueError("expected genotype to be one of "
                             "'ref', 'het', 'alt', but got '%s'" %
                             geno_class)

        matching_inds = []

        for i in range(len(individuals)):
            geno_prob = genotypes[i*3 + geno_idx]

            if geno_prob > 0.95:
                # this individual's genotype is the requested one
                matching_inds.append(individuals[i])

        return matching_inds


    def get_tracks(self, individuals, options):
        gdb = options['gdb']

        track_list = []
        for ind in individuals:
            track_name = options['track'].replace("@INDIVIDUAL@", ind)

            if gdb.has_track(track_name):
                track_list.append(gdb.open_track(track_name))

        return track_list









