import sys

import numpy as np
import scipy

import genome.trackstat
import genome.track

from .continuoustrack import ContinuousTrack

class ReadDepthTrack(ContinuousTrack):     
    def __init__(self, region, options):
        track_name = options['track']
        
        track = genome.track.Track(track_name)
        
        values = track.get_nparray(region.chrom, 
                                   start=region.start,
                                   end=region.end)


        if "scale_factor" in options or "downsample" in options:
            total_reads = self.get_total_reads(gdb, track)
        else:
            total_reads = None

        if total_reads and ("downsample" in options):
            desired_total = int(options['downsample'])

            if desired_total >= total_reads:
                sys.stderr.write("Not downsampling: total "
                                 "reads (%d) > desired reads"
                                 "(%d)\n" % (total_reads, 
                                             desired_total))
            else:                
                # perform in-place downsampling of reads
                self.downsample_reads(values, total_reads, 
                                      desired_total)
                total_reads = desired_total       

        if total_reads and ("scale_factor" in options):
            scale_factor = float(options['scale_factor'])
            scale = scale_factor / float(total_reads)
            values = values * scale
            sys.stderr.write("  total reads %d, using "
                             "scale %.3f\n" % (total_reads, scale))


        track.close()
          
        log_scale = False
        if "log_scale" in options:
            log_scale = self.parse_bool_str(options['log_scale'])

        if log_scale:
            # add one to values, but avoid possible overflow of
            # 8 bit values
            f = values < 255
            values[f] += 1
            values = np.log2(values)

        super_init = super(ReadDepthTrack, self).__init__
        super_init(values, region, options)


    def get_total_reads(self, gdb, track):
        try:
            stat = gdb.get_track_stat(track)
            return stat.sum
        except ValueError as err:
            sys.stderr.write("  WARNING: cannot scale or resample "
                             "values for "
                             "track %s because track stats have "
                             "not been calculated. run "
                             "set_track_stats.py first.\n" % track.name)

        return None

    
    def downsample_reads(self, read_counts, total_reads, desired_total):
        nonzero = np.where(read_counts > 0)[0]

        p = float(desired_total) / float(total_reads)        
        downsamp = scipy.random.binomial(read_counts[nonzero], p)
        sys.stderr.write("  total: %d, desired_total: %d -- "
                         "downsampled %d reads to %d\n" % 
                         (total_reads, desired_total, 
                          np.sum(read_counts[nonzero]), 
                          np.sum(downsamp)))
        
        read_counts[nonzero] = downsamp.astype(read_counts.dtype)
        
        
        
        
          

