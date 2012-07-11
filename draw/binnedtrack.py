
import Track from track

class BinnedData(object):
    def __init__(self, values, region, bin_sz):
        vals = values[(region.start-1):region.end]
        vals = np.nan_to_num(vals)
        
        self.region = region
        self.bin_starts = np.arange(region.start, region.end, bin_sz)
        self.bin_ends = self.bin_starts + bin_sz
        self.bin_means = np.empty(self.bin_starts.size)

        # calculate mean value in each bin
        for i in range(0, self.bin_starts.size):
            start = self.bin_starts[i] - region.start
            end = self.bin_ends[i] - region.start + 1
            self.bin_means[i] = np.mean(vals[start:end])



class BinnedTrack(Track):
    
    def __init__(self, bins, region, height=1.0,
                 border=False, color="black", max_val=None, min_val=None):
        
        # call superclass constructor
        super(BinnedTrack, self).__init__(height=height, border=border)
        
        self.bins = bins
        self.color = color

        if max_val is None:
            self.max_val = np.max(bins.bin_means)
        else:
            self.max_val = max_val

        if min_val is None:
            self.min_val = np.min(bins.bin_means)
        else:
            self.min_val = np.min(bins.bin_means)

    
    def draw_track(self, r):
        yscale = self.height / (self.max_val - self.min_val)
        
        xleft = self.bins.bin_starts-1.0
        xright = self.bins.bin_ends

        n_bin = self.bins.bin_starts.size
        ybottom = np.zeros(n_bin) + self.bottom
        ytop = (self.bins.bin_means - self.min_val)*yscale + self.bottom
                
        r.rect(robjects.FloatVector(xleft),
               robjects.FloatVector(ybottom),
               robjects.FloatVector(xright),
               robjects.FloatVector(ytop), col=self.color,
               border=robjects.r("NA"))
