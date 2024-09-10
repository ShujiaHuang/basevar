"""
==============================================
Same with GATK
==============================================

Author : Shujia Huang
Date   : 2014-05-20 17:49:27
"""

class VariantDatum(object):

    def __init__ (self):
        self.annotations         = [] # Will be normalize and use for VQSR
        self.raw_annotations     = [] # Keep the raw value of each variant
        self.lod                 = None 
        self.prior               = 2.0
        self.atTrainingSite      = False
        self.atAntiTrainingSite  = False
        self.failingSTDThreshold = False
        self.worstAnnotation     = None
        self.variantOrder        = None

