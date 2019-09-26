"""
Author : Shujia Huang
Date   : 2014-05-20 17:49:27
"""


class VariantDatum(object):

    def __init__(self):
        self.annotations = []  # Will be normalize and use for VQSR
        self.raw_annotations = {}  # Keep the raw value of each variant
        self.lod = None
        self.prior = 2.0
        self.at_training_site = False
        self.at_anti_training_site = False
        self.failing_STD_threshold = False
        self.worst_annotation = None
        self.variant_order = None
