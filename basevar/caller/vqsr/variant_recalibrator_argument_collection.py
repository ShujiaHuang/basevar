"""
====================================
====================================
Author : Shujia Huang
Date   : 2014-05-21 18:03:28

"""
class VariantRecalibratorArgumentCollection(object):

    def __init__ (self):
        self.NITER = 150
        self.NINIT = 100
        self.STD_THRESHOLD         = 8.0
        self.MIN_NUM_BAD_VARIANTS  = 1000
        self.MAX_NUM_TRAINING_DATA = 50000
        self.TRAIN_SIZE_RATE       = 0.6 # The ratio of Traing date set 
        self.CV_SIZE_RATE          = 0.2 # The ratio of cross validation data set 
        self.TEST_SIZE_RATE        = 0.2 # The ratio of test data set
        self.MAX_GAUSSIANS         = 6

        # The threshold that the positive trainingset -> negative
        self.POSITIVE_TO_NEGATIVE_RATE = 0.05
        self.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL = 6

