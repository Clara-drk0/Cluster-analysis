"""
Horrifying function that sum the signal of an area from a rough CutOut2D object

Libraries: from astropy.nddata import Cutout2D
           from astropy.wcs import WCS 

Inputs: cutout -> CutOut2D object of the area

Outputs: signal_sum -> float. Sum in count/sec? 
        
To improve: -> Is this function really the most efficient way to estimate the signal?
            -> Once the ellipse_fitting is ready, pass the masks and sum only over the masks
"""
def signal(cutout):
    im = cutout.data
    wcs = cutout.wcs
    
    signal_sum = im.sum() 
        
    return (signal_sum)