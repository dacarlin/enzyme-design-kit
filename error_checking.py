from numpy import zeros 
from scipy.optimize import curve_fit 

def fit( x, y ):

  '''
  Input: 

    x: array 

    y: array, same length as x 

  Output: 
  
    kcat: 
    
    km:  

  '''

  def v( x, kcat, km ):
    return (kcat*s)/(km+s) 

  p0 = ( 1e3, 5e-3 ) 
  popt, perr = curve_fit( v, x, y, p0 ) or 2 * zeros( 2 ) 
  
  checks = [
    ( kcat < 0.1 or kcat > 1e5, 'kcat out of range' ), 
    ( km < 5e-5 or km > 200, 'km out of range' ), 
    ( kcat_err > kcat, 'kcat error > kcat' ), 
    ( km_err > km, 'km error > km' ), 
  ]


