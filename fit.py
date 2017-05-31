import numpy as np
from scipy.optimize import curve_fit

def kobs( s, kcat, km ):
    return (kcat*s)/(km+s)

def kobs_with_substrate_inhibition( s, kcat, km, ks ):
    return (kcat*s)/(km+s*(1+(s/ks)))

def r( x, x0, k ):
    return 1 / ( 1 + np.exp( -k * ( x - x0 ) ) )

def my_curve_fit( f, xdata, ydata, p0 ):
    null_return = [ [ np.nan ] * len(p0) ] * 2
    try:
        curve_return = curve_fit( f, xdata, ydata, p0 )
        if len( curve_return ) == 2:
            errors = np.sqrt( np.diag( curve_return[1] ) )
            return curve_return[0], errors
        else:
            return null_return
    except:
        return null_return

def do_substrate_inhibition_fit( df ):
    try:
        p0 = ( df.kobs.max(), df.s.mean(), 0.04 )
        popt, pcov = curve_fit( kobs_with_substrate_inhibition, df.s, df.kobs, p0=p0 )
        perr = np.sqrt( np.diag( pcov ) ) / popt * 100
        return popt, perr
    except Exception as e:
        #print( e )
        empty3 = np.array( [np.nan, np.nan, np.nan] )
        return empty3, empty3

def do_fit( df ):
    try:
        p0 = ( df.kobs.max(), df.s.mean() )
        popt, pcov = curve_fit( kobs, df.s, df.kobs, p0=p0 )
        perr = np.sqrt( np.diag( pcov ) ) / popt * 100
        return popt, perr
    except Exception as e:
        #print( e )
        return np.array( [] ), np.array( [] )
