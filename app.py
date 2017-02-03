import pandas
import datetime
from numpy import diag, sqrt, linspace, array, nan, empty
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import use; use( 'Agg' )
import matplotlib.pyplot as plt
from io import StringIO
from flask import Flask, request, render_template

app = Flask( __name__ )
app.debug = True

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
        perr = sqrt( diag( pcov ) ) / popt * 100
        return popt, perr
    except Exception as e:
        #print( e )
        empty3 = array( [nan, nan, nan] )
        return empty3, empty3

def do_fit( df ):
    try:
        p0 = ( df.kobs.max(), df.s.mean() )
        popt, pcov = curve_fit( kobs, df.s, df.kobs, p0=p0 )
        perr = sqrt( diag( pcov ) ) / popt * 100
        return popt, perr
    except Exception as e:
        #print( e )
        return array( [] ), array( [] )

# BglB-specific values below!
s = [ 0.075, 0.01875, 0.0047, 0.0012, 0.0003, 0.000075, 0.000019, 0 ]
extcoef = 113000

# kinetics
@app.route( '/', methods=['GET', 'POST'] )
def simple():
    if request.method == 'POST':

        # read in form
        clean_dat = request.form[ 'data' ].replace('Max V [420]', 'rate').replace(' ', '\n').lower()

        # turn into pandas df
        df = pandas.read_csv( StringIO( clean_dat ), sep='\t' )

        # map the form values to the DataFrame
        samplemap = { str(i+1): request.form[ 'mut{}-name'.format( (i//3)+1 ) ] for i in range(12) }
        yieldmap = { str(i+1): request.form[ 'mut{}-yield'.format( (i//3)+1 ) ] for i in range(12) }
        dilutionmap = { str(i+1): request.form[ 'mut{}-dilution'.format( (i//3)+1 ) ] for i in range(12) }

        # construct a useful dataframe from the raw data
        df.dropna( inplace=True )
        df = df[( df.rate > 0 )]
        df[ 's' ] = df['well'].str[0].map( dict( zip( 'abcdefgh', s ) ) )
        df[ 'sample' ] = df['well'].str[1:].map( samplemap )
        df[ 'yield' ] = df['well'].str[1:].map( yieldmap ).astype( 'float' )
        df[ 'dilution' ] = df['well'].str[1:].map( dilutionmap ).astype( 'float' )
        df[ 'kobs' ] = df.rate * 0.0002 / ( df[ 'yield' ] * df[ 'dilution' ] * 0.25 / extcoef )

        # save this data set
        #df.to_csv( '/data/bagel/uploads/submitted_{}.csv'.format( datetime.datetime.now() ) )

        # group df by sample
        grouped = df.groupby( 'sample', sort=False )
        samples = [ ]

        # iterate over 4 samples by name, in entered order (see sort=False above)
        for name, df in grouped:

            # collect metadata
            conc = '{0:.2f}'.format( df['yield'].mean() )
            dilution = df['dilution'].mean()
            popt, perr = do_fit( df )
            popt_si, perr_si = do_substrate_inhibition_fit( df )

            # set up plot
            fig, ax = plt.subplots( figsize=(3,3) )
            ax.scatter( df.s, df.kobs, color='orange' )
            ax.set_title( name )
            ax.set_xlabel( '[pNPG] (mM)' )
            ax.set_ylabel( 'Rate (min$^{-1}$)' )
            ax.set_xticks([ 0, 0.025, 0.050, 0.075 ])
            ax.set_xticklabels([0, 25, 50, 75])
            yticks = ax.get_yticks()
            ax.set_yticks( yticks[1:-1] )
            plt.tight_layout()

            xvals = linspace( df.s.min(), df.s.max(), 100 )
            # check if we have evidence of substrate inhibition
            if popt_si[2] < 0.075:
                ax.plot( xvals, kobs_with_substrate_inhibition( xvals, *popt_si ), color='g' )
            # check if we have MM
            elif popt[0] > 0.1: # dumb check
                ax.plot( xvals, kobs( xvals, *popt ), color='k' )
            # neither fit acceptable
            else:
              popt = perr = array([ nan, nan ])

            # encode plot as a string
            img = StringIO()
            fig.savefig( img, format='svg' )
            plt.close( fig )
            img.seek( 0 )

            notes = []

            # collect everything
            samples.append((
                name, conc, dilution, popt, perr, popt_si, perr_si, img.read(), notes
            ))
        return render_template( 'results.html', samples=samples )

    else: # not a POST request, is GET request
        return render_template( 'plate.html' )

# thermal stability
@app.route( '/thermal', methods=['GET', 'POST'] )
def thermal():
    if request.method == 'POST':

        clean_dat = request.form[ 'data' ].replace('Max V [420]', 'rate')

        # turn into pandas df
        df = pandas.read_csv( StringIO( clean_dat ), sep='\t' )

        # map the form values to the DataFrame
        samplemap = { str(i+1): request.form[ 'mut{}-name'.format( (i//3)+1 ) ] for i in range(12) }

        # construct a useful dataframe from the raw data
        df[ 'sample' ] = df['Well'].str[1:].map( samplemap )
        df[ 'temp' ] = list( np.linspace( 50, 30, 8) ) * 12
        #df.to_csv( '/data/bagel/uploads/thermal_submitted_{}.csv'.format( datetime.datetime.now() ) )

        # group df by sample
        grouped = df.groupby( 'sample', sort=False )
        result = []
        #iterate over 4 samples by name, in entered order (see sort=False above)
        for name, df in grouped:
            rate = df.rate / df.rate.max()
            params, errors = my_curve_fit( r, df.temp, rate, p0=(40,-1) )
            fig, ax = plt.subplots( figsize=(3,3) )
            ax.scatter( df.temp, rate, color='red' )
            x_vals = np.linspace( 30, 50, 50 )
            ax.plot( x_vals, r(x_vals, *params), color='black' )
            ax.set_xticks([30, 35, 40, 45, 50])
            ax.set_yticks([0, .25, .5, .75, 1])
            ax.set_xlabel( 'Incubation temp. (˚C)' )
            ax.set_ylabel( 'Normalized activity' )
            fig.tight_layout()

            # render into HTML
            img = StringIO()
            fig.savefig( img, format='svg' )
            plt.close( fig )
            img.seek( 0 )
            result.append( ( name, params, errors, img.read() ) )

        return render_template( 'thermal_results.html', result=result )
    else:
        return render_template( 'thermal.html' )
