import pandas
import datetime

from numpy import diag, sqrt, linspace
from scipy.optimize import curve_fit

from matplotlib import use; use( 'Agg' )
import matplotlib.pyplot as plt

from StringIO import StringIO

from flask import Flask, request, render_template

app = Flask( __name__ )
app.debug = True 

def kobs( s, kcat, km ):
  return (kcat*s)/(km+s)

def do_fit( df ):
  try:
    p0 = ( df.kobs.max(), df.s.mean() )
    popt, pcov = curve_fit( kobs, df.s, df.kobs, p0=p0 )
    perr = sqrt( diag( pcov ) ) / popt * 100
    return popt, perr
  except Exception as e:
    print e
    return [], []

# BglB-specific values below!
s = [ 0.075, 0.01875, 0.0047, 0.0012, 0.0003, 0.000075, 0.000019, 0 ]
smap = dict( zip( 'ABCDEFGH', s ) )
extcoef = 113000

# URLs 
@app.route( '/', methods=['GET', 'POST'] )
def simple():
  if request.method == 'POST':
    # read in form
    clean_dat = request.form.get( 'data' ).replace('Max V [420]', 'rate').replace(' ', '\n').lower()

    # turn into pandas df
    df = pandas.read_csv( StringIO( clean_dat ), sep='\t' )

    # map the form values to the DataFrame 
    samplemap = { str(i+1): request.form.get( 'mut{}-name'.format( (i/3)+1 ) ) for i in range(12) }
    yieldmap = { str(i+1): request.form.get( 'mut{}-yield'.format( (i/3)+1 ) ) for i in range(12) }
    dilutionmap = { str(i+1): request.form.get( 'mut{}-dilution'.format( (i/3)+1 ) ) for i in range(12) }

    # construct a useful dataframe from the raw data
    df.dropna( inplace=True )
    df = df[( df.rate > 0 )]
    df[ 's' ] = df['well'].str[0].map( dict( zip( 'abcdefgh', s ) ) )
    df[ 'sample' ] = df['well'].str[1:].map( samplemap )
    df[ 'yield' ] = df['well'].str[1:].map( yieldmap ).astype( 'float' )
    df[ 'dilution' ] = df['well'].str[1:].map( dilutionmap ).astype( 'float' )
    df[ 'kobs' ] = df.rate * 0.0002 / ( df[ 'yield' ] * df[ 'dilution' ] * 0.25 / extcoef )

    # save this data set
    df.to_csv( 'saved_runs/submitted_{}.csv'.format( datetime.datetime.now() ) )

    # group df by sample
    grouped = df.groupby( 'sample', sort=False )
    samples = [ ]

    # iterate over 4 samples by name, in entered order (see sort=False above)
    for name, df in grouped:

      conc = '{0:.2f}'.format( df['yield'].mean() )
      popt, perr = do_fit( df )

      # make notes of possible errors 
      notes = []
      if float( conc ) < 0.2:
          notes.append( 'Protein yield is below 0.2 mg/mL' )

      if perr[0] > 25:
          notes.append( '<em>k</em><sub>cat</sub> error is greater than 25%' )
      if perr[1] > 25:
          notes.append( 'K<sub>M</sub> error is greater than 25%' )

      # set up plot 
      fig, ax = plt.subplots( figsize=(3,3) )
      ax.scatter( df.s, df.kobs, color='cornflowerblue', alpha=0.9, marker='.' )
      xvals = linspace( df.s.min(), df.s.max(), 100 )

      # check for MM params 
      if popt[0] > 1:
        ax.plot( xvals, kobs( xvals, *popt ), alpha=0.7, color='k' )

      # finish up plots 
      ax.set_title( name )
      ax.set_xlabel( '[pNPG] (M)' )
      ax.set_ylabel( 'Rate (min$^{-1}$)' )
      ax.set_xticks( [ 0, 0.04, 0.08 ] )
      yticks = ax.get_yticks()
      ax.set_yticks( yticks[1:-1] )
      plt.tight_layout()

      # encode plot as a string
      img = StringIO()
      fig.savefig( img, format='svg' )
      plt.close( fig )
      img.seek( 0 )

      # collect everything
      samples.append( ( name, conc, popt, perr, img.read(), notes ) )

    return render_template( 'results.html', samples=samples )

  else:
    # not a POST request, is GET request
    return render_template( 'simple.html' )

if __name__ == '__main__':
  app.run()
