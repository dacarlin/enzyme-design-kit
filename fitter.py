import pandas
from numpy import diag, sqrt, linspace
from scipy.optimize import curve_fit

from matplotlib import use; use( 'Agg' )
import matplotlib.pyplot as plt
from StringIO import StringIO
from base64 import b64encode

from flask import Flask, request, render_template

app = Flask( __name__ )
app.config[ 'UPLOAD_FOLDER' ] = 'uploads'
app.debug = True

def allowed_file( filename ):
  return '.' in filename and filename.rsplit('.', 1)[1] in [ 'csv' ]

def kobs( s, kcat, km ):
  return (kcat*s)/(km+s)

def do_fit( df ):
  p0 = ( df.kobs.max(), df.s.mean() )
  popt, pcov = curve_fit( kobs, df.s, df.kobs, p0=p0 )
  perr = sqrt( diag( pcov ) )
  for i in [ 0, 1 ]:
    if not popt[ i ] or perr[ i ] > popt[ i ]:
      popt[ i ] = perr[ i ] = None
    else:
      perr[ i ] = perr[ i ] / popt [ i ] * 100
  return popt, perr

# BglB-specific values below!
s = [ 0.075, 0.01875, 0.0047, 0.0012, 0.0003, 0.000075, 0.000019, 0 ]
smap = dict( zip( 'ABCDEFGH', s ) )
extcoef = 113000

def assign_groups( df ):
  df[ 's' ] = df['well'].str[0].map( smap )
  df[ 'kobs' ] = df.rate * 0.0002 / ( df[ 'yield' ] * df[ 'dilution' ] * 0.25 / extcoef )
  return df

@app.route( '/' )
def index():
  return render_template( 'help.html' )

@app.route( '/batch', methods=['GET', 'POST'] )
def batch():
  if request.method == 'POST':
    df = pandas.read_csv( request.files[ 'file' ] )
    df = assign_groups( df )
    grouped = df.groupby( 'sample' )
    plots = [ ]
    for name, df in grouped:
      conc = '{0:.2f}'.format( df['yield'].mean() )
      popt, perr = do_fit( df )
      fig, ax = plt.subplots()
      fig.suptitle( name )
      ax = plt.scatter( df.s, df.kobs )
      img = StringIO()
      fig.savefig( img )
      plt.close( fig )
      img.seek( 0 )
      plots.append( ( name, conc, popt, perr, b64encode( img.read() ) ) )
    return render_template( 'results.html', plots=plots )
  else:
    return render_template( 'batch.html' )

@app.route( '/plate', methods=['GET', 'POST'] )
def simple():
  if request.method == 'POST':
    clean_dat = request.form.get( 'data' ).replace('Max V [420]', 'rate').replace(' ', '\n').lower()
    df = pandas.read_csv( StringIO( clean_dat ), sep='\t' )

    samplemap = { str(i+1): request.form.get( 'mut{}-name'.format( (i/3)+1 ) ) for i in range(12) }
    yieldmap = { str(i+1): request.form.get( 'mut{}-yield'.format( (i/3)+1 ) ) for i in range(12) }
    dilutionmap = { str(i+1): request.form.get( 'mut{}-dilution'.format( (i/3)+1 ) ) for i in range(12) }

    df[ 's' ] = df['well'].str[0].map( dict( zip( 'abcdefgh', s ) ) )
    df[ 'sample' ] = df['well'].str[1:].map( samplemap )
    df[ 'yield' ] = df['well'].str[1:].map( yieldmap ).astype( 'float' )
    df[ 'dilution' ] = df['well'].str[1:].map( dilutionmap ).astype( 'float' )
    df[ 'kobs' ] = df.rate * 0.0002 / ( df[ 'yield' ] * df[ 'dilution' ] * 0.25 / extcoef )

    grouped = df.groupby( 'sample' )
    plots = [ ]

    from matplotlib.ticker import FuncFormatter
    def substrateFormatter(x, pos):
        if x < 0 or x > 0.090:
            return ''
        else:
            return int( x * 1000 )

    def rateFormatter(x, pos):
        if x < 0:
            return ''
        elif 0 < x < 5:
            return '{:0.1f}'.format( x )
        else:
            return int( x )

    for name, df in grouped:
      conc = '{0:.2f}'.format( df['yield'].mean() )
      popt, perr = do_fit( df )
      fig, ax = plt.subplots( figsize=(4,4) )
      ax.set_title( name )
      ax.scatter( df.s, df.kobs, color='burlywood', alpha=0.99 )
      xvals = linspace( df.s.min(), df.s.max(), 100 )
      ax.plot( xvals, kobs( xvals, *popt ), alpha=0.45, color='#222222' )
      ax.xaxis.set_major_formatter(FuncFormatter(substrateFormatter))
      ax.yaxis.set_major_formatter(FuncFormatter(rateFormatter))
      ax.set_xlabel( 'Substrate concentration (mM)' )
      ax.set_ylabel( 'Rate (min$^{-1}$)' )
      plt.tight_layout()
      img = StringIO()
      fig.savefig( img )
      plt.close( fig )
      img.seek( 0 )
      plots.append( ( name, conc, popt, perr, b64encode( img.read() ) ) )

    return render_template( 'simple_results.html', plots=plots )

  else:
    return render_template( 'simple.html' )

if __name__ == '__main__':
  app.run()
