import pandas
from numpy import diag, sqrt
from scipy.optimize import curve_fit

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from StringIO import StringIO
from base64 import b64encode

from flask import Flask, request, render_template

app = Flask( __name__ )
app.config[ 'UPLOAD_FOLDER' ] = 'uploads'

def allowed_file( filename ):
  return '.' in filename and filename.rsplit('.', 1)[1] in [ 'txt', 'csv', 'xlsx' ]

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
  # 0.25 from the procedure, 0.0002 from standard curve
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

@app.route( '/simple', methods=['GET', 'POST'] )
def simple():
  if request.method == 'POST':
    clean_dat = request.form.get( 'data' ).replace('Max V [420]', 'rate').replace(' ', '\n').lower()
    df = pandas.read_csv( StringIO( clean_dat ), sep='\t' )

    samplemap = {
        '1': request.form.get( 'mut1-name' ),
        '2': request.form.get( 'mut1-name' ),
        '3': request.form.get( 'mut1-name' ),
        '4': request.form.get( 'mut2-name' ),
        '5': request.form.get( 'mut2-name' ),
        '6': request.form.get( 'mut2-name' ),
        '7': request.form.get( 'mut3-name' ),
        '8': request.form.get( 'mut3-name' ),
        '9': request.form.get( 'mut3-name' ),
        '10': request.form.get( 'mut4-name' ),
        '11': request.form.get( 'mut4-name' ),
        '12': request.form.get( 'mut4-name' ),
    }

    yieldmap = {
        '1': request.form.get( 'mut1-yield' ),
        '2': request.form.get( 'mut1-yield' ),
        '3': request.form.get( 'mut1-yield' ),
        '4': request.form.get( 'mut2-yield' ),
        '5': request.form.get( 'mut2-yield' ),
        '6': request.form.get( 'mut2-yield' ),
        '7': request.form.get( 'mut3-yield' ),
        '8': request.form.get( 'mut3-yield' ),
        '9': request.form.get( 'mut3-yield' ),
        '10': request.form.get( 'mut4-yield' ),
        '11': request.form.get( 'mut4-yield' ),
        '12': request.form.get( 'mut4-yield' ),
    }

    dilutionmap = {
        '1': request.form.get( 'mut1-dilution' ),
        '2': request.form.get( 'mut1-dilution' ),
        '3': request.form.get( 'mut1-dilution' ),
        '4': request.form.get( 'mut2-dilution' ),
        '5': request.form.get( 'mut2-dilution' ),
        '6': request.form.get( 'mut2-dilution' ),
        '7': request.form.get( 'mut3-dilution' ),
        '8': request.form.get( 'mut3-dilution' ),
        '9': request.form.get( 'mut3-dilution' ),
        '10': request.form.get( 'mut4-dilution' ),
        '11': request.form.get( 'mut4-dilution' ),
        '12': request.form.get( 'mut4-dilution' ),
    }

    df[ 's' ] = df['well'].str[0].map( dict( zip( 'abcdefgh', s ) ) )
    df[ 'sample' ] = df['well'].str[1:].map( samplemap )
    df[ 'yield' ] = df['well'].str[1:].map( yieldmap ).astype( 'float' )
    df[ 'dilution' ] = df['well'].str[1:].map( dilutionmap ).astype( 'float' )
    df[ 'kobs' ] = df.rate * 0.0002 / ( df[ 'yield' ] * df[ 'dilution' ] * 0.25 / extcoef )

    print df

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
    return render_template( 'simple.html' )

if __name__ == '__main__':
  app.run( debug=True )
