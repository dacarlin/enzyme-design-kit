from flask import Flask, render_template, request, redirect, url_for, send_from_directory
from numpy import diag, sqrt  
import pandas 
from scipy.optimize import curve_fit
import os 
from werkzeug import secure_filename

UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = set( [ 'txt', 'csv', 'xlsx' ] ) # would be nice to accept xlsx 

pandas.options.display.max_colwidth = 10000 
pandas.options.display.max_rows = 10000 

app = Flask( __name__ )
app.config[ 'UPLOAD_FOLDER' ] = UPLOAD_FOLDER 

def allowed_file( filename ):
  return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def kobs( s, kcat, km ): return (kcat*s)/(km+s) 

def fit( df ):
  p0 = ( df.kobs.max(), df.s.mean() ) # guesses  
  popt, pcov = curve_fit( kobs, df.s, df.kobs, p0=p0 ) 
  perr = sqrt(diag(pcov))
  return 'kcat: %.0f +- %.0f, km: %.4f +- %.5f' % ( popt[0], perr[0], popt[1], perr[1] )  

@app.route( '/', methods=['GET', 'POST'] )
def upload_file():
  if request.method == 'POST':
    file = request.files[ 'file' ]
    df = pandas.read_csv( file ) 
    grouped = df.groupby( by='sample' ).apply( fit ) 
    return '<pre>' + str( grouped ) + '</pre>'

  return render_template( 'index.html' ) 
  
if __name__ == '__main__':
  app.run( debug=True ) 
