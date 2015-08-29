# -*- coding: utf-8 -*-
# the friendly way to enzyme kinetics 

import pandas 
from numpy import diag, sqrt  
from scipy.optimize import curve_fit

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt

from flask import Flask, request 

app = Flask( __name__ )
app.config[ 'UPLOAD_FOLDER' ] = 'uploads' 

def allowed_file( filename ):
  return '.' in filename and filename.rsplit('.', 1)[1] in [ 'txt', 'csv', 'xlsx' ]

def kobs( s, kcat, km ): 
  return (kcat*s)/(km+s)
  #return (kcat*s)/(km+s+(s**2/ki)) 
  #Vmax/E = ( kcat * S ) / ( Km + S + ( S^2 / Ki ) ) 

def fit( df ):
  p0 = ( df.kobs.max(), df.s.mean() ) #df.s.mean() ) # guesses  
  popt, pcov = curve_fit( kobs, df.s, df.kobs, p0=p0 ) 
  perr = sqrt( diag( pcov ) )
  return pandas.Series( { 
    'kcat': '%.0f ± %.0f' % (popt[0], perr[0]), 
    'km': '%.4f ± %.4f' % (popt[1], perr[1]),
    #'ki': '%.4f ± %.4f' % (popt[2], perr[2]),
  } )

# BglB-specific values below! 
smap = { 'A': 0.75, 'B': 0.1875, 'C': 0.047, 'D': 0.012, 'E': 0.003, 'F': 0.00075, 'G': 0.00019, 'H': 0 } 
extcoef = 113000

def assign_groups( df, smap=smap, extcoef=extcoef ):
  'Adds substrate concentrations and calculates kobs'

  df[ 's' ] = df['well'].str[0].map( smap ) 
  df[ 'kobs' ] = df.rate * 0.0002 / ( df[ 'yield' ] * df[ 'dilution' ] * 0.25 / extcoef ) 
  # 0.25 from the procedure, 0.0002 from standard curve 
  
  return df 
  
@app.route( '/', methods=['GET', 'POST'] )
def upload_file():
  if request.method == 'POST':
    df = pandas.read_csv( request.files[ 'file' ] ) 
    df = assign_groups( df ) 
    grouped = df.groupby( by='sample' ).apply( fit ) 
    return grouped.to_html() 
  else: #GET request
    return '''
      <!doctype html><title>Michaelis-Menten fitter</title><h1>Upload your data</h1>
      <form action="" method=post enctype=multipart/form-data>
        <input type=file name=file id=file>
        <input type=submit value=Fit>
      </form>
      '''

if __name__ == '__main__':
  app.run( debug=True ) 
