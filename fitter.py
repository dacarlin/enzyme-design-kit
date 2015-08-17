# -*- coding: utf-8 -*-
# the friendly way to enzyme kinetics 

import pandas 
import os 
import random
import StringIO
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, make_response
from numpy import diag, sqrt  
from scipy.optimize import curve_fit
from werkzeug import secure_filename
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt

app = Flask( __name__ )
app.config[ 'UPLOAD_FOLDER' ] = 'uploads' 

pandas.options.display.max_colwidth = 10000 
pandas.options.display.max_rows = 10000 

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
    'plot': plt.subplots()
    #'ki': '%.4f ± %.4f' % (popt[2], perr[2]),
  } )
  
@app.route( '/', methods=['GET', 'POST'] )
def upload_file():
  if request.method == 'POST':
    file = request.files[ 'file' ]
    df = pandas.read_csv( file ) 
    grouped = df.groupby( by='sample' ).apply( fit ) 
    return '<pre>' + str( grouped ) + '</pre>'
  #return render_template( 'index.html' ) 

  return '''
    <!doctype html>
    <title>fitter.py</title>
    <style>
      div { 
        font: normal 18px 'Source Serif Pro'; 
        background-color: #eee; 
        width: 400px; 
        padding: 3em; 
        margin: 200px auto; 
      }  
      form input { font-size: 18px; } 
      #file { width: 200px; height: 200px; background-color: #222; } 
    </style>
    <div>
    <h1>Upload your data</h1>
    <form action="" method=post enctype=multipart/form-data>
      <input type=file name=file id=file>
      <input type=submit value=Fit>
    </form>
    </div>
    '''

if __name__ == '__main__':
  app.run( debug=True ) 
