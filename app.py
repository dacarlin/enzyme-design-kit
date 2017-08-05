# -*- coding: utf-8 -*-

from flask import Flask, request, render_template
import pandas
import datetime
import numpy as np
from io import StringIO

from Bio import SeqIO
from skbio import DNA
from glob import glob

from matplotlib import use; use( 'Agg' )
import matplotlib.pyplot as plt

from fit import kobs, kobs_with_substrate_inhibition, r, my_curve_fit, do_substrate_inhibition_fit, do_fit

app = Flask( __name__ )
app.debug = True

@app.route('/')
def home():
  return render_template('home.htm')

# BglB-specific values below!
s = [ 0.075, 0.01875, 0.0047, 0.0012, 0.0003, 0.000075, 0.000019, 0 ]
extcoef = 113000

def get_data_sets():
    data_sets = []
    for g in glob('/static/data_sets/*.csv'):
        with open(g) as fn:
            data_sets.append(g, fn.read())
    print(data_sets)
    return dict(data_sets)

# kinetics
@app.route( '/kinetics', methods=['GET', 'POST'] )
def kinetics():
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
        # the 0.0002 is an empirically-determined scale factor for our plate reader

        # save this data set
        #df.to_csv( '/data/bagel/uploads/submitted_{}.csv'.format( datetime.datetime.now() ) )

        # group df by sample
        grouped = df.groupby( 'sample', sort=False )
        samples = [ ]

        # iterate over 4 samples by name, in entered order (see sort=False above)
        for name, df in grouped:

            # collect metadata
            conc = '{0:.2f}'.format(df['yield'].mean())
            dilution = df['dilution'].mean()
            popt, perr = do_fit(df)
            popt_si, perr_si = do_substrate_inhibition_fit(df)

            # set up plot
            fig, ax = plt.subplots(figsize=(3,3))
            ax.scatter(df.s, df.kobs, color='orange')
            ax.set_title( name )
            ax.set_xlabel( '[pNPG] (mM)' )
            ax.set_ylabel( 'Rate (min$^{-1}$)' )
            ax.set_xticks([ 0, 0.025, 0.050, 0.075 ])
            ax.set_xticklabels([0, 25, 50, 75])
            yticks = ax.get_yticks()
            ax.set_yticks( yticks[1:-1] )
            plt.tight_layout()

            notes = []
            xvals = np.linspace( df.s.min(), df.s.max(), 100 )
            # check if we have evidence of substrate inhibition
            if popt_si[2] < 0.075:
                ax.plot( xvals, kobs_with_substrate_inhibition( xvals, *popt_si ), color='g' )
                notes.append("Likely substrate inhibition (green curve)")
            # check if we have MM
            elif popt[0] > 0.1: # dumb check
                ax.plot( xvals, kobs( xvals, *popt ), color='k' )
                # check on errors
                if perr[0] > (0.5 * popt[0]):
                    notes.append('Note high error for parameter {k_cat}')
                if perr[1] > (0.5 * popt[1]):
                    notes.append('Note high error for parameter {K_M}')
            # neither fit acceptable
            else:
              popt = perr = array([ np.nan, np.nan ])

            # encode plot as a string
            img = StringIO()
            fig.savefig( img, format='svg' )
            plt.close( fig )
            img.seek( 0 )

            # collect everything
            # would be better to decide here if we report 1) kcat and km plot, 2) SI plot, 3) neither
            samples.append((
                name, conc, dilution, popt, perr, popt_si, perr_si, img.read(), notes
            ))
        return render_template( 'kinetics/results.htm', samples=samples )

    else: # not a POST request, is GET request
        data_sets = get_data_sets()
        return render_template('kinetics/input.htm', data_sets=data_sets)

# thermal stability
@app.route( '/thermal', methods=['GET', 'POST'] )
def thermal():
    if request.method == 'POST':

        clean_dat = request.form['data'].replace('Max V [420]', 'rate')

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

        return render_template( 'thermal/results.htm', result=result )
    else:
        return render_template( 'thermal/input.htm' )


@app.route('/oligo_design', methods=['GET', 'POST'])
def oligo_design():
    if request.method == 'POST':
        sequence_text = request.form['sequence_text']
        record = next(SeqIO.parse(StringIO(sequence_text), 'fasta'))
        print(record)

        ecoli_favorite = {
            'G':'GGC', 'A':'GCG', 'V':'GTG', 'F':'TTT', 'E':'GAA',
            'D':'GAT', 'N':'AAC', 'C':'TGC', 'K':'AAA', 'L':'CTG',
            'H':'CAT', 'P':'CCG', 'Q':'CAG', 'W':'TGG', 'Y':'TAT',
            'I':'ATT', 'M':'ATG', 'R':'CGT', 'T':'ACC', 'S':'AGC',
        }

        dna = DNA.read(StringIO(sequence_text))
        kmers = [dna[i:i+33] for i in range(0, len(dna), 3)]

        my_oligos = []
        for i, k in enumerate(kmers):
            for aa, codon in ecoli_favorite.items():
                my_str = str( k[:15] ) + codon + str( k[18:] )
                my_dna = DNA( my_str )
                my_oligo = my_dna.reverse_complement()
                my_name = str( k[15:18].translate() ) + str( i + 6 ) + aa
                if len( my_oligo ) == 33:
                    my_oligos.append('>{}\n{}\n'.format(my_name, my_oligo))
        my_oligos = ''.join(my_oligos)
        # path = os.path.join('/data/bagel/uploads/oligos_{}.fa'.format(datetime.datetime.now()))
        # with open(path, 'w') as fn: fn.write(my_oligos)

        preview_string = my_oligos[0:1000]
        return render_template('oligo_design/output.htm', preview_string=preview_string, fasta_file=my_oligos)
    else: # get
        things_to_render = []
        return render_template('oligo_design/input.htm', things_to_render=things_to_render)

@app.route('/uploads')
def uploads():
    uploads = glob('/data/bagel/uploads/*csv')
    return render_template('uploads.htm', uploads=uploads)
