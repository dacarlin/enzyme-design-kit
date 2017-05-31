from skbio import DNA

ecoli_favorite = {
    'G':'GGC', 'A':'GCG', 'V':'GTG', 'F':'TTT', 'E':'GAA', 'D':'GAT', 'N':'AAC', 'C':'TGC', 'K':'AAA', 'L':'CTG',
    'H':'CAT', 'P':'CCG', 'Q':'CAG', 'W':'TGG', 'Y':'TAT', 'I':'ATT', 'M':'ATG', 'R':'CGT', 'T':'ACC', 'S':'AGC',
}

dna = DNA.read( 'bglb_model_coding.fa' ) 
kmers = [ dna[i:i+33] for i in range( 0, len( dna ), 3 ) ]

my_oligos = []
for i, k in enumerate( kmers ):
    for aa, codon in ecoli_favorite.items():
        my_str = str( k[:15] ) + codon + str( k[18:] )
        my_dna = DNA( my_str )
        my_oligo = my_dna.reverse_complement()
        my_name = str( k[15:18].translate() ) + str( i + 6 ) + aa
        if len( my_oligo ) == 33:
            my_oligos.append( '>{}\n{}\n'.format( my_name, my_oligo ) )

with open( 'antisense_oligos.fa', 'w' ) as fn:
    fn.write( ''.join( my_oligos ) )
