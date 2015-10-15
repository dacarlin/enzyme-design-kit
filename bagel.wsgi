import sys
sys.stdout = sys.stderr
sys.path.insert( 0, '/data/bagel/' )
from fitter import app as application
