import sys
sys.stdout = sys.stderr
sys.path.insert( 0, '/data/bagel/' )
from app import app as application
